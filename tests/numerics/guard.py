"""Guard that reports every semantic change of the CELMA numerics.

Compares two commits of a git repository and lists the changes in the guarded
content: C++ sources (as token streams, so formatting, comments and header guards do
not count; the set of ``#include`` lines does, their order does not), every
``BOUT.inp`` (as a map of options, so comments, spacing, key order and number spelling do not count) and a set of
files that are compared byte for byte (references, the build image, compiler
flags, every workflow and the guard itself).

Modified, deleted or added guarded content is a change.
The files are read from the git object store only, so nothing of either commit
is checked out or executed.

Run as ``python3 tests/numerics/guard.py BASE HEAD``. The exit code is 1 if
anything changed, else 0.
"""

import argparse
import difflib
import re
import subprocess  # nosec
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

CPP_EXT = r"(cxx|hxx|cpp|hpp|cc|hh|h|c|ixx|ipp|tpp|txx|inl|inc)"
CPP_PATTERNS = tuple(
    re.compile(rf"{folder}/.+\.{CPP_EXT}", re.IGNORECASE)
    for folder in (
        "celma",
        "celmaWithBoussinesqApprox",
        "common/BOUTExtensions",
        "MES",
    )
)
INP_PATTERN = re.compile(r"(.+/)?BOUT\.inp")
RAW_PATTERNS = (
    re.compile(r"tests/regression/[^/]*reference[^/]*"),
    re.compile(r"tests/numerics/.+"),
    re.compile(r"docker/Dockerfile\.bout5"),
    re.compile(r"(.+/)?CMakeLists\.txt"),
    re.compile(r".+\.cmake"),
    # The scripts hold the tolerances of the references
    re.compile(r"tests/regression/.+\.py"),
    re.compile(r"\.github/workflows/.+"),
)
OPERATORS = (
    "<<=",
    ">>=",
    "...",
    "->*",
    "<=>",
    "::",
    "->",
    "++",
    "--",
    "<<",
    ">>",
    "<=",
    ">=",
    "==",
    "!=",
    "&&",
    "||",
    "+=",
    "-=",
    "*=",
    "/=",
    "%=",
    "&=",
    "|=",
    "^=",
    ".*",
    "##",
)
TOKEN = re.compile(
    r'"(?:\\.|[^"\\])*"'
    r"|'(?:\\.|[^'\\])*'"
    r"|\d(?:[eEpP][+-]|[\w.'])*"
    r"|[A-Za-z_]\w*"
    r"|" + "|".join(re.escape(op) for op in OPERATORS) + r"|\S"
)
INCLUDE = re.compile(r"\s*#\s*include\b")
Token = Tuple[str, int]


class GitError(RuntimeError):
    """A git command failed."""


def git(repo: Path, *args: str) -> bytes:
    """Run a git command.

    Parameters
    ----------
    repo : Path
        Repository to run in.
    *args : str
        Arguments of the command.

    Returns
    -------
    bytes
        Standard output.

    Raises
    ------
    GitError
        If git exits with an error.
    """
    done = subprocess.run(  # nosec
        ["git", "-C", str(repo), *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if done.returncode != 0:
        raise GitError(f"git {' '.join(args)}: {done.stderr.decode().strip()}")
    return done.stdout


def list_tree(repo: Path, rev: str) -> Dict[str, str]:
    """List the files of a commit.

    Parameters
    ----------
    repo : Path
        Repository.
    rev : str
        Commit.

    Returns
    -------
    dict
        Path to blob id.
    """
    out = git(repo, "ls-tree", "-r", "-z", rev)
    files = {}
    for entry in out.split(b"\0"):
        if not entry:
            continue
        meta, path = entry.split(b"\t", 1)
        files[path.decode("utf-8", "surrogateescape")] = meta.split()[2].decode()
    return files


SPLICE = re.compile(r"\\[ \t]*\r?\n")


def splice(text: str) -> str:
    r"""Join backslash-newline pairs, as the first phases of C++ translation do.

    Every path that reads C++ (tokens, includes) goes through this first, so a
    directive split over lines (``#\`` newline ``include``) is seen whole.

    Parameters
    ----------
    text : str
        Source code or an input file.

    Returns
    -------
    str
        The text with every continued line joined.
    """
    return SPLICE.sub("", text)


def is_digit_separator(text: str, i: int) -> bool:
    """Tell whether the quote at ``text[i]`` is a C++14 digit separator.

    Parameters
    ----------
    text : str
        Source code.
    i : int
        Position of a single quote.

    Returns
    -------
    bool
        True if the quote sits inside a number such as ``100'000``, not at the
        start of a character literal (``L'a'`` starts with a letter).
    """
    start = i
    while start > 0 and (text[start - 1].isalnum() or text[start - 1] in "_'."):
        start -= 1
    return start < i and text[start].isdigit()


def strip_comments(text: str) -> str:
    """Remove C++ comments.

    Parameters
    ----------
    text : str
        Source code, already spliced.

    Returns
    -------
    str
        The code with every comment replaced by one space, as the compiler
        does, so a directive keeps its line.
    """
    out: List[str] = []
    i, n = 0, len(text)
    while i < n:
        char = text[i]
        two = text[i : i + 2]
        if char == "'" and is_digit_separator(text, i):
            out.append(char)
            i += 1
        elif char in "\"'":
            j = i + 1
            while j < n and text[j] != char and text[j] != "\n":
                j += 2 if text[j] == "\\" else 1
            out.append(text[i : j + 1])
            i = j + 1
        elif two == "//":
            while i < n and text[i] != "\n":
                i += 1
            out.append(" ")
        elif two == "/*":
            end = text.find("*/", i + 2)
            i = n if end < 0 else end + 2
            out.append(" ")
        else:
            out.append(char)
            i += 1
    return "".join(out)


def preprocess(text: str) -> str:
    """Splice lines, strip comments and spell the digraph ``%:`` as ``#``.

    Parameters
    ----------
    text : str
        Source code.

    Returns
    -------
    str
        The text every reader of C++ in this module works on.
    """
    code = strip_comments(splice(text))
    return code.replace("%:%:", "##").replace("%:", "#")


def cpp_tokens(text: str) -> List[Token]:
    """Tokenise C++ without comments, whitespace, includes and header guards.

    Preprocessor lines other than ``#include`` and the header guard stay, each
    closed by an end-of-line token, as their line structure carries meaning.
    Line numbers are those of the spliced text.

    Parameters
    ----------
    text : str
        Source code.

    Returns
    -------
    list
        Tuples of the token and its line number.
    """
    lines = list(enumerate(preprocess(text).split("\n"), 1))
    lines = [x for x in lines if not INCLUDE.match(x[1])]
    lines = drop_header_guard(lines)
    tokens: List[Token] = []
    for number, line in lines:
        found = [(m.group(0), number) for m in TOKEN.finditer(line)]
        tokens.extend(found)
        if found and found[0][0] == "#":
            tokens.append(("<eol>", number))
    return tokens


def cpp_includes(text: str) -> List[str]:
    """List the ``#include`` lines of a source, without order or spacing.

    A new header included from a model changes what the model compiles, so the
    set of includes is guarded, though their order is not.

    Parameters
    ----------
    text : str
        Source code.

    Returns
    -------
    list
        The sorted, whitespace-normalised include lines, without duplicates.
    """
    found = set()
    for line in preprocess(text).split("\n"):
        if INCLUDE.match(line):
            found.add(re.sub(r"\s+", "", line))
    return sorted(found)


def drop_header_guard(lines: List[Tuple[int, str]]) -> List[Tuple[int, str]]:
    """Remove a header guard (``#ifndef NAME``, ``#define NAME``, ``#endif``).

    Parameters
    ----------
    lines : list
        Logical lines with their numbers.

    Returns
    -------
    list
        The lines without the guard, unchanged if there is none.
    """
    real = [i for i, (_, text) in enumerate(lines) if text.strip()]
    if len(real) < 3:
        return lines
    first, second, last = real[0], real[1], real[-1]
    opening = re.fullmatch(r"\s*#\s*ifndef\s+(\w+)\s*", lines[first][1])
    defining = re.fullmatch(r"\s*#\s*define\s+(\w+)\s*", lines[second][1])
    closing = re.fullmatch(r"\s*#\s*endif\b.*", lines[last][1])
    if opening and defining and closing and opening.group(1) == defining.group(1):
        return [x for i, x in enumerate(lines) if i not in (first, second, last)]
    return lines


def snippet(tokens: List[Token]) -> str:
    """Join tokens for a report line.

    Parameters
    ----------
    tokens : list
        Tokens with line numbers.

    Returns
    -------
    str
        The tokens as text, shortened if long.
    """
    text = " ".join(t for t, _ in tokens if t != "<eol>") or "(nothing)"
    return text if len(text) <= 160 else text[:157] + "..."


def compare_cpp(old: str, new: str) -> List[str]:
    """Compare two C++ sources by their token streams.

    Parameters
    ----------
    old : str
        Source at the base commit.
    new : str
        Source at the head commit.

    Returns
    -------
    list
        One line per changed stretch of tokens.
    """
    a, b = cpp_tokens(old), cpp_tokens(new)
    matcher = difflib.SequenceMatcher(
        a=[t for t, _ in a], b=[t for t, _ in b], autojunk=False
    )
    changes = []
    old_inc, new_inc = cpp_includes(old), cpp_includes(new)
    changes.extend(f"include removed: {x}" for x in old_inc if x not in new_inc)
    changes.extend(f"include added: {x}" for x in new_inc if x not in old_inc)
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            continue
        where = a[i1][1] if i1 < len(a) else (a[-1][1] if a else 0)
        where_new = b[j1][1] if j1 < len(b) else (b[-1][1] if b else 0)
        changes.append(
            f"code near line {where} (head line {where_new}): "
            f"{snippet(a[i1:i2])} -> {snippet(b[j1:j2])}"
        )
    return changes


Value = Union[float, str]


def parse_inp(text: str, joined: bool = True) -> Dict[str, List[Value]]:
    """Parse a BOUT.inp into a map of ``section:key`` to its values.

    Comments, spacing and the order of keys are dropped, and a value that
    reads as a number is a float, so ``1.0`` and ``1e0`` are equal. Nothing is
    dropped silently: a key given twice keeps every value in order, and a line
    without ``=`` that is no section header continues the value of the option
    before it (or, with none before it in the section, is kept as a positional
    entry), so a multi-line expression is compared whole.

    Parameters
    ----------
    text : str
        Content of the input file.
    joined : bool
        Join lines that end in a backslash first. Whether a reader does is not
        certain, so the comparison reads the file both ways.

    Returns
    -------
    dict
        Option name to the list of its values, in order of appearance.
    """
    options: Dict[str, List[Value]] = {}
    section = ""
    last = ""
    for raw in (splice(text) if joined else text).split("\n"):
        line = strip_inp_comment(raw).strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip()
            last = ""
            continue
        if "=" not in line:
            if last:
                options[last][-1] = f"{options[last][-1]} {line}"
            else:
                options.setdefault(f"{section}:<positional>", []).append(line)
            continue
        key, value = (x.strip() for x in line.split("=", 1))
        last = f"{section}:{key}" if section else key
        try:
            options.setdefault(last, []).append(float(value))
        except ValueError:
            options.setdefault(last, []).append(value)
    return options


def strip_inp_comment(line: str) -> str:
    r"""Remove a ``#`` comment that is not inside quotes.

    A ``;`` is not taken as a comment start: if a reader does, the text after it
    shows as a changed value, which is the safe way round.

    A backslash inside quotes escapes the next character, so ``\"`` does not
    end the string.

    Parameters
    ----------
    line : str
        One line of an input file.

    Returns
    -------
    str
        The line without its comment.
    """
    quote = ""
    skip = False
    for i, char in enumerate(line):
        if skip:
            skip = False
        elif quote and char == "\\":
            skip = True
        elif quote:
            quote = "" if char == quote else quote
        elif char in "\"'":
            quote = char
        elif char == "#":
            return line[:i]
    return line


def show(value: Optional[List[Value]]) -> str:
    """Format the values of an option for a report line.

    Parameters
    ----------
    value : list or None
        The values, ``None`` for an absent option.

    Returns
    -------
    str
        Text for the report.
    """
    if value is None:
        return "(absent)"
    return str(value[0]) if len(value) == 1 else "[" + ", ".join(map(str, value)) + "]"


def compare_inp(old: str, new: str) -> List[str]:
    """Compare two input files by their options.

    Parameters
    ----------
    old : str
        Content at the base commit.
    new : str
        Content at the head commit.

    Returns
    -------
    list
        One line per added, removed or changed option.
    """
    found: List[str] = []
    for joined in (True, False):
        a, b = parse_inp(old, joined), parse_inp(new, joined)
        for name in sorted(set(a) | set(b)):
            line = f"{name} {show(a.get(name))} -> {show(b.get(name))}"
            if a.get(name) != b.get(name) and line not in found:
                found.append(line)
    return found


def compare_raw(old: str, new: str) -> List[str]:
    """Compare two files byte for byte.

    Parameters
    ----------
    old : str
        Content at the base commit.
    new : str
        Content at the head commit.

    Returns
    -------
    list
        One line if the content differs, else none.
    """
    if old == new:
        return []
    sizes = f"{len(old.encode())} -> {len(new.encode())} bytes"
    diff = difflib.unified_diff(old.splitlines(), new.splitlines(), lineterm="", n=0)
    first = [x for x in diff if x[:1] in "+-" and x[:3] not in ("+++", "---")][:6]
    return [f"content changed ({sizes})"] + ["    " + x[:160] for x in first]


def blob_text(repo: Path, blob: str) -> str:
    """Read a blob as text.

    Parameters
    ----------
    repo : Path
        Repository.
    blob : str
        Blob id.

    Returns
    -------
    str
        Content decoded as UTF-8, undecodable bytes kept as escapes.
    """
    return git(repo, "cat-file", "blob", blob).decode("utf-8", "surrogateescape")


def kind_of(path: str) -> Optional[str]:
    """Classify a path.

    Parameters
    ----------
    path : str
        Path in the repository.

    Returns
    -------
    str or None
        ``"cpp"``, ``"inp"``, ``"raw"`` or ``None`` if unguarded.
    """
    if any(p.fullmatch(path) for p in CPP_PATTERNS):
        return "cpp"
    if INP_PATTERN.fullmatch(path):
        return "inp"
    if any(p.fullmatch(path) for p in RAW_PATTERNS):
        return "raw"
    return None


def compare(repo: Path, base: str, head: str) -> Dict[str, List[str]]:
    """Find every semantic change in the guarded content between two commits.

    Parameters
    ----------
    repo : Path
        Repository holding both commits.
    base : str
        The commit before.
    head : str
        The commit after.

    Returns
    -------
    dict
        Path to the lines describing what changed in it; empty if nothing did.
    """
    old_tree, new_tree = list_tree(repo, base), list_tree(repo, head)
    report: Dict[str, List[str]] = {}
    for path in sorted(set(old_tree) | set(new_tree)):
        kind = kind_of(path)
        if kind is None or old_tree.get(path) == new_tree.get(path):
            continue
        if path not in new_tree:
            changes = ["file deleted"]
        elif path not in old_tree:
            # A new file in guarded scope is a change too: a new local header can
            # shadow an include, a new input or workflow changes what runs
            changes = ["file added"]
        else:
            old = blob_text(repo, old_tree[path])
            new = blob_text(repo, new_tree[path])
            if kind == "cpp":
                changes = compare_cpp(old, new)
            elif kind == "inp":
                changes = compare_inp(old, new)
            else:
                changes = compare_raw(old, new)
        if changes:
            report[path] = changes
    return report


def main(argv: Optional[List[str]] = None) -> int:
    """Print the changes between two commits.

    Parameters
    ----------
    argv : list, optional
        Command line arguments, ``sys.argv[1:]`` by default.

    Returns
    -------
    int
        1 if guarded content changed, 0 if not.
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("base", help="commit before the change")
    parser.add_argument("head", help="commit after the change")
    parser.add_argument("--repo", type=Path, default=Path("."))
    args = parser.parse_args(argv)
    report = compare(args.repo, args.base, args.head)
    for path, changes in report.items():
        print(path)
        for change in changes:
            print(f"  {change}")
    if report:
        print(f"\nNumerics changed in {len(report)} file(s): consent needed.")
        return 1
    print("No change of the numerics.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
