"""Tests of the numerics guard on the real history of the repository.

Needs the full history (``fetch-depth: 0`` in CI). Synthetic changes are made
as commits built with git plumbing on top of ``HEAD``, so the work tree and the
refs are never touched.
"""

import os
import subprocess  # nosec
import tempfile
from pathlib import Path
from typing import Dict

import approval
import guard
import pytest

ROOT = Path(__file__).resolve().parents[2]
ENV = {
    "GIT_AUTHOR_NAME": "t",
    "GIT_AUTHOR_EMAIL": "t@example.com",
    "GIT_COMMITTER_NAME": "t",
    "GIT_COMMITTER_EMAIL": "t@example.com",
}
FORMATTING_ONLY = "b359f167"
PORT = "f812d33e"


def run_git(*args: str, env: Dict[str, str] = ENV, data: bytes = b"") -> str:
    """Run git in the repository and return its output."""
    done = subprocess.run(  # nosec
        ["git", "-C", str(ROOT), *args],
        input=data,
        stdout=subprocess.PIPE,
        check=True,
        env={**os.environ, **env},
    )
    return done.stdout.decode().strip()


def commit_with(edits: Dict[str, str], parent: str = "HEAD") -> str:
    """Make a commit on top of HEAD with some files replaced.

    Parameters
    ----------
    edits : dict
        Path to the new content.
    parent : str
        Commit to build on.

    Returns
    -------
    str
        Id of the new commit.
    """
    with tempfile.TemporaryDirectory() as tmp:
        env = {**ENV, "GIT_INDEX_FILE": str(Path(tmp) / "index")}
        run_git("read-tree", parent, env=env)
        for path, content in edits.items():
            blob = run_git("hash-object", "-w", "--stdin", data=content.encode())
            run_git(
                "update-index", "--add", "--cacheinfo", f"100644,{blob},{path}", env=env
            )
        tree = run_git("write-tree", env=env)
    return run_git("commit-tree", tree, "-p", parent, "-m", "synthetic")


def read_head(path: str) -> str:
    """Read a file at HEAD."""
    return run_git("show", f"HEAD:{path}")


def have(commit: str) -> bool:
    """Tell whether a commit is in the history."""
    return (
        subprocess.run(  # nosec
            ["git", "-C", str(ROOT), "cat-file", "-e", commit + "^{commit}"],
            check=False,
        ).returncode
        == 0
    )


def test_port_commit_fails_and_reports_the_ddz_change(capsys):
    """The commit that moved [ddz] to [mesh:ddz] is a numerics change."""
    assert have(PORT), "the history is shallow: fetch with depth 0"
    report = guard.compare(ROOT, PORT + "^", PORT)
    assert report
    lines = report["celma/CSDXNyScan/BOUT.inp"]
    assert "ddz:first FFT -> (absent)" in " ".join(lines)
    assert "mesh:ddz:first (absent) -> C2" in " ".join(lines)
    assert guard.main(["--repo", str(ROOT), PORT + "^", PORT]) == 1
    assert "celma/CSDXNyScan/BOUT.inp" in capsys.readouterr().out


def test_formatting_only_commit_passes():
    """The clang-format commit changes whitespace and include order only."""
    assert have(FORMATTING_ONLY), "the history is shallow: fetch with depth 0"
    assert guard.compare(ROOT, FORMATTING_ONLY + "^", FORMATTING_ONLY) == {}
    assert (
        guard.main(["--repo", str(ROOT), FORMATTING_ONLY + "^", FORMATTING_ONLY]) == 0
    )


def test_boundary_condition_change_fails():
    """Neumann to Dirichlet at the outer x boundary is reported."""
    path = "celma/CSDXNyScan/BOUT.inp"
    text = read_head(path)
    assert "bndry_xout = neumann_o4(0.0)" in text
    new = text.replace(
        "bndry_xout = neumann_o4(0.0)", "bndry_xout = dirichlet_o4(0.0)", 1
    )
    report = guard.compare(ROOT, "HEAD", commit_with({path: new}))
    assert list(report) == [path]
    assert "bndry_xout" in report[path][0]
    assert "neumann_o4(0.0) -> dirichlet_o4(0.0)" in report[path][0]


def test_stencil_edit_fails():
    """Changing a coefficient of a stencil in ownOperators.cxx is reported."""
    path = "common/BOUTExtensions/src/ownOperators.cxx"
    text = read_head(path)
    old = "4.0 * DDXPhi(ghostIndX - 3, yInd, zInd)"
    assert old in text
    report = guard.compare(
        ROOT,
        "HEAD",
        commit_with(
            {path: text.replace(old, "5.0 * DDXPhi(ghostIndX - 3, yInd, zInd)", 1)}
        ),
    )
    assert list(report) == [path]
    assert "4.0" in report[path][0] and "5.0" in report[path][0]


@pytest.mark.parametrize(
    ("path", "edit"),
    [
        (
            "celma/CSDXNyScan/BOUT.inp",
            lambda t: t.replace("first = FFT", "first   =   FFT # x", 1),
        ),
        ("celma/CSDXNyScan/BOUT.inp", lambda t: "# new comment\n" + t),
        (
            "common/BOUTExtensions/src/ownOperators.cxx",
            lambda t: "// a comment\n" + t.replace("(", " ( ", 5),
        ),
    ],
)
def test_harmless_edits_pass(path, edit):
    """Comments, spacing and spelling of the same value are not changes."""
    text = read_head(path)
    assert guard.compare(ROOT, "HEAD", commit_with({path: edit(text)})) == {}


def test_header_guard_name_is_ignored_but_code_is_not():
    """The name of the header guard may change, a statement may not."""
    path = "common/BOUTExtensions/include/ownOperators.hxx"
    text = read_head(path)
    name = "COMMON_BOUTEXTENSIONS_INCLUDE_OWNOPERATORS_HXX_"
    assert name in text
    renamed = text.replace(name, "OTHER_GUARD_")
    assert guard.compare(ROOT, "HEAD", commit_with({path: renamed})) == {}
    assert guard.compare(ROOT, "HEAD", commit_with({path: text + "\nint x;\n"}))


def test_deleted_and_added_files_in_guarded_scope_fail():
    """A deleted guarded file trips, and so does a new one."""
    added = commit_with({"celma/newFile.cxx": "int x;\n"})
    assert guard.compare(ROOT, "HEAD", added) == {"celma/newFile.cxx": ["file added"]}
    for path in (
        "MES/shadow/ownOperators.hxx",
        "tests/numerics/conftest.py",
        ".github/workflows/extra.yml",
        "celma/NewScan/BOUT.inp",
        "newdir/CMakeLists.txt",
        "tests/regression/new-reference.json",
    ):
        report = guard.compare(ROOT, "HEAD", commit_with({path: "x\n"}))
        assert report == {path: ["file added"]}
    assert guard.compare(ROOT, "HEAD", commit_with({"README.md": "x\n"})) == {}
    with tempfile.TemporaryDirectory() as tmp:
        env = {**ENV, "GIT_INDEX_FILE": str(Path(tmp) / "index")}
        run_git("read-tree", "HEAD", env=env)
        run_git("update-index", "--force-remove", "celma/celma.cxx", env=env)
        tree = run_git("write-tree", env=env)
    deleted = run_git("commit-tree", tree, "-p", "HEAD", "-m", "synthetic")
    assert guard.compare(ROOT, "HEAD", deleted) == {"celma/celma.cxx": ["file deleted"]}


def test_every_workflow_and_compiler_flags_are_guarded():
    """Any edit of a workflow, the pinned image included, or of a CMakeLists."""
    for workflow in (".github/workflows/test.yml", ".github/workflows/lint.yml"):
        text = read_head(workflow)
        assert list(
            guard.compare(ROOT, "HEAD", commit_with({workflow: text + "# x\n"}))
        ) == [workflow]
    assert guard.kind_of(".github/workflows/brand-new.yml") == "raw"
    flags = read_head("CMakeLists.txt") + "\nadd_compile_options(-ffast-math)\n"
    assert "CMakeLists.txt" in guard.compare(
        ROOT, "HEAD", commit_with({"CMakeLists.txt": flags})
    )


@pytest.mark.parametrize(
    "path",
    [
        "celma/CSDXNyScan/deep/model.cxx",
        "celma/sub/a.hpp",
        "celmaWithBoussinesqApprox/x/y/z.h",
        "common/BOUTExtensions/src/new.cpp",
        "common/BOUTExtensions/include/new.cc",
        "MES/a/b/c/case.cxx",
        "MES/a/header.hxx",
    ],
)
def test_cpp_in_subdirectories_and_all_extensions_is_guarded(path):
    """Every C++ extension in any depth of a guarded folder is compared."""
    assert guard.kind_of(path) == "cpp"
    base = commit_with({path: "int x = 1;\n"})
    assert guard.compare(ROOT, base, commit_with({path: "int x = 2;\n"}, base))
    assert not guard.compare(ROOT, base, commit_with({path: "int  x=1; // c\n"}, base))


def test_new_include_is_a_change_but_reordering_is_not():
    """A new header from a model is a change; the order of includes is not."""
    path = "celma/celma.cxx"
    text = read_head(path)
    assert guard.compare(
        ROOT, "HEAD", commit_with({path: '#include "extra/newScheme.hxx"\n' + text})
    )
    lines = text.split("\n")
    includes = [i for i, x in enumerate(lines) if x.startswith("#include")]
    assert len(includes) >= 2
    lines[includes[0]], lines[includes[1]] = lines[includes[1]], lines[includes[0]]
    assert guard.compare(ROOT, "HEAD", commit_with({path: "\n".join(lines)})) == {}


def test_escaped_quote_does_not_end_a_string_in_an_input_file():
    """A comment sign after an escaped quote is still inside the string."""
    line = 'name = "a\\"b # c"  # the comment'
    assert guard.strip_inp_comment(line).strip() == 'name = "a\\"b # c"'
    old = 'name = "a\\"b # c"\n'
    new = 'name = "a\\"b # d"\n'
    assert guard.compare_inp(old, new) == ['name "a\\"b # c" -> "a\\"b # d"']
    assert guard.compare_inp(old, old.rstrip() + " # more\n") == []


ALLOWED = ["loeiten"]
AGENT = "r3becc4-agent"
HEAD = "a" * 40
BASE = "b" * 40


def event(action, sender, label=approval.LABEL):
    """Build a pull request event payload."""
    return {
        "action": action,
        "sender": {"login": sender},
        "label": {"name": label},
        "pull_request": {
            "head": {"sha": HEAD},
            "base": {"sha": BASE, "ref": "main"},
        },
        "statuses": [{"context": "numerics-approval", "state": "success"}],
    }


def live_pr(**changes):
    """Build the pull request as the API shows it now, by default as approved."""
    live = {
        "state": "open",
        "labels": [{"name": approval.LABEL}],
        "head": {"sha": HEAD},
        "base": {"sha": BASE, "ref": "main"},
    }
    live.update(changes)
    return live


def decided(payload, actor="loeiten", live=None):
    """Decide a pull request run."""
    live = live_pr() if live is None else live
    return approval.pull_request_approved(payload, actor, live, ALLOWED)


def test_fresh_label_by_the_approver_passes():
    """The approver applies the label and starts the run: approved."""
    assert decided(event("labeled", "loeiten"))
    assert decided(event("labeled", "Loeiten"), actor="LOEITEN")


@pytest.mark.parametrize(
    ("actor", "live"),
    [
        ("loeiten", live_pr(labels=[])),  # replay after the label was removed
        ("loeiten", live_pr(labels=[{"name": "bug"}])),
        ("loeiten", live_pr(state="closed")),
        ("loeiten", live_pr(head={"sha": "c" * 40})),  # replay after a push
        ("loeiten", live_pr(base={"sha": "c" * 40, "ref": "main"})),  # base moved
        ("loeiten", live_pr(base={"sha": BASE, "ref": "release"})),  # retargeted
        ("loeiten", {}),  # no live data
        (AGENT, live_pr()),  # the agent re-ran the old approval run
        ("", live_pr()),
    ],
)
def test_replay_of_the_approval_event_fails_when_the_pull_request_moved(actor, live):
    """A re-run replays the old event, so the live pull request and actor decide."""
    assert not decided(event("labeled", "loeiten"), actor=actor, live=live)


def test_queued_run_fails_after_an_earlier_run_deleted_the_label():
    """The label was dropped by an earlier run: a queued labeled run fails."""
    dropped = live_pr(labels=[])
    assert not decided(event("labeled", "loeiten"), live=dropped)


@pytest.mark.parametrize(
    "payload",
    [
        event("labeled", AGENT),  # the agent applies the label
        event("labeled", "loeiten", label="bug"),  # another label
        event("synchronize", "loeiten"),  # a push after the label
        event("reopened", "loeiten"),  # close and reopen
        event("edited", "loeiten"),  # retarget
        event("opened", "loeiten"),
        event("unlabeled", "loeiten"),
        event("labeled", ""),
    ],
)
def test_every_other_event_fails_and_a_forged_status_is_not_read(payload):
    """Nothing but the approver applying the label approves."""
    assert not decided(payload)


def run(**changes):
    """Build a workflow run as the API lists it, by default an approving one."""
    title = (
        f"numerics-guard event=pull_request_target action=labeled pr=7 "
        f"head={HEAD} base={BASE} approval=true"
    )
    entry = {
        "event": "pull_request_target",
        "status": "completed",
        "conclusion": "success",
        "run_number": 10,
        "display_title": title,
        "actor": {"login": "loeiten"},
        "triggering_actor": {"login": "loeiten"},
    }
    entry.update(changes)
    return entry


def later(action, number=11, conclusion="failure", **changes):
    """Build a later run of the workflow for the same pull request."""
    title = run()["display_title"].replace("action=labeled", f"action={action}")
    return run(run_number=number, conclusion=conclusion, display_title=title, **changes)


MERGED = "f" * 40
PULLS = [
    {
        "number": 7,
        "head": {"sha": HEAD},
        "base": {"ref": "main"},
        "merge_commit_sha": MERGED,
        "labels": [{"name": approval.LABEL}],
    }
]


def pushed(runs, pulls=PULLS, base=BASE, after=MERGED, branch="main"):
    """Decide a push with the approvers of these tests."""
    return approval.push_approved(
        runs, pulls, ALLOWED, approval.Push(base, after, branch)
    )


def test_push_passes_with_the_approvers_latest_run_for_head_and_base():
    """The approval run is the latest for the merged head and base."""
    assert pushed([run()])
    assert pushed([run(run_number=3, conclusion="failure"), run()])


@pytest.mark.parametrize(
    "runs",
    [
        [run(actor={"login": AGENT}, triggering_actor={"login": AGENT})],
        [run(triggering_actor={"login": AGENT})],  # the agent re-ran it
        [run(actor={"login": AGENT})],
        [run(event="pull_request")],  # a copy of the workflow on the branch
        [run(event="push")],
        [run(conclusion="failure")],
        [run(status="in_progress", conclusion=None)],
        [run(display_title=run()["display_title"].replace(HEAD, "c" * 40))],
        [run(display_title=run()["display_title"].replace(BASE, "c" * 40))],
        [run(display_title=run()["display_title"].replace("pr=7", "pr=8"))],
        [run(display_title=run()["display_title"].replace("labeled", "synchronize"))],
        [run(display_title=run()["display_title"] + " head=" + HEAD)],
        [run(display_title="numerics-guard approval=true")],
        # a crafted label name cannot reach the title: only true or false does
        [run(display_title=run()["display_title"] + " label=" + approval.LABEL)],
        [run(display_title=run()["display_title"].replace("true", "numerics"))],
        [],
        # events after the approval void it, whatever their outcome
        [run(), later("unlabeled")],
        [run(), later("synchronize")],
        [run(), later("edited")],
        [run(), later("reopened", conclusion="success")],
    ],
)
def test_push_fails_without_the_approvers_own_latest_run(runs):
    """A run by the agent, of another event, head, base or label, or a later run."""
    assert not pushed(runs)


def test_push_fails_if_head_base_or_label_differ_at_the_merge():
    """A push after the approval changed the head; main moved; the label was off."""
    other_head = [
        {"number": 7, "head": {"sha": "d" * 40}, "labels": PULLS[0]["labels"]}
    ]
    assert not pushed([run()], other_head)
    assert not pushed([run()], [])
    assert not pushed([run()], base="e" * 40)
    no_label = [{"number": 7, "head": {"sha": HEAD}, "labels": []}]
    assert not pushed([run()], no_label)
    assert not pushed([run()], [{"number": 7, "head": {"sha": HEAD}}])


def test_push_fails_for_a_pull_request_that_was_not_merged_into_this_push():
    """The pull request must be merged into the pushed branch, as the pushed commit."""
    elsewhere = [{**PULLS[0], "base": {"ref": "release"}}]
    assert not pushed([run()], elsewhere)
    assert not pushed([run()], branch="release")
    other_commit = [{**PULLS[0], "merge_commit_sha": "9" * 40}]
    assert not pushed([run()], other_commit)
    assert not pushed([run()], after="9" * 40)
    assert not pushed([run()], after="")
    assert not pushed([run()], [{k: v for k, v in PULLS[0].items() if k != "base"}])


def test_run_name_has_no_free_text_field():
    """The title is event, action, numbers, shas and a boolean, nothing a user types."""
    text = (ROOT / ".github/workflows/numerics-guard.yml").read_text()
    assert "github.event.label.name }}" not in text
    assert (
        "approval=${{ github.event.label.name == 'numerics-change-approved' }}" in text
    )


def test_the_scripts_of_the_push_path_come_from_the_base():
    """The workflow extracts guard.py and approval.py from git and runs them -I."""
    text = (ROOT / ".github/workflows/numerics-guard.yml").read_text()
    assert 'git show "$TRUSTED:tests/numerics/$script"' in text
    assert text.count("python3 ") == text.count("python3 -I")
    assert "python3 tests/" not in text
    assert "cd tests/numerics" not in text


def test_line_splicing_cannot_hide_an_include_or_a_token():
    """A directive split by a backslash is seen whole, a splice alone is no change."""
    path = "celma/celma.cxx"
    text = read_head(path)
    for hidden in ('#\\\ninclude "evil.hxx"\n', '# \\ \n include "evil.hxx"\n'):
        report = guard.compare(ROOT, "HEAD", commit_with({path: hidden + text}))
        assert any("include added" in x for x in report[path]), hidden
    assert guard.compare_cpp("int x = 1;\n", "int x = \\\n1;\n") == []
    assert guard.compare_cpp("int x = 1;\n", "int x = \\\n2;\n")
    assert guard.compare_cpp("#define A 1\n", "#def\\\nine A 2\n")
    assert guard.compare_cpp("#include <a>\n", "%:include <b>\n")


def test_input_file_continuations_are_not_dropped():
    """A continued line is read both ways, so a hidden option shows."""
    old = "x = 1\ny = 5\n"
    new = "x = 1 \\\ny = 5\n"
    assert guard.compare_inp(old, new)
    assert guard.compare_inp(old, old) == []
    assert guard.compare_inp("a = 1\n", "a = 1 \\\n\n")
    path = "celma/CSDXNyScan/BOUT.inp"
    edited = read_head(path).replace("nx = 32", "nx = 32 \\\n+ 0", 1)
    assert guard.compare(ROOT, "HEAD", commit_with({path: edited}))


def test_digit_separators_do_not_open_character_literals():
    """A comment after 100'000 is a comment, a prefixed char literal is one."""
    plain = "int n = 100'000;\nint k = 1;\n"
    commented = "int n = 100'000; // it's a comment\nint k = 1; /* don't */\n"
    assert guard.compare_cpp(plain, commented) == []
    assert guard.compare_cpp(plain, plain.replace("k = 1", "k = 2"))
    assert guard.compare_cpp("char c = L'a';\n", "char c = L'a'; // it's\n") == []
    assert guard.compare_cpp("char c = L'a';\n", "char c = L'b';\n")
    assert guard.compare_cpp("int n = 0x1'F;\n", "int n = 0x1'E;\n")


def test_input_lines_without_equals_sign_and_duplicate_keys_count():
    """A multi-line expression and a repeated key are compared whole."""
    base = "[a]\nx = (1 +\n  2)\ny = 1\n"
    assert guard.compare_inp(base, base) == []
    assert guard.compare_inp(base, base.replace("  2)", "  3)"))  # continuation
    assert guard.compare_inp(base, base.replace("  2)", ""))  # line removed
    # a positional line before any option of the section
    assert guard.compare_inp("[a]\ny = 1\n", "[a]\nstray\ny = 1\n")
    # a duplicate key added, changed, removed: the last value no longer hides
    assert guard.compare_inp("[a]\ny = 1\n", "[a]\ny = 1\ny = 1\n")
    assert guard.compare_inp("[a]\ny = 1\ny = 2\n", "[a]\ny = 1\ny = 3\n")
    assert guard.compare_inp("[a]\ny = 1\ny = 2\n", "[a]\ny = 2\n")
    assert guard.compare_inp("[a]\ny = 9\ny = 2\n", "[a]\ny = 2\n")
    assert guard.compare_inp("[a]\ny = 1\n", "[a]\ny = 1\n[a]\ny = 1\n")
    assert guard.parse_inp("[a]\ny = 1\ny = 2\n") == {"a:y": [1.0, 2.0]}
    path = "celma/CSDXNyScan/BOUT.inp"
    text = read_head(path)
    doubled = text.replace("first = FFT", "first = FFT\nfirst = C2", 1)
    assert guard.compare(ROOT, "HEAD", commit_with({path: doubled}))


def test_semicolon_is_not_a_comment_in_an_input_file():
    """Text after a semicolon is part of the value, so it cannot hide a change."""
    assert guard.compare_inp("y = 1\n", "y = 1 ; z = 2\n")
    assert guard.compare_inp("y = 1 # c\n", "y = 1 # d\n") == []


@pytest.mark.parametrize(
    "value",
    ["1_000", " inf", "inf", "nan", "1e3 x", "0x10", "1,5", "١٢٣", "1e", "--1"],
)
def test_only_plain_decimals_are_numbers_in_an_input_file(value):
    """float() accepts spellings a reader may not, so those stay text and trip."""
    assert guard.compare_inp("n = 1000\n", f"n = {value}\n")
    assert guard.parse_inp(f"n = {value}\n")["n"] == [value.strip()]


def test_equal_decimal_spellings_are_the_same_number():
    """Decimal and exponent spellings of one number are no change."""
    for other in ("1000", "1000.0", "1e3", "1.0E+3", "+1000", "10e2"):
        assert guard.compare_inp("n = 1000\n", f"n = {other}\n") == []


def test_cpp_is_guarded_anywhere_and_mode_changes_trip():
    """C++ in any folder is guarded; a mode or symlink change on a guarded path trips."""
    for path in ("tools/x.cxx", "tests/radialLowPass/a.hxx", "top.cpp", "a/b/c/d.H"):
        assert guard.kind_of(path) == "cpp", path
    assert guard.kind_of("tools/readme.md") is None
    path = "celma/celma.cxx"
    blob = run_git("rev-parse", f"HEAD:{path}")
    with tempfile.TemporaryDirectory() as tmp:
        env = {**ENV, "GIT_INDEX_FILE": str(Path(tmp) / "index")}
        run_git("read-tree", "HEAD", env=env)
        run_git("update-index", "--cacheinfo", f"100755,{blob},{path}", env=env)
        executable = run_git(
            "commit-tree", run_git("write-tree", env=env), "-p", "HEAD", "-m", "mode"
        )
        link = run_git("hash-object", "-w", "--stdin", data=b"../elsewhere.cxx")
        run_git("update-index", "--cacheinfo", f"120000,{link},{path}", env=env)
        symlink = run_git(
            "commit-tree", run_git("write-tree", env=env), "-p", "HEAD", "-m", "link"
        )
    mode = guard.compare(ROOT, "HEAD", executable)
    assert mode == {path: ["file mode 100644 -> 100755"]}
    assert any("120000" in x for x in guard.compare(ROOT, "HEAD", symlink)[path])
