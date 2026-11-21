"""Tests of the numerics guard on the real history of the repository.

Needs the full history (``fetch-depth: 0`` in CI). Synthetic changes are made
as commits built with git plumbing on top of ``HEAD``, so the work tree and the
refs are never touched.
"""

import os
import subprocess  # nosec
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

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
FORMATTING_ONLY_2 = "4b1ad55f"
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


def test_formatting_only_commits_pass():
    """Whitespace, comments and include order are no change of the C++."""
    assert have(FORMATTING_ONLY), "the history is shallow: fetch with depth 0"
    assert have(FORMATTING_ONLY_2), "the history is shallow: fetch with depth 0"
    # The clang-format commit also added a hook to .pre-commit-config.yaml, which
    # is guarded now; every C++ file in it is formatting only
    report = guard.compare(ROOT, FORMATTING_ONLY + "^", FORMATTING_ONLY)
    assert list(report) == [".pre-commit-config.yaml"]
    # A commit that only reformats one C++ file passes as a whole
    assert guard.compare(ROOT, FORMATTING_ONLY_2 + "^", FORMATTING_ONLY_2) == {}
    args = ["--repo", str(ROOT), FORMATTING_ONLY_2 + "^", FORMATTING_ONLY_2]
    assert guard.main(args) == 0
    assert (
        guard.main(["--repo", str(ROOT), FORMATTING_ONLY + "^", FORMATTING_ONLY]) == 1
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


ALLOWED = ["approver"]
AGENT = "r3becc4-agent"
HEAD = "a" * 40
BASE = "b" * 40


def event(action: str, sender: str, label: str = approval.LABEL) -> Dict[str, Any]:
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


def live_pr(**changes: Any) -> Dict[str, Any]:
    """Build the pull request as the API shows it now, by default as approved."""
    live = {
        "state": "open",
        "labels": [{"name": approval.LABEL}],
        "head": {"sha": HEAD},
        "base": {"sha": BASE, "ref": "main"},
    }
    live.update(changes)
    return live


def decided(
    payload: Dict[str, Any],
    actor: str = "approver",
    live: Optional[Dict[str, Any]] = None,
) -> bool:
    """Decide a pull request run."""
    live = live_pr() if live is None else live
    return approval.pull_request_approved(payload, actor, live, ALLOWED)


def test_fresh_label_by_the_approver_passes():
    """The approver applies the label and starts the run: approved."""
    assert decided(event("labeled", "approver"))
    assert decided(event("labeled", "Approver"), actor="APPROVER")


@pytest.mark.parametrize(
    ("actor", "live"),
    [
        ("approver", live_pr(labels=[])),  # replay after the label was removed
        ("approver", live_pr(labels=[{"name": "bug"}])),
        ("approver", live_pr(state="closed")),
        ("approver", live_pr(head={"sha": "c" * 40})),  # replay after a push
        ("approver", live_pr(base={"sha": "c" * 40, "ref": "main"})),  # base moved
        ("approver", live_pr(base={"sha": BASE, "ref": "release"})),  # retargeted
        ("approver", {}),  # no live data
        (AGENT, live_pr()),  # the agent re-ran the old approval run
        ("", live_pr()),
    ],
)
def test_replay_of_the_approval_event_fails_when_the_pull_request_moved(actor, live):
    """A re-run replays the old event, so the live pull request and actor decide."""
    assert not decided(event("labeled", "approver"), actor=actor, live=live)


def test_queued_run_fails_after_an_earlier_run_deleted_the_label():
    """The label was dropped by an earlier run: a queued labeled run fails."""
    dropped = live_pr(labels=[])
    assert not decided(event("labeled", "approver"), live=dropped)


@pytest.mark.parametrize(
    "payload",
    [
        event("labeled", AGENT),  # the agent applies the label
        event("labeled", "approver", label="bug"),  # another label
        event("synchronize", "approver"),  # a push after the label
        event("reopened", "approver"),  # close and reopen
        event("edited", "approver"),  # retarget
        event("opened", "approver"),
        event("unlabeled", "approver"),
        event("labeled", ""),
    ],
)
def test_every_other_event_fails_and_a_forged_status_is_not_read(payload):
    """Nothing but the approver applying the label approves."""
    assert not decided(payload)


def run(**changes: Any) -> Dict[str, Any]:
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
        "actor": {"login": "approver"},
        "triggering_actor": {"login": "approver"},
    }
    entry.update(changes)
    return entry


def later(
    action: str, number: int = 11, conclusion: str = "failure", **changes: Any
) -> Dict[str, Any]:
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
        assert report[path], hidden
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


def test_rebase_and_merge_of_a_multi_commit_pull_request():
    """BEFORE..AFTER spans all rebased commits; the approval is for that base."""
    first = commit_with(
        {
            "celma/CSDXNyScan/BOUT.inp": read_head("celma/CSDXNyScan/BOUT.inp")
            + "\nextra = 1\n"
        }
    )
    second = commit_with({"celma/extra.cxx": "int x;\n"}, first)
    third = commit_with({"README.md": "x\n"}, second)
    report = guard.compare(ROOT, "HEAD", third)
    assert set(report) == {"celma/CSDXNyScan/BOUT.inp", "celma/extra.cxx"}
    before, after = BASE, MERGED
    pull = {**PULLS[0], "merge_commit_sha": after}
    assert pushed([run()], [pull], base=before, after=after)
    # main moved between the approval and the merge: the base differs
    assert not pushed([run()], [pull], base="e" * 40, after=after)
    # the rebase made the last commit the merge commit, not an earlier one
    assert not pushed([run()], [pull], base=before, after="1" * 40)


def test_the_push_path_reads_every_page_of_associated_pull_requests():
    """A commit with many associated pull requests is read whole."""
    text = (ROOT / ".github/workflows/numerics-guard.yml").read_text()
    assert 'gh api --paginate "repos/$REPO/commits/$AFTER/pulls' in text


@pytest.mark.parametrize(
    "path",
    [
        ".pre-commit-config.yaml",
        "pyproject.toml",
        "tests/sub/pyproject.toml",
        "setup.cfg",
        "pytest.ini",
        "tox.ini",
        "conftest.py",
        "tests/regression/conftest.py",
        ".github/dependabot.yml",
        ".github/actions/local/action.yml",
        ".github/CODEOWNERS",
    ],
)
def test_what_lint_and_pytest_run_is_guarded(path):
    """Hooks, tool settings, conftest files and everything under .github trip."""
    assert guard.kind_of(path) == "raw"
    report = guard.compare(ROOT, "HEAD", commit_with({path: "x\n"}))
    assert list(report) == [path]  # added, or the content changed if it exists
    assert guard.kind_of(".clang-format") is None


def test_raw_string_literals_are_one_token():
    """Comment signs and quotes inside R"(...)" do not start comments or strings."""
    plain = 'const char *s = R"(")"; int k = 1;\n'
    assert guard.compare_cpp(plain, plain.rstrip() + " // c\n") == []
    assert guard.compare_cpp(plain, plain.replace("k = 1", "k = 2"))
    slash = 'auto s = R"x(a // b /* c)x"; int k = 1; /* gone */\n'
    assert guard.compare_cpp(slash, slash.replace(" /* gone */", "")) == []
    assert guard.compare_cpp(slash, slash.replace("a //", "b //"))
    assert guard.compare_cpp(slash, slash.replace("k = 1", "k = 2"))
    for prefix in ("u8R", "uR", "UR", "LR"):
        text = f'auto s = {prefix}"(a /* b)"; int k = 1; // c\n'
        assert guard.compare_cpp(text, text.replace("k = 1", "k = 2"))
        assert guard.compare_cpp(text, text.replace(" // c", "")) == []
    # a multi-line literal with an include-looking line is data, not an include
    multi = 'auto s = R"(\n#include <x>\n)";\nint k = 1;\n'
    assert not any(t == "include" for t, _ in guard.cpp_tokens(multi))
    assert guard.compare_cpp(multi, multi.replace("<x>", "<y>"))
    # an identifier ending in R is not a prefix
    assert guard.cpp_tokens('x = fooR"(a)";')[2] == ("fooR", 1)


def test_include_blocks_sort_but_an_include_cannot_move_past_code():
    """clang-format sorts adjacent includes; a move past a macro changes meaning."""
    block = "#include <b>\n#include <a>\nint k;\n"
    assert guard.compare_cpp(block, "#include <a>\n#include <b>\nint k;\n") == []
    moved = "#include <b>\n#define X 1\n#include <a>\nint k;\n"
    swapped = "#include <b>\n#include <a>\n#define X 1\nint k;\n"
    assert guard.compare_cpp(moved, swapped)
    # a blank line splits blocks, so an include moved across it is a change
    split = "#include <a>\n\n#include <b>\n"
    assert guard.compare_cpp(split, "#include <b>\n\n#include <a>\n")
    assert guard.compare_cpp("int k;\n", "#include <new.hxx>\nint k;\n")


def test_header_guards_are_not_ignored():
    """A guard whose name collides with a defined macro switches a file off."""
    body = "int k;\n"
    guarded = "#ifndef A_H\n#define A_H\nint k;\n#endif\n"
    assert guard.compare_cpp(body, guarded)
    assert guard.compare_cpp(guarded, guarded.replace("A_H", "STRING_H"))
    assert guard.compare_cpp(guarded, body)
    assert guard.compare_cpp(guarded, guarded.replace("#endif", "#endif  // A_H")) == []


def test_numbers_that_overflow_stay_text():
    """1e400 and 1e999 are both inf as floats, so they must not compare equal."""
    assert guard.compare_inp("n = 1e400\n", "n = 1e999\n")
    assert guard.compare_inp("n = 1e400\n", "n = 1e400\n") == []


def test_another_label_event_is_not_an_event_of_the_approval():
    """Labels other than the approval label are recognised by their payload."""
    assert approval.is_other_label_event(event("labeled", "someone", label="bug"))
    assert approval.is_other_label_event(event("unlabeled", "someone", label="bug"))
    assert not approval.is_other_label_event(event("labeled", "approver"))
    assert not approval.is_other_label_event(event("unlabeled", "approver"))
    assert not approval.is_other_label_event(event("synchronize", "approver"))


def test_runs_for_another_label_do_not_void_the_push_approval():
    """A label other than the approval label leaves the approval the latest run."""
    other = {"approval": "false"}
    for action in ("labeled", "unlabeled"):
        title = run()["display_title"].replace("action=labeled", f"action={action}")
        title = title.replace("approval=true", "approval=false")
        newer = run(run_number=11, conclusion="failure", display_title=title, **other)
        assert pushed([run(), newer])
    assert not pushed([run(), later("unlabeled")])


def live_now(**changes: Any) -> Dict[str, Any]:
    """Build the live pull request with its number, by default approved."""
    return live_pr(number=7, **changes)


def test_read_of_the_live_pull_request_keeps_the_approval():
    """A run for another label only re-reads: it passes while the approval holds."""
    assert approval.pull_request_read([run()], live_now(), ALLOWED)
    assert not approval.pull_request_read([], live_now(), ALLOWED)
    assert not approval.pull_request_read([run()], live_now(labels=[]), ALLOWED)
    assert not approval.pull_request_read([run()], live_now(state="closed"), ALLOWED)
    moved = live_now(head={"sha": "c" * 40})
    assert not approval.pull_request_read([run()], moved, ALLOWED)
    assert not approval.pull_request_read(
        [run(), later("synchronize")], live_now(), ALLOWED
    )


def test_a_base_that_lags_main_is_compared_by_the_change_it_makes():
    """The approval holds when the base sha differs but the change is the same."""
    seen = []

    def covers(approved: str, now: str, head: str) -> bool:
        seen.append((approved, now, head))
        return True

    lagging = "d" * 40
    assert not pushed([run()], base=lagging)
    assert approval.push_approved(
        [run()], PULLS, ALLOWED, approval.Push(lagging, MERGED, "main"), covers
    )
    assert seen == [(BASE, lagging, HEAD)]
    assert not approval.push_approved(
        [run()],
        PULLS,
        ALLOWED,
        approval.Push(lagging, MERGED, "main"),
        lambda *_: False,
    )
    other_head = [{**PULLS[0], "head": {"sha": "c" * 40}}]
    assert not approval.push_approved(
        [run()], other_head, ALLOWED, approval.Push(BASE, MERGED, "main"), covers
    )


def test_git_covers_follows_the_fork_point_of_the_head():
    """The same change against a newer main is covered; a grown change is not."""
    work = ROOT
    older = run_git("rev-parse", "HEAD~3")
    newer = run_git("rev-parse", "HEAD")
    side = commit_with({"README.md": "side\n"}, older)
    head = commit_with({"README.md": "head\n"}, side)
    assert approval.git_covers(side, side, head)
    assert approval.git_covers(older, older, head)
    assert approval.git_covers(older, side, head) is True
    assert not approval.git_covers(side, older, head)
    assert not approval.git_covers(older, newer, "0" * 40)
    assert work


@pytest.mark.parametrize(
    "number", [".5e-3", ".5", "1.5e+3", "0x1.8p-2", ".5f", "1'000.5e-3"]
)
def test_a_float_literal_is_one_token(number):
    """A literal that starts with the point is not split into two tokens."""
    tokens = [t for t, _ in guard.cpp_tokens(f"x = {number};")]
    assert tokens == ["x", "=", number, ";"]
    assert guard.compare_cpp(f"x = {number};", f"x = {number}1;")
    assert guard.cpp_tokens("a.b;")[1][0] == "."
    assert guard.cpp_tokens("a .* b;")[1][0] == ".*"


def test_a_raw_string_delimiter_cannot_hold_a_carriage_return():
    """R"a<CR>b( is not a raw string literal, as in the standard's list."""
    assert not guard.raw_string_start('R"a\rb(x)a\rb"', 1)
    assert guard.raw_string_start('R"ab(x)ab"', 1)
    tokens = [t for t, _ in guard.cpp_tokens('f(R"a\rb(x)a\rb");')]
    assert not any(t.startswith('R"a') and t.endswith('b"') for t in tokens)


def test_digraphs_are_not_replaced_inside_literals():
    """%: is # in code only, never inside a string or a character literal."""
    assert guard.compare_cpp('s = "%:";\n', 's = "#";\n')
    assert guard.compare_cpp('s = "%:%:";\n', 's = "##";\n')
    assert guard.compare_cpp('s = R"(%:)";\n', 's = R"(#)";\n')
    assert guard.compare_cpp("c = '%:';\n", "c = '#';\n")
    assert guard.compare_cpp("%:define A 1\n", "#define A 1\n") == []
    assert guard.compare_cpp("%:%:\n", "##\n") == []


def test_the_workflow_does_not_delete_the_label_for_another_label():
    """The label is dropped only on the path of an event that is not another label."""
    text = (ROOT / ".github/workflows/numerics-guard.yml").read_text()
    other = text.index('approval.py" other-label "$GITHUB_EVENT_PATH"')
    read = text.index('approval.py" read')
    delete = text.index("-X DELETE")
    assert other < read < delete
    assert text[read:delete].count("exit 1") == 1


def test_the_push_check_reads_the_runs_since_the_pull_request_was_made():
    """The runs query is bounded by a date, not the whole history."""
    text = (ROOT / ".github/workflows/numerics-guard.yml").read_text()
    assert text.count('"created=>=') == 2
    assert "runs?event=pull_request_target&per_page=100" not in text


def test_a_lagging_base_does_not_fail_the_live_check():
    """live.base.sha moved on after the event: the change is the same."""
    moved = live_pr(base={"sha": "d" * 40, "ref": "main"})
    assert not decided(event("labeled", "approver"), live=moved)
    assert approval.pull_request_approved(
        event("labeled", "approver"),
        "approver",
        moved,
        ALLOWED,
        lambda *_: True,
    )
    assert not approval.pull_request_approved(
        event("labeled", "approver"),
        "approver",
        live_pr(base={"sha": "d" * 40, "ref": "release"}),
        ALLOWED,
        lambda *_: True,
    )


def test_a_digit_separator_in_a_float_that_starts_with_the_point():
    """.5'000 is one number, so its quote does not open a character literal."""
    assert guard.is_digit_separator(".5'000", 2)
    tokens = [t for t, _ in guard.cpp_tokens("x = .5'000; // it's\ny = 1;")]
    assert tokens == ["x", "=", ".5'000", ";", "y", "=", "1", ";"]
