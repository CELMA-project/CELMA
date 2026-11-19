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


def commit_with(edits: Dict[str, str]) -> str:
    """Make a commit on top of HEAD with some files replaced.

    Parameters
    ----------
    edits : dict
        Path to the new content.

    Returns
    -------
    str
        Id of the new commit.
    """
    with tempfile.TemporaryDirectory() as tmp:
        env = {**ENV, "GIT_INDEX_FILE": str(Path(tmp) / "index")}
        run_git("read-tree", "HEAD", env=env)
        for path, content in edits.items():
            blob = run_git("hash-object", "-w", "--stdin", data=content.encode())
            run_git(
                "update-index", "--add", "--cacheinfo", f"100644,{blob},{path}", env=env
            )
        tree = run_git("write-tree", env=env)
    return run_git("commit-tree", tree, "-p", "HEAD", "-m", "synthetic")


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


def test_include_and_header_guard_are_ignored_but_code_is_not():
    """Includes and the guard may change, a statement may not."""
    path = "common/BOUTExtensions/include/ownOperators.hxx"
    text = read_head(path)
    new = text.replace("#include", "#include <map>\n#include", 1)
    assert guard.compare(ROOT, "HEAD", commit_with({path: new})) == {}
    assert guard.compare(ROOT, "HEAD", commit_with({path: new + "\nint x;\n"}))


def test_deleted_file_fails_and_added_file_passes():
    """A deleted guarded file trips; a new one does not."""
    added = commit_with({"celma/newFile.cxx": "int x;\n"})
    assert guard.compare(ROOT, "HEAD", added) == {}
    with tempfile.TemporaryDirectory() as tmp:
        env = {**ENV, "GIT_INDEX_FILE": str(Path(tmp) / "index")}
        run_git("read-tree", "HEAD", env=env)
        run_git("update-index", "--force-remove", "celma/celma.cxx", env=env)
        tree = run_git("write-tree", env=env)
    deleted = run_git("commit-tree", tree, "-p", "HEAD", "-m", "synthetic")
    assert guard.compare(ROOT, "HEAD", deleted) == {"celma/celma.cxx": ["file deleted"]}


def test_image_pin_and_compiler_flags_are_guarded():
    """The pinned image and a CMakeLists.txt are compared as they are."""
    workflow = ".github/workflows/test.yml"
    text = read_head(workflow)
    other = text.replace("name: Test", "name: Tests", 1)
    assert guard.compare(ROOT, "HEAD", commit_with({workflow: other})) == {}
    pinned = text.replace("sha256:8ac5", "sha256:9ac5", 1)
    assert list(guard.compare(ROOT, "HEAD", commit_with({workflow: pinned}))) == [
        workflow
    ]
    flags = read_head("CMakeLists.txt") + "\nadd_compile_options(-ffast-math)\n"
    assert "CMakeLists.txt" in guard.compare(
        ROOT, "HEAD", commit_with({"CMakeLists.txt": flags})
    )
