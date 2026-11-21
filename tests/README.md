# tests

[Root](../README.md)

Tests the continuous integration runs after the
[CMake build](../README.md#build-and-test), from
[test.yml](../.github/workflows/test.yml).
Operators and boundary conditions are checked by the cases in
[../MES](../MES/README.md), which the second regression test below also
runs.

* [radialLowPass](radialLowPass/conclusion.txt) - Checks that the radial low
  pass filter works as expected: up to round-off error it does.
  CI runs `./radialLowPass -d data` in this directory.
  `driverTest.py` and `postProcess/checkModes.py` are the older driver and the
  mode check; the driver needs `bout_runners`, which CI does not install.
* [regression](regression) - Two scripts, run from the repository root:
  * `regression.py` runs both models on a small mesh and compares the time
      derivatives and the evolved fields with `reference.json.gz`.
  * `mes.py` runs every `BOUT.inp` under `MES/` at two small resolutions and
      compares the largest absolute error with `mes-reference.json`.

  The references depend on the `BOUT++` version and the compiler, so CI
  records them inside its image.
  Each script writes the numbers of the run to `recorded-*.json` in the
  current directory.
* [numerics](numerics) - Guards against an unintended change of the numerics:
  * `fingerprint.py` runs both models for exactly one right-hand-side evaluation
    (`nout=0`: one evaluation, no solver step, asserted in every run) from z-,
    x- and y-dependent initial fields and compares the time derivative of
    every evolved field and every saved term, plus a few operator cases of
    `MES/`, with `fingerprint-reference.json` at relative tolerance 1e-9; it
    names the term that moved. It reruns with a changed z derivative, x
    derivative, upwind scheme and boundary condition and requires each to fail
    by a wide margin, and runs the unchanged case twice and requires every
    value bit for bit the same. CI runs it in the image, which also records
    the reference (`recorded-fingerprint.json`).
  * `guard.py BASE HEAD` compares two commits and lists every semantic change in
    the guarded content: C++ (any folder, any C++ extension) as token streams
    after joining continued lines, with raw string literals as one token
    (formatting and comments ignored; `#include` lines may be reordered only
    among lines directly after one another; a header guard is a normal
    change), every `BOUT.inp` as a map of options (comments, spacing, order
    and number spelling ignored; continued lines read both ways; duplicate
    keys and lines without `=` all count), and byte for byte the references,
    the scripts of `tests/regression`, `docker/Dockerfile.bout5`, every
    `CMakeLists.txt` and `.cmake`, everything under `.github/` and
    `tests/numerics/`, `.pre-commit-config.yaml`, every `pyproject.toml`,
    `setup.cfg`, `pytest.ini`, `tox.ini` and `conftest.py`, and any change of
    the file mode or a symlink on guarded paths. Modified, deleted or added
    guarded content is a change (a new local header can shadow an include). It
    exits 1 on a change.
  * `test_guard.py` proves the guard on the history: it fails on the commit that
    moved `[ddz]` to `[mesh:ddz]`, passes on a clang-format commit (and on the
    C++ of the first one) and fails on a synthetic boundary condition and
    stencil change. The `Lint` workflow runs it with the full history
    (`python3 -m pytest tests/numerics/test_guard.py`).

  The workflow [numerics-guard.yml](../.github/workflows/numerics-guard.yml)
  (check `numerics-guard`) runs the guard on every pull request and every push
  to `main`, with the scripts of the base commit (extracted from git, run with
  `python3 -I`), and never runs code of the pull request or of the pushed
  commit. Nothing is stored. A numerics change passes only in the run for the
  event "an approver applied the label `numerics-change-approved`" (GitHub sets
  the sender), which approves the head and base that event carries, provided the
  run was started by an approver and the live pull request is still open,
  labelled and at that head and base (a re-run replays the old event). Every
  other event with a change (push, reopen, retarget, another label or sender)
  fails and drops the label. A push to `main` passes only if the latest
  `pull_request_target` run for the merged pull request is a successful one of
  that kind, started and re-run by an approver, for the merged head and the base
  the push started from, for a pull request merged into `main` as exactly the
  pushed commit (`merge_commit_sha`), which still had the label; a copy of the
  workflow on a branch cannot make one, and any later event voids it, except
  another label being added or removed: that run only re-reads whether the
  approval holds and is not counted as a later event. The base of the approval
  counts as the same when the fork point with the head is the same or a
  descendant of the approved one (a branch that moved on).
  `approval.py` holds the decisions and `test_guard.py` tests them. Every file
  under `.github/workflows/` is guarded too. The label is a maintainer's consent
  for that specific change, see
  [AGENTS.md](../AGENTS.md#numerics-need-a-maintainers-consent).
