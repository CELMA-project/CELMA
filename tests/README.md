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
    * `fingerprint.py` runs both models for exactly one right-hand-side
      evaluation (explicit Euler, one step) from z-, x- and y-dependent initial
      fields and compares the time derivative of every evolved field and every
      saved term, plus a few operator cases of `MES/`, with
      `fingerprint-reference.json` at relative tolerance 1e-9; it names the term
      that moved.
      It reruns with a changed z derivative, x derivative, upwind scheme and
      boundary condition and requires each to fail by a wide margin, and runs
      the unchanged case twice to show the numbers are reproducible.
      CI runs it in the image, which also records the reference
      (`recorded-fingerprint.json`).
    * `guard.py BASE HEAD` compares two commits and lists every semantic change
      in the guarded content: C++ as token streams (formatting, comments,
      includes and header guards ignored), every `BOUT.inp` as a map of options
      (comments, spacing, order and number spelling ignored), and byte for byte
      the references, `docker/Dockerfile.bout5`, the image line of `test.yml`,
      every `CMakeLists.txt` and the guard's own files. Modified or deleted
      guarded content is a change, a new file is not. It exits 1 on a change.
    * `test_guard.py` proves the guard on the history: it fails on the commit
      that moved `[ddz]` to `[mesh:ddz]`, passes on the clang-format commit and
      fails on a synthetic boundary condition and stencil change.
      The `Lint` workflow runs it with the full history
      (`python3 -m pytest tests/numerics/test_guard.py`).

  The workflow [numerics-guard.yml](../.github/workflows/numerics-guard.yml)
  (check `numerics-guard`) runs the guard on every pull request, with the
  base branch's copy of `guard.py`, and never runs code of the pull request.
  A change passes only if the pull request carries the label
  `numerics-change-approved` and `loeiten` applied it last; a push to the pull
  request removes the label. A push to `main` is checked against its first
  parent and passes only if the merged pull request had the approved label.
  The label is Michael's consent for that specific change, see
  [AGENTS.md](../AGENTS.md#numerics-need-michaels-consent).
