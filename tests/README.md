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
