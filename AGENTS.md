# AGENTS.md

Entry point for coding agents. Humans: see [README.md](README.md).

CELMA is a `BOUT++` 5.2.0 model of plasma turbulence in a cylinder, plus the
cases that verify its operators.

## Layout

| Directory | What it is | Index |
| --- | --- | --- |
| `celma/` | The CELMA model (`celma.cxx`, `celma.hxx`), scan inputs, legacy run scripts | [README](celma/README.md) |
| `celmaWithBoussinesqApprox/` | The model with the Boussinesq approximation (`celmaWBA.cxx`) | [README](celmaWithBoussinesqApprox/README.md) |
| `common/` | `BOUTExtensions/`: own operators, boundary conditions, filters, monitors, Laplace inversion, parameters | [README](common/README.md) |
| `MES/` | Verification cases (method of manufactured solutions) | [README](MES/README.md) |
| `derivations/` | Notebooks and scripts that derive stencils and boundary conditions | [README](derivations/README.md) |
| `tests/` | Filter test, the regression tests and the numerics guard CI runs | [README](tests/README.md) |
| `docker/` | `Dockerfile.bout5`: the build image of CI. `Dockerfile` is old and unused | [Dockerfile.bout5](docker/Dockerfile.bout5) |
| `.github/` | Workflows `lint.yml`, `test.yml`, `image.yml` | [workflows](.github/workflows) |
| `.codeCleaner/` | Old helper scripts, not used by CI | [README](.codeCleaner/README.md) |

Every model and case has its own `CMakeLists.txt`, listed in the root
[CMakeLists.txt](CMakeLists.txt).
A new case needs an `add_subdirectory` line there.
The sources in `common/BOUTExtensions/` are not compiled on their own: the
headers include them, and the models include the headers.

## Build, test and lint as CI does

Needs `BOUT++` 5.2.0 on `CMAKE_PREFIX_PATH`.
The image `ghcr.io/celma-project/celma-env` (built from
[docker/Dockerfile.bout5](docker/Dockerfile.bout5)) has it.
[test.yml](.github/workflows/test.yml) pins the image by digest.

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
(cd tests/radialLowPass && ./radialLowPass -d data)
python3 tests/regression/regression.py
python3 tests/regression/mes.py
python3 tests/numerics/fingerprint.py
python3 -m pytest tests/numerics/test_guard.py
pre-commit run --all-files
```

* The executables land next to their sources, not in `build/`.
* The regression references in [tests/regression](tests/regression) depend on
  the `BOUT++` version and the compiler.
  CI records them in the image, so a mismatch on another machine is not
  necessarily a bug.
  Both scripts write the numbers of the run to `recorded-*.json` in the
  current directory.
* `pre-commit` needs `cppcheck` from the system package manager.

## Conventions

* Lint is [.pre-commit-config.yaml](.pre-commit-config.yaml): the standard
  `pre-commit-hooks`, `ruff-check` and `ruff-format`, `clang-format`
  (reads [.clang-format](.clang-format)), `cpplint` and `cppcheck`.
  Fix findings in the code; do not loosen a hook or its filters.
* Ruff is configured in [pyproject.toml](pyproject.toml): line length 88,
  Python 3.12, numpy docstring convention.
  The `ignore` list there is known debt, not a model to copy.
* C++ follows `.clang-format` (LLVM based).
* Python code in new files needs numpy-style docstrings.
* Notebooks are committed without output (see
  [.codeCleaner/README.md](.codeCleaner/README.md)).

## Numerics need Michael's consent

Never change any of these without Michael's explicit consent for that specific
change:

* model equations, operators, derivative schemes or stencils
* boundary conditions, the Jacobian or the metric
* solver or inversion settings, numeric input options
* the `BOUT++` version (the pinned image and the compiler flags)
* a rerecorded numerics reference

Ask for it as its own question, with the numbers or the line before and after.
Consent shows as the label `numerics-change-approved`, which Michael applies to
the pull request; a push to the branch removes it, so it binds to the head he
saw.
The check `numerics-guard` ([numerics-guard.yml](.github/workflows/numerics-guard.yml),
[tests/README.md](tests/README.md)) fails any pull request that changes
the guarded content without it.
It ignores formatting, comments and include order, so a formatting commit
passes.
Do not edit `tests/numerics/` or the guard workflow to get past it.

Where the numerics live:

* Model equations: [celma/](celma/celma.cxx),
  [celmaWithBoussinesqApprox/](celmaWithBoussinesqApprox/celmaWBA.cxx) and
  [common/BOUTExtensions/src/ownOperators.cxx](common/BOUTExtensions/src/ownOperators.cxx).
  The MES cases and the tests exist to catch changes here; run them after any
  edit.
* `BOUT.inp` values: they set the physics of a scan and the inputs of the MES
  cases.
  `tests/regression/mes.py` reads every `BOUT.inp` under `MES/`.
  The z derivatives are FFT, as in the thesis (v1.0.0).
* [tests/regression](tests/regression) and [tests/numerics](tests/numerics)
  reference files: rerecord them only with consent, for a reason you can
  state, in the CI image.
* `derivations/`: the results are hard-coded in the C++, so changing a
  notebook does not change the simulation.

## Python scripts that do not run here

`PBS*.py`, `testScan.py`, `testPostProcessing.py`, the `driver*.py` files in
`MES/` and `tests/radialLowPass/driverTest.py` import `bout_runners` and
`CELMAPy`.
`CELMAPy` was moved to a separate repository (CELMA-py) and neither package
is installed by the CI image.
The scripts were written for the earlier makefile build.
Do not treat them as a working entry point.

## Glossary

* **BOUT++**: framework for fluid plasma simulations in curvilinear
  geometry; the models are `PhysicsModel`s.
* **CELMA**: the cylinder model in this repository; the name of the project.
* **CSDX**: the linear plasma device whose setup the scan inputs
  (`CSDX*`) follow.
* **Boussinesq approximation**: the variant in `celmaWithBoussinesqApprox/`
  (directories named `Bous*`); the model in `celma/` does not use it.
* **MES**: method of manufactured solutions.
  Pick an analytic solution, compute the source it needs, run the code with
  that source and check that the error falls at the expected order as the
  grid is refined.
  The notebooks in `MES/*/*/calculations/` produce the sources.
* **Own operators / own BCs**: the differential operators and boundary
  conditions in `common/BOUTExtensions/`, written to handle the singularity
  at the cylinder axis.
* **Naulin solver**: Laplace inversion that `celma/` uses to get the potential
  from the vorticity (`NaulinSolver` in `ownLaplacianInversions`).
* **Cauchy BC**: boundary condition that sets both the value and the
  derivative at the boundary.
* **Scan**: a set of runs with one parameter varied; the `CSDX*Scan*`
  directories hold the input of each.
