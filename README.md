# CELMA

[![Lint](https://github.com/CELMA-project/CELMA/actions/workflows/lint.yml/badge.svg?branch=main)](https://github.com/CELMA-project/CELMA/actions/workflows/lint.yml)
[![Test](https://github.com/CELMA-project/CELMA/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/CELMA-project/CELMA/actions/workflows/test.yml)
[![Image](https://github.com/CELMA-project/CELMA/actions/workflows/image.yml/badge.svg?branch=main)](https://github.com/CELMA-project/CELMA/actions/workflows/image.yml)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FCELMA-project%2FCELMA.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FCELMA-project%2FCELMA?ref=badge_shield)

Repository for the [CELMA](https://celma-project.github.io/) code: a
[`BOUT++`](https://boutproject.github.io/) model of plasma turbulence in a
cylinder, with the code that verifies its operators.

> **NOTE**: The C++ builds with CMake against `BOUT++` 5.2.0, and the
> continuous integration builds every target and runs the tests in
> [Build and test](#build-and-test).
> The Python drivers (`PBS*.py`, `testScan.py`, `driver*.py`) still use the
> old `bout_runners` workflow and need the `CELMAPy` package, which lives in
> a separate repository (CELMA-py).
> They are not run by the continuous integration and were not updated for the
> CMake build.
> The repository is not actively maintained; open an
> [issue](https://github.com/CELMA-project/CELMA/issues) to ask for a fix.
> See the [releases](https://github.com/CELMA-project/CELMA/releases) for
> older, stable code.

Agents: start at [AGENTS.md](AGENTS.md).

* [Repository structure](#repository-structure)
* [Build and test](#build-and-test)
* [Contribute](#contribute)
* [Issues](#issues)
* [License](#license)

## Repository structure

* [`celma/`](celma/README.md)
    * The CELMA model: source, input files of the scans, run scripts.
* [`celmaWithBoussinesqApprox/`](celmaWithBoussinesqApprox/README.md)
    * The CELMA model with the Boussinesq approximation.
* [`common/`](common/README.md)
    * Own `BOUT++` extensions in C++ (operators, boundary conditions,
      filters, monitors, Laplace inversion, parameters) shared by the models
      and the cases below.
* [`derivations/`](derivations/README.md)
    * Derivations of stencils, boundary conditions, collisionalities and
      dispersion relations, mostly in
      [jupyter notebooks](http://jupyter.org/).
      The results are copied into the code by hand.
* [`MES/`](MES/README.md)
    * Cases that verify the operators and boundary conditions with the method
      of manufactured solutions (MES).
* [`tests/`](tests/README.md)
    * The filter test and the regression tests the continuous integration
      runs.
* [`docker/`](docker/Dockerfile.bout5)
    * `Dockerfile.bout5` builds the `BOUT++` 5.2.0 environment the continuous
      integration uses.
      `Dockerfile` is the old Ubuntu 18.04 image from before the CMake build
      and is not used by any workflow.
* [`.github/`](.github/workflows)
    * Workflows: `lint.yml`, `test.yml` and `image.yml`.
* [`.codeCleaner/`](.codeCleaner/README.md)
    * Old helper scripts for formatting, converting and cleaning the code.

The build is defined in [CMakeLists.txt](CMakeLists.txt) and one
`CMakeLists.txt` per model and case.
The lint configuration is in [.pre-commit-config.yaml](.pre-commit-config.yaml),
[pyproject.toml](pyproject.toml) (ruff) and [.clang-format](.clang-format).

## Build and test

Build against an installed `BOUT++` 5.2.0, or use the image the continuous
integration uses (see [docker/Dockerfile.bout5](docker/Dockerfile.bout5)):

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/bout-install
cmake --build build
```

Each executable is written next to its source and input files, for example
`celma/celma` with `celma/CSDXNyScan/BOUT.inp`.
Run a model as described in the
[`BOUT++` user manual](https://bout-dev.readthedocs.io/en/latest/), for
example `./celma -d CSDXNyScan` from `celma/`.

The tests, as [test.yml](.github/workflows/test.yml) runs them, after the
build:

```bash
cd tests/radialLowPass && ./radialLowPass -d data   # filter test
cd ../..
python3 tests/regression/regression.py              # both models
python3 tests/regression/mes.py                     # every MES case
```

The lint, as [lint.yml](.github/workflows/lint.yml) runs it:

```bash
pre-commit run --all-files
```

## Contribute

Contributions are welcome through pull requests.
Run `pre-commit run --all-files` first.

## Issues

Bugs, issues and questions are handled at the
[issue tracker](https://github.com/CELMA-project/CELMA/issues).

## License

`CELMA` is licensed under the MIT license, see [LICENSE](LICENSE).

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FCELMA-project%2FCELMA.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FCELMA-project%2FCELMA?ref=badge_large)
