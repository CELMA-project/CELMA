# CELMA

[![Build Status](https://travis-ci.org/CELMA-project/CELMA.svg?branch=master)](https://travis-ci.org/CELMA-project/CELMA)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FCELMA-project%2FCELMA.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FCELMA-project%2FCELMA?ref=badge_shield)

> **NOTE**: This repo is needs some refactorings to work out of the box.
I'm currently not actively maintaining this repo, but could implement some of the fixes upon request.
Feel free to do so by opening a new [issue](https://github.com/CELMA-project/CELMA/issues)
Refer to [releases](https://github.com/CELMA-project/CELMA/releases) for more
stable code

Repository for the [CELMA](https://celma-project.github.io/) code.

* [Repository structure](#repository-structure)
* [Install](#install)
* [Usage](#usage)
* [Contribute](#contribute)
* [Issues](#issue)

## Repository structure

* `celma/`
    * The code for solving the CELMA model.
* `celmaWBoussinesq/`
    * The code for solving the CELMA model using the Boussinesq approximation.
* `common/`
    * Common files for extra implementations (not yet included in BOUT++)
      written in `c++` and pre- and post-proccesing functions written in `python3`.
* `derivations/`
    * Derivations of the numerics used in the code, namely:
        * Boundary polynomials
        * Own operators
        * Collisionalities
        * Analytic dispersion relations
    * Mostly written using [jupyter-notebooks](http://jupyter.org/).
* `MES/`
    * Method of exact solutions used to verify own operators.
* `tests/`
    * Test used to check that the implemntations are working.

## Install

An install script can be found in [install/boutPPInstall.sh](install/boutPPInstall.sh).
This will install the master branch of `BOUT-dev` in `$HOME`, whereas the
dependencies will be installed in `local`.

**NOTE** : The `makefiles` in this repository assumes that `BOUT++` is located
           in the `$HOME` directory.

## Usage

The models can be simulated in the way described in the [`BOUT++` user-manual](http://bout-dev.readthedocs.io/en/latest/).
For more advanced jobs, such as parameter scans using supercomputers, one can
use the `PBS<name>.py`-scripts located in `celma` and `celmaWBoussinesq`.

## Contribute

Contributions are more than welcome through pull requests.

## Issues

Bugs, issues and questions are being handled at the [issue-tracker](https://github.com/CELMA-project/CELMA/issues).


## License
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FCELMA-project%2FCELMA.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FCELMA-project%2FCELMA?ref=badge_large)
