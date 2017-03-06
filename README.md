# CELMA

Repository for the [CELMA](https://celma-project.github.io/) code.

* [Repository structure](#Repository-structure)
* [Install](#Install)
* [Contribute](#Contribute)
* [Issues](#Issue)

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

Requirements:

* The same requirements needed for [BOUT++-3.1](https://github.com/boutproject/BOUT-dev/releases/tag/v3.1)-release
* `bout_runners.py` found [here](https://github.com/CELMA-project/CELMA/releases/download/v0.1beta/bout_runners.py)
* For the pre- and post-processing:
    * `python3`. [Installation guide](https://github.com/loeiten/usingLinux/blob/master/installationProcedures/python.md)
      for an easy installation guide
    * `matplotlib` (version `2.0.0` has been used in the thesis).
    * `ffmpeg` for animaions. [Installation guide](https://github.com/loeiten/usingLinux/blob/master/installationProcedures/ffmpeg.md).

Installation

1. Install BOUT++ as explained in the `user_manual` of BOUT++, or in [this manual](https://github.com/loeiten/usingLinux/blob/master/installationProcedures/BOUT-dev.md).
   **NOTE** : The `makefiles` assumes that `BOUT++` is located in the `home`
   directory.
2. Copy the attached `bout_runners.py` to `<yourBOUT++Installation>/tools/pylib/bout_runners`.
3. `make` using the `makefiles` located in `celmaRC2`.

The models can be simulated in the way described in the `BOUT++`-`user_manual`.
For more advanced jobs, such as parameter scans using supercomputers, one can
use the `PBS<name>.py`-scripts located in `celma` and `celmaWBoussinesq`.

## Contribute

Contributions are more than welcome through pull requests.

## Issues

Bugs, issues and questions are being handled at the [issue-tracker](https://github.com/CELMA-project/CELMA/issues).
