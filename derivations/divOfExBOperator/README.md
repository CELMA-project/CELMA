# divOfExBOperator

[Root](../../README.md) | [derivations](../README.md)

This folder contains the derivation of the term
`div(u_E dot grad(grad_perp(phi/B) n))`.

**NOTE**: The derivation can also be found in appendix J of the
[thesis](https://github.com/CELMA-project/dissertation/releases/download/v1.1.x/17_PhD_Loeiten.pdf).

**NOTE**: The operator is expanded in single terms in the notebooks, however,
the expanded form is never used in the celma code due to stability issues.

* [clebschVector.py](clebschVector.py) - Handles the coordinates
* [common.py](common.py) - Contains common variables and functions for all
  coordinate systems
* [derivationByHand.ipynb](derivationByHand.ipynb) - Doing the derivation by
  "hand" (checks that own implementation of the vector classes is working as
  expected)
* [divOfVectorAdvectionWithN.ipynb](divOfVectorAdvectionWithN.ipynb) -
  Derivation using own implementation
* [manual_simplification.txt](manual_simplification.txt) - Output from `divOfVectorAdvectionWithN.ipynb`
