# common

Contains own extensions to `BOUT++` together with pre- and post-processing
modules written in `python`.

* [BOUTExtensions](BOUTExtensions) - Own extensions to BOUT++.
* [CELMAPy](CELMAPy) - `python` modules used for pre- and post-processing of
  the simulations.
* [standardPlots](standardPlots) - Contains the `PlotSubmitter`-class, which
  easen the post-processing calls for plotting.
  See [`testPostProcessing.py`](../celma/testPostProcessing.py) for an example
  of usage.
  This folder also contains the standarized caller-functions used by the
  `PlotSubmitter` to the post-processing modules using standard parameters.
