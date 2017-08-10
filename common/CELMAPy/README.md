# CELMAPy

`python` modules used for pre- and post-processing of the simulations.
The usual structure for these files is

```
__init__.py         # Simplifies import outside of the module
collectAndCalcX.py  # Procedures for collecting and calculating something to be plotted
driverX.py          # The driver which calls the procedures in collecAndCalcX and parse the result to plotX.py
plotX.py            # The plotting routine
```

[blobs](blobs) - Procedures which collects and plots conditional averaging and counting of blobs and holes
[calcVelocities](calcVelocities) - Procedures for calculating ExB velocities
[collectAndCalcHelpers](collectAndCalcHelpers) - Procedures which are commonly used in the `collectAndCalc*.py` files
[combinedPlots](combinedPlots) - Procedures which aquires and plots the fluctuation, the PDF and PSD for three different positions
[driverHelpers](driverHelpers) - Procedures which are commonly used in (all kinds of) drivers (such as the `PBSSubmitter`)
[energy](energy) - Procedures which collects and plots the energy
[fields1D](fields1D) - Procedures which collects and plots profiles
[fields2D](fields2D) - Procedures which collects and plots 2D cuts of the data
[fourierModes](fourierModes) - Procedures which collects and plots the Fourier modes
[growthRates](growthRates) -  Procedures which collects and plots the growth rates (both analytically and from the simulations)
[logReader](logReader) - Contains a module which reads the log files
[MES](MES) - Procedures used in the `MES` routines
[modelSpecific](modelSpecific) - Defines which fields to be collected, and in which order for profile plots
[PDF](PDF) - Procedures which collects and plots the Probability Distribution Function
[performance](performance) - Procedures which collects and plots the performance
[plotHelpers](plotHelpers) - Procedures which are commonly used in the `plot*.py` files
[poloidalFlow](poloidalFlow) - Procedures which collects and plots the poloidal flows
[PSD](PSD) - Procedures which collects and plots the Power Spectra Density
[radialFlux](radialFlux) - Procedures which collects and plots the radial flux
[radialProfile](radialProfile) - Procedures which collects and plots the comparison between steady state and turbulence profiles
[repairBrokenExit](repairBrokenExit) - Procedure which repairs the dump file if they are exited badly (for example the cluster stops the job in the middle of a write)
[scanDriver](scanDriver) - The driver used for running scans with `bout_runners`
[scripts](scripts) - Contains non-callable scripts
[skewnessKurtosis](skewnessKurtosis) - Procedures which collects and plots the poloidal skewness and kurtosis as a function of radius
[superClasses](superClasses) - Contains all the parent classes (note that this is maybe not the best structure)
[timeTrace](timeTrace) - Procedures which collects and plots a single time trace
[totalFlux](totalFlux) - Procedures which collects and plots the total particle flux out of the system
[unitsConverter](unitsConverter) - Procedures for converting between normalized and physical units
