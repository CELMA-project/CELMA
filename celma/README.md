# celma

Folder with the different versions of the code.

Each run of each group has its own driver. `bout_runners` is used for running
and post-processing.

Original forumlation is used here

* 1-CELMA - Original code
* 1.0-CELMATopDown - CELMA with cheats to make it run
* 1.1-CELMACauchy - Original code with cauchy BC on the density
* 2-CELMACollisionDamper - Restarting from 1-CELMA, and let the collisionality
  fall slowly to $0$
* 3.0.0-CELMAIMEXOnlyDiffusion - CELMA using the IMEX solver (cheating to make work)
* 4-CELMARadialNeumann/ - Using only neumann conditions
* 4.1-CELMARadialNeumannIMEX - Formulated in IMEX formulation
* 5-CELMAWNoise - Adding noise to the runs
* 5.1-CELMANoiseFiltering1-3 - Adding spectral filtering (accidentially filtering 2/3 of the modes)
* 5.1.1-CELMANoiseFiltering2-3 - Spectral filter away 1/3 of the top modes
* 5.2-CELMAWNoisFilterRho - Spectral filter which filter away the top 1/3 modes resolvable on the outermost radius
* 6.0-StripCELMAVelAdvection - Stripping the system to identify fast waves
* 6.1-StripCELMAVelAdvectionVecAdv - Stripping the system to identify fast waves
* 7-CELMAAnnulus - Checking for difference if inner is treated as annulus
* 8-CELMASplit - New IMEX formulation
* 8.1-CELMASplitPerpArt0.25 - Compare start-up time with 8 using artPerp = 0.25
* 8.2-CELMASplitNoNInArtVisc - Compare start-up time with 8 when changing artificial viscosity for parallel currents
* 8.3-CELMASplitCleanUp - Clean-up code after finding good parameters
* 8.4-CELMASplitHypervisc - As 8.3, but made viscosities changable to hyper viscosities
* 9-CELMAUpdateArtViscCleanUp - Not damp on evolved variables, but rather at density, added artificial viscosity on OmegaD
* 9.1-CELMAUpdateKeepDivN - Same as 9, but changed artificial viscosity for parallel currents
* common - python post processing and own implementations to BOUT++
* initialProfiles - A trial to analytically calculate the initial profiles
