# celma

Folder with the different versions of the code.
Each run of group of run has its own driver, and uses `bout_runners` for
running and post-processing.

* 1-CELMA - Original code
* 1.0-CELMATopDown - CELMA with cheats to make it run
* 1.1-CELMACauchy - Original code with cauchy BC on the density
* 2-CELMACollisionDamper - Restarting from 1-CELMA, and let the collisionality
  fall slowly to $0$
* 3.0.0-CELMAIMEXOnlyDiffusion - CELMA using the IMEX solver (cheating to make
  work)
* common - python post processing and own implementations to BOUT++
* initialProfiles - A trial to analytically calculate the initial profiles
