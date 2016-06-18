
run-01

  Large terms in uEpar equation
   uEParRes  (max ddt = 2.6)
   gradPhiLnN (3.9)
   uEparParArtVisc (2.6)

  Vorticity equation terms small (1e-3)

run-02
  Removed:
    vortDParArtVisc
    uEParParArtVisc
    uEParPerpArtVisc

  -> Still small (perhaps worse?)

run-03
  Changing to nz=1 using same input as run-01

  -> Runs quite well. Erratic time step, but much faster

run-04
  Expanding run-03 to 129 points (same as run-01). Steps of 0.1

  1.000e+03        150       7.38e+00    66.9    2.1    7.0    1.7   22.4
  1.000e+03         40       2.64e+00    72.0    1.2    2.8   10.3   13.7
  1.000e+03         32       1.90e+00    59.6    1.7    3.0    7.1   28.7
  1.000e+03         31       1.95e+00    43.9    4.1    5.3    6.3   40.4
  1.001e+03         23       1.46e+00    42.1    4.6    6.9    1.9   44.6
 
run-05
  Adding noise to the expanded run-03. Step size 0.1
  >>> addnoise(var="vortD")

  1.000e+03        291       1.39e+01    62.4    2.3    6.8    2.2   26.3
  1.000e+03        110       6.22e+00    51.1    1.1   26.3    4.1   17.5
  1.000e+03        127       6.28e+00    68.7    1.5   12.0    1.6   16.3
  1.000e+03        102       5.09e+00    58.2    2.3    5.8    5.9   27.8
  1.001e+03         93       5.13e+00    52.1    2.1   14.2    5.9   25.7
 
run-06
  Same as run-05, but removing filtering (maxmode=-1 in source)

  1.000e+03          2       3.82e-01    23.4    0.7    7.8   83.5  -15.4
  1.000e+03        285       1.25e+01    64.7    1.4    9.8    1.0   23.0
  1.000e+03        144       6.96e+00    62.9    2.0    6.6    2.6   26.0
  1.000e+03        165       9.13e+00    52.2    1.0   23.4    3.0   20.4
  1.000e+03         88       4.50e+00    62.6    1.5    2.1    6.9   27.0
  1.001e+03         94       4.80e+00    58.1    0.8   17.0    4.6   19.5

run-07
  Same as run-06, but larger output timestep of 1

  1.001e+03       1449       7.73e+01    51.9    4.0    4.8    0.3   38.9
  1.002e+03       1428       7.72e+01    57.9    1.8   23.1    0.8   16.4
  1.003e+03       1215       6.42e+01    57.3    2.5   14.3    0.6   25.2
  1.004e+03        944       5.14e+01    55.0    4.3    4.9    0.6   35.1

run-08
  Same as run-07, adding D4DZ4 term to vorticity:

    -0.01*D4DZ4(vortD)*SQ(SQ(mesh->dz))

  1.001e+03       1397       7.49e+01    53.2    1.5   24.7    0.5   20.1
  1.002e+03       1406       7.82e+01    54.3    1.8   27.1    0.4   16.3
  1.003e+03       1263       6.75e+01    55.5    0.9   26.1    0.3   17.2
  1.004e+03        794       4.35e+01    68.8    1.7   16.0    0.5   13.0

  Runs, growing modes, but still 500 - 1000 iterations per /wci


run-09
  Same as run-08. Changing nuEI to 1 (from 25), 
  removing parallel viscosity from vortD

run-10
 Starting from scratch with nz=1. Removing parallel viscosity from uEpar
 -> Failed before first output (t=5)

 Changing j term in vorticity:

 //     + nGradUiUe
 //     + uiUeGradN
  + Div_par(n*(uIPar - uEPar))

 -> Runs, but still takes many iterations after reaching steady state

run-11
  Restarting from run-10. Adding ddt() vars to output
  Printing amax(abs("ddt(*)")) 
  
  ddt(lnN)   : 1.86067591418e-07
  ddt(vortD) : 2.20509021116e-08
  ddt(uEPar) : 8.74309369414e-05
  ddt(uIPar) : 4.52466732229e-05

  but number of iterations is ~800 per 1/wci

run-12
  Restarting from run-10, using PETSc

  mpirun -np 16 nice -n 10 ./celma -d a-data/ restart solver:type=petsc solver:start_timestep=1 -ts_type theta -ts_theta_theta 1 -ksp_atol 1e-7 -snes_type newtonls -{ksp,snes}_monitor

  Runs very fast when it runs (~ 5 - 20 iterations per 5 wci), but quite fragile

  snes_type "newtonls" or "qn" both work well

run-13
  Restarting from run-10, using imexbdf2

  order 2: "qn" snes type runs with dt=1 or 10, but fails to converge with dt=0.1


  # 2nd order, timestep=1, not matrix free, not using coloring
  # Jacobian only computed once at the start, then re-used

  mpirun -np 16 nice -n 10 ./celma -d a-data/ restart -snes_lag_jacobian -2

  3.280e+03       1548       4.55e+00    20.1    2.1   14.6    0.2   63.0
  3.285e+03         11       9.56e-01     1.2    0.1    0.1    1.7   96.9
  3.290e+03         11       9.45e-01     1.2    0.1    0.3    2.4   96.0
  3.295e+03         11       2.09e+00     0.6    0.1    0.4    0.8   98.2
  3.300e+03         10       8.31e-01     1.3    0.1    0.2    2.8   95.7
  3.305e+03          7       2.75e-01     2.5    0.5    0.5    5.6   90.9
  3.310e+03          7       2.59e-01     3.0    0.1    0.2    8.3   88.2

run-14
  Starting from scratch, nz=1, using imexbdf2 method.
  
  mpirun -np 16 nice -n 10 ./celma -d a-data/

  Runs, but needs a large number of evaluations to calculate the Jacobian e.g.:

  2.000e+00      15419       1.82e+02    21.0    5.1   12.3    0.0   61.6


run-15
  Like run-10, using pvode with nz=1. Trying with nuEI = 0.1
  -> Runs ok

  5.000e+00       3701       5.15e+00    67.1    5.7   14.7    0.5   11.9
  1.000e+01       2985       3.84e+00    78.4    5.2    7.9    0.7    7.8
  1.500e+01       3167       4.09e+00    47.5    4.8   27.7    1.5   18.5
  2.000e+01       3208       4.15e+00    83.6    5.0    5.0    0.3    6.1
  2.500e+01       3431       4.37e+00    55.3    4.9   22.3    0.7   16.8

run-16
  Like run-15, reducing nuEI to 0.01. Also runs ok.

  5.000e+00       3564       4.80e+00    41.5    7.9   24.3    0.3   26.1
  1.000e+01       3181       3.88e+00    67.9    5.9    6.1    0.4   19.7
  1.500e+01       3072       4.09e+00    42.3    5.9   13.1    0.5   38.2
  2.000e+01       3106       3.91e+00    46.7    6.1   16.5    0.4   30.3
  2.500e+01       3838       4.70e+00    46.1    5.1   25.2    0.5   23.0

  Checking magnitude of the variables

  lnN     2.9
  uEPar   1.0
  uIPar   0.98
  vortD   0.007

  so vorticity seems to be ~200 times smaller than other variables.

run-17

  Dividing vortD by 200 at start of rhs(), multiplying ddt(vortD) by 200
  at the end. 

  5.000e+00       3940       5.15e+00    82.9    5.0    5.9    0.2    6.0
  1.000e+01       3237       4.08e+00    48.1    4.7    7.2    0.4   39.6
  1.500e+01       3124       4.05e+00    47.3    5.5   12.5    0.5   34.1
  2.000e+01       3219       4.09e+00    70.5    4.4   12.2    0.5   12.4
  2.500e+01       3517       4.39e+00    55.4    4.8   18.8    0.7   20.3

  -> Seems to make little difference

run-18
  Starting from run-16, expanding and adding noise
  IDL> expand_restarts, 129, path="run-16"
  
  >>> from boutdata.restart import addnoise
  >>> addnoise(var="vortD")

  -> run-16 vorticity is vortD[-1,4,3,0] =  -0.0065059924682432137

  -> Starting vorticity is ok: vortD[0,4,3,0] = -0.0065325484611093998
     (noise ~ 1e-5 added)
  
  -> Very rapidly changes

run-19
  -> Restarting from run-16, nz=1

   vortD[:,4,3,0] = [-0.00650599, -0.00650599, -0.00650599, -0.00650599,..

   so not changing

run-20
  Starting from run-16, expanding but not adding noise
  -> vortD changes 

  vortD[:,4,3,0] = [-0.00650599, -0.00849473, -0.01292647, -0.01568719, -0.0158599 ,
       -0.01532865, -0.01953092, -0.0204429 , -0.0196575 ]

  (output timestep 1/wci)

run-21
  Starting from run-16, not expanding. Adding some additional outputs

  -> Difference between run-20 and run-21 seems to be in term:
   vortDPerpArtVisc           =   artViscPerpVortD*Laplace_perp(vortD);

  Much smaller in case with nz=128 than when nz=1!

run-22:
  Changing Laplace_perp to Delp2
  -> Same result: vortDPerpArtVisc much smaller for nz=128
  At x=4, almost constant ratio of 332.7
  Vorticity is the same to within ~1e-7
  -> Difference in metric?

  -> Difference in artViscPerpVortD 
  
  line 96:
  artViscPerpVortD *= SQ(mesh->dx(0,0) + mesh->dz);

  ** Changing to 

  artViscPerpVortD *= SQ(mesh->dx(0,0));

run-23
  Like run-16, starting from scratch with nz=1

  Runs, but slows down:
  
  5.000e+00       3958       6.39e+00    76.6    7.6    7.2    0.2    8.3
  1.000e+01       3423       8.86e+00    44.5   10.0   18.4    0.2   26.9
  1.500e+01       3198       7.55e+00    46.0   11.1   17.5    0.3   25.1
  2.000e+01       3041       7.19e+00    52.0    7.9   11.8    0.2   28.1
  ...
  4.990e+03      13491       1.85e+01    76.3    3.8    7.9    0.1   11.8
  4.995e+03      12227       1.73e+01    78.3    3.7    5.0    0.1   12.8
  5.000e+03      14255       2.13e+01    63.8    3.7   11.2    0.1   21.2

  Slower than comparable run-16:

  run-16: 1.630e+03       3214       3.87e+00    73.0    3.3    3.8    0.6   19.3
  run-23: 1.630e+03      12439       1.75e+01    73.2    3.8    8.4    0.1   14.5

  -> Difference is in magnitude of perpendicular diffusion, due to dz factor

Checking sensitivity to parameters. Run for dt=1, record # iterations
   
  No changes: 3216
  artPar reduced (5 -> 2.5):  1164
  artPerp reduced (0.5 -> 0.25): 1473
  mu reduced (1836.36 to 900) :  978
  nuEI reduced (0.01 to 0.005) :  1497
  artViscParLnN reduced (5 to 2.5) :  1436
  artViscParUIPar reduced (5 to 2.5) : 1179
  artPar turned off (5 -> 0) : 1308
  Adding artPar to artViscParUEPar : 552
     Reducing artPerp (0.5 to 0.25) : 542
  Changing parallel viscosity terms, removing "/n" : 523
    Removing "/n" from perpendicular viscosity : 607  (**)

run-24: Running with modified perp. & par. viscosity terms
  restarting from run-23
  
  Starts much faster:

  5.005e+03       2887       4.25e+00    42.2    8.3   31.3    0.5   17.7
  5.010e+03       2342       3.35e+00    42.2    9.8   22.8    0.6   24.6
  5.015e+03       1732       2.69e+00    39.1   10.0   32.6    1.2   17.0

  ... but slows down

  5.490e+03       4740       6.45e+00    43.2    6.8   27.8    0.3   21.9
  5.495e+03       4987       6.82e+00    45.7    5.0   32.1    0.3   16.7
  5.500e+03       4916       6.62e+00    81.5    7.2    6.5    0.3    4.6

Checking sensitivity

  No change: 1100
  Turning off artViscPerpVortD: 458
  Turning off all artPerp : 527
  artperp changed (0.25 -> 0.1) : 516
  Reducing artPar (5 -> 1) : 1196

run-25: 
  Restarting from run-24, removing artViscPerpVortD

run-26:
  Restarting from run-25, expanding and adding noise

run-27
  Continuing from run-26

run-28
  Continuing from run-26, changing Z derivatives from FFT to C2

run-29
  Continuing from run-26, increasing D4DZ4 coefficient from 0.01 to 0.1

run-30
  Continuing from run-26, removing ExB advection of vorticity 

run-31
  run-26 restart, changing ExB advection to Arakawa
  -> Runs well, saturated turbulent state

run-32
  nz=1 simulation with double size machine

run-33

  IDL> expand_restarts, 257, path="run-32"
  >>> from boutdata.restart import addnoise
  >>> addnoise(var="vortD")  


