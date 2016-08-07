#ifndef __HELPERFUNCTIONS_CXX__
#define __HELPERFUNCTIONS_CXX__

#include "../include/helperFunctions.hxx"


/*!
 * The energy monitor is a wrapper which calls functions which
 * calculates the monitored energy. The monitors accept only the input
 * arguments listed.
 *
 * \param[in] solver   Pointer to the solver
 * \param[in] simtime  Current time
 * \param[in] iter     Current iteration
 * \param[in] NOUT     Number of outputs
 *
 * \returns 0 On success
 */
int energyIntMon(Solver *solver, BoutReal simtime, int iter, int NOUT)
{
    TRACE("Halt in energyIntMon");

    Field3D polAvgN;
    Field3D polAvgPhi;
    Field3D polAvgUEPar;
    Field3D polAvgUIPar;

    // Call the functions
    // FIXME: THESE ARE NOT GLOBAL VARIABLES...HOW TO MAKE THE MONITOR SEE
    // THEM?
    // IDEA: CONSTRUCTOR WHICH TAKES THE CLASS/OBJECT AS INPUT
    // IDEA: fullE must be in the model
    fullEnergy(fullE, n, phi, uEPar, uIPar);
    // Calculate averages
    polAvgN = polAvg(n);
    polAvgN = polAvg(phi);
    polAvgN = polAvg(uEPar);
    polAvgN = polAvg(uIPar);
    // Calculate energy in fluctuations
    polAvgEnergy(polAvgE, n, phi, uEPar, uIPar);

    return 0;
}

/*!
 * Returns the poloidal average of a field
 *
 * \param[in] f     The field to take the average of
 * \param[in] xInd  The rho index to take the average at
 * \param[in] yInd  The parallel index to take the average at
 *
 * \returns result The poloidal average of the field
 */
Field3D const polAvg(const Field3D &f)
{
    TRACE("Halt in polAvg");

    Field3D result = 0.0;
    BoutReal avg;

    for(int xInd = mesh->xstart+1; xInd <= mesh->xend-1; xInd++){
        for(int yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
            // Find the average
            avg = 0.0;
            for(int zInd = 0; zInd < mesh->ngz -1; zInd ++){
                avg += f(xInd, yInd, zInd);
            }
            avg /= (mesh->ngz - 1);

            // Subtract the average from the field
            for(int zInd = 0; zInd < mesh->ngz -1; zInd ++){
                result(xInd, yInd, zInd) = f(xInd, yInd, zInd) - avg ;
            }
        }
    }

    return result;
}

#endif
