#ifndef __OWNMONITORS_CXX__
#define __OWNMONITORS_CXX__

#include "../include/ownMonitors.hxx"

/*!
 * This function is used instead of a constructor as a PhysicsModel object is
 * used as an input.
 *
 * \param[in] model Object with the simulated model
 */
template<class PhysicsClass>
void ownMonitors<PhysicsClass>::create(PhysicsClass *inputModel)
{

    TRACE("Halt in ownMonitors::create");

    model = inputModel;
    //std::cout << inputModel->phi(0,0,0) << std::endl;
    //std::cout << model->phi(0,0,0) << std::endl;
}


/*!
 * The energy monitor is a wrapper which calls functions which
 * calculates the monitored energy. The monitors accept only the input
 * arguments listed.
 *
 * \note This must be defined as static in order to get the correct type
 *       http://stackoverflow.com/questions/15841338/c-unresolved-overloaded-function-type
 *
 * \param[in] solver   Pointer to the solver
 * \param[in] simtime  Current time
 * \param[in] iter     Current iteration
 * \param[in] NOUT     Number of outputs
 *
 * \returns 0 On success
 */
template<class PhysicsClass>
int ownMonitors<PhysicsClass>::energyIntMon(Solver *solver, BoutReal simtime, int iter, int NOUT)
{
    TRACE("Halt in ownMonitors::energyIntMon");

    std::cout << model->phi(0,0,0) << std::endl;
    std::cout << "foo" << std::endl;
//    // Call the functions
//    // FIXME:
//    See 6.2 for previous monitors
//    fullE must be in the model
//    fullEnergy(fullE, n, phi, uEPar, uIPar);
//    // Calculate averages
//    polAvgN = polAvg(n);
//    polAvgN = polAvg(phi);
//    polAvgN = polAvg(uEPar);
//    polAvgN = polAvg(uIPar);
//    // Calculate energy in fluctuations
//    polAvgEnergy(polAvgE, n, phi, uEPar, uIPar);

    return 0;
}

#endif
