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
}

/*!
 * The energy monitor is a wrapper which calls functions which
 * calculates the monitored energy. The monitors accept only the input
 * arguments listed.
 *
 * \note This must be defined as static in order to get the correct type
 *       (remember that static members belong to the class rather than the
 *       object)
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
    TRACE("Halt in ownMonitors<PhysicsClass>::energyIntMon");

    // Call the fillEnergy function
    int lol;
    kinEnergy(model->kinE, model->n, model->phi, model->uEPar, model->uIPar);

    std::cout << model->phi(0,0,0) << std::endl;
//    // Call the functions
//    // FIXME:
//    See 6.2 for previous monitors
//    kinE must be in the model
//    fullEnergy(kinE, n, phi, uEPar, uIPar);
//    // Calculate averages
//    polAvgN = polAvg(n);
//    polAvgN = polAvg(phi);
//    polAvgN = polAvg(uEPar);
//    polAvgN = polAvg(uIPar);
//    // Calculate energy in fluctuations
//    polAvgEnergy(polAvgE, n, phi, uEPar, uIPar);

    return 0;
}

/*!
 * The energy monitor is a wrapper which calls functions which
 * calculates the monitored energy. The monitors accept only the input
 * arguments listed.
 *
 * \note This must be defined as static in order for the static energyIntMon to
 *       see it
 *       http://stackoverflow.com/questions/15235526/the-static-keyword-and-its-various-uses-in-c
 *
 * ## Derivation
 * For a single particle, we have \f$E_{kin} = \frac{1}{2}m\mathbf{v}^2\f$,
 * whereas for a fluid, we have \f$E_{kin} = \frac{1}{2}m\int n\mathbf{u}^2 dV\f$
 *
 * \param[in] kinE    Variable where the kinetic energy will be stored
 * \param[in] n       The density
 * \param[in] phi     The potential
 * \param[in] uEPar   The parallel electron velocity
 * \param[in] uIPar   The parallel ion velocity
 *
 * \param[out] kinE   Variable where the kinetic energy is stored
 */
template<class PhysicsClass>
void ownMonitors<PhysicsClass>::kinEnergy(BoutReal &kinE      ,
                                          Field3D const &n    ,
                                          Field3D const &phi  ,
                                          Field3D const &uEPar,
                                          Field3D const &uIPar)
{
    TRACE("Halt in ownMonitors<PhysicsClass>::kinEnergy");
}

#endif
