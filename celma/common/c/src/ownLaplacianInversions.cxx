#ifndef __OWNLAPLACIANINVERSIONS_CXX__
#define __OWNLAPLACIANINVERSIONS_CXX__

#include "../include/ownLaplacianInversions.hxx"

/*!
 * This function is used instead of a constructor as a OwnOperators and a
 * OwnBCs object is needed as input. The function sets the private member data.
 * Options belonging to the phiSolver is read from the input file.
 *
 * \param[in] opObj Object containing own operators
 * \param[in] BCObj Object containing own boundary conditions
 * \param[in] section Section to read the input data from
 */
void OwnLaplacianInversions::create(OwnOperators &opObj,
                                    OwnBCs &BCObj,
                                    const string &section)
{

    TRACE("Halt in OwnLaplacianInversions::OwnLaplacianInversions");

    // Set the operator and boundary object
    ownOp = &opObj;
    ownBC = &BCObj;

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Create the laplace object
    // ************************************************************************
    // The laplace object will look in the section phiSolver in the BOUT.inp
    // file
    Options *phiSol_opt = options->getSection(section);
    phiSolver = Laplacian::create(phiSol_opt);
    // Set the coefficients manually (should also be set to this by default)
    phiSolver->setCoefD(1.0);
    phiSolver->setCoefC(1.0);
    phiSolver->setCoefA(0.0);
    // Get the tolerances
    phiSol_opt->get("atol",  atol,  1.0e-10);
    phiSol_opt->get("rtol",  rtol,  1.0e-5);
    phiSol_opt->get("maxit", maxit, 300);
    // Check if the Naulin iterator should be monitored
    phiSol_opt->get("monitor", monitor, false);
    // ************************************************************************
}

// Member function
/*!
 * Implementation of the Naulin solver.
 *
 * \param[in] gradPerpLnN \f$\nabla_\perp \ln(n) \f$
 * \param[in] n           The density
 * \param[in] vortD       The current \f$\nabla\cdot(n\nabla_\perp \phi)\f$
 * \param[in] phiInit     The initial \f$\phi\f$
 *
 * \param[out] vort The vorticity
 *
 * \return phiNext A phi which matches the tolerances
 *
 * ## Explanation of the procedure:
 * An way to invert the equation \f$\Omega^D = \nabla\cdot(n\nabla_\perp \phi)\f$
 * invented by Naulin, V.
 * In an orthogonal system, we have that:
 *
 * \f{eqnarray}{
 * \Omega^D &=& \nabla\cdot(n\nabla_\perp \phi)\\
 *       &=& n \nabla_\perp^2 \phi + \nabla n\cdot\nabla_\perp \phi\\
 *       &=& n \Omega + \nabla n\cdot\nabla_\perp \phi\\
 *       &=& n \Omega + \nabla_\perp n\cdot\nabla_\perp \phi
 * \f}
 *
 * Rearranging gives
 *
 * \f{eqnarray}{
 * \Omega  &=& \frac{\Omega^D}{n} - \nabla_\perp \ln(n)\cdot\nabla_\perp \phi\\
 * \nabla_\perp^2 \phi
 * &=& \frac{\Omega^D}{n} - \nabla_\perp \ln(n)\cdot\nabla_\perp \phi
 * \f}
 *
 * The iteration now works as follows:
 *      1. Get the voritcity from
 *         \code{.cpp}
 *         vort = (vortD/n) - grad_perp(ln_n)*grad_perp(phiCur)
 *         \endcode
 *         where phiCur is phi of the current iteration
 *      2. Invert \f$phi\f$ to find the voricity using
 *         \code{.cpp}
 *         phiNext = invert_laplace_perp(vort)
 *         \endcode
 *         where phiNext is the newly obtained \f$phi\f$
 *      3. Calculate
 *         \code{.cpp}
 *         EAbsLInf = phiCur - phi_new
 *         ERelLInf = (phiCur - phi_new)/phiCur
 *         \endcode
 *      4. Calculate the infinity norms of the Es
 *         \code{.cpp}
 *         EAbsLInf = max(abs(EAbsLInf))
 *         ERelLInf = max(abs(ERelLInf))
 *         \endcode
 *      5. Check whether
 *         \code{.cpp}
 *         EAbsLInf > atol
 *         \endcode
 *          * If yes
 *              * Check whether
 *                \code{.cpp}
 *                ERelLInf > rtol
 *                \endcode
 *              * If yes
 *                  * Set
 *                    \code{.cpp}
 *                    phiCur = phiNext
 *                    \endcode
 *                    increase curCount and start from step 1
 *                  * If number of iteration is above maxit, trhow exception
 *              * If no
 *                  * Stop: Function returns
 *          * if no
 *              * Stop: Function returns
 *
 * \note The vorticity is an output together with the return value (phi)
 * \note Needs orthogonal basis
 */
Field3D OwnLaplacianInversions::NaulinSolver(const Vector3D &gradPerpLnN,
                                             const Field3D  &n,
                                             const Field3D  &vortD,
                                             const Field3D  &phiInit,
                                                   Field3D  &vort)
{
    TRACE("Halt in OwnLaplacianInversions::NaulinSolver");

    curCount = 0;             // Resetting the counter
    phiCur   = copy(phiInit); // Init phiCur (no shared memory due to copy)

    // Infinite loop
    for(;;){
        /* For this to work, the following fields must be communicated:
         * 1. phiCur
         */
        mesh->communicate(phiCur);
        vort = (vortD/n) - gradPerpLnN*ownOp->Grad_perp(phiCur);
        // Solve for phi (solve takes care of the perp boundaries)
        /* No need for a second argument here because:
         * 1. This is a direct solver (no need for initial guess)
         * 2. Boundaries are set through a flag
         */
        /* For this to work, the following fields must be communicated:
         * 1. vort
         */
        mesh->communicate(vort);
        phiNext = phiSolver->solve(vort);

        // Calculate the errors (the true in the end means the max will be
        // taken over all processors)
        EAbsLInf = (abs( phiCur - phiNext)        ).max(true);
        ERelLInf = (abs((phiCur - phiNext)/phiCur)).max(true);

        if(EAbsLInf > atol){
            if(ERelLInf > rtol){
                // Solve for phi has taken care of perpendicular BC
                // Just to be sure, we force the boundary on inner rho
                ownBC->innerRhoCylinder(phiNext);
                // Update phiCur (no shared memory due to copy)
                phiCur = copy(phiNext);
            }
            else{
                break; // Condition fulfilled, break out of the loop
            }
        }
        else{
            break; // Condition fulfilled, break out of the loop
        }

        // Increase the counter
        ++curCount;

        if (curCount > maxit){
            throw BoutException("NaulinSolver reached max iterations with\n"
                                "iterations=%d, atol=%.3e and rtol=%.3e\n"
                                "aerr=%.3e, rerr=%.3e",
                                curCount, atol, rtol, EAbsLInf, ERelLInf);
        }
    }

    if(monitor){
        output << "NaulinSolver success:"<<
                  " Iterations=" << curCount <<
                  " abserr=" << EAbsLInf <<
                  " relerr=" << ERelLInf <<
                  std::endl;
    }

    // Return phiNext, vort is also changed, as it is a call by
    // reference
    return phiNext;
 }

// Destructor
/*!
 * Destroys the laplace object
 */
OwnLaplacianInversions::~OwnLaplacianInversions()
{
    TRACE("Halt in OwnLaplacianInversions::~OwnLaplacianInversions");

    delete phiSolver;
}
#endif
