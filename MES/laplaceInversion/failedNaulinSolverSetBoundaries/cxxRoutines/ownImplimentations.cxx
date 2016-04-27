#include "../NaulinSolver.hxx"

// Own laplacian solver
// ############################################################################
 Field3D NaulinSolver::NSolver(const BoutReal &atol,
                               const BoutReal &rtol,
                               const int &maxit,
                               const Field3D &ln_n,
                               const Field3D &n,
                               const Field3D &vortD,
                               Field3D &phi,
                               Field3D &vort,
                               bool const &monitor){
    /* Explanation:
     * An way to invert the equation vortD = div(n(grad_perp(phi)) invented by
     * Naulin, v.
     * We have that:
     *
     * vortD = div(n(grad_perp(phi))
     *       = n laplace_perp(phi) + grad(n)*grad_perp(phi)
     *       = n vort + grad(n)*grad_perp(phi)
     *       = n vort + grad_perp(n)*grad_perp(phi)
     *
     * The last equality comes from the b field being perpendicular on the
     * other basis vectors. This gives:
     *
     * vort = (vortD/n) + grad_perp(ln_n)*grad_perp(phi)
     * laplacian_perp(phi) = (vortD/n) + grad_perp(ln_n)*grad_perp(phi)
     *
     * The iteration now works as follows:
     * Step 1   - Get the voritcity from
     *            vort = (vortD/n) + grad_perp(ln_n)*grad_perp(phiCur)
     *            where phiCur is phi of the current iteration
     * Step 2   - Invert phi to find the voricity using
     *            vort = laplace_perp(phi_new)
     *            where phi_new is the newly obtained phi
     * Step 3.1 - Calculate
     *            E_abs = phiCur - phi_new
     *            E_rel = (phiCur - phi_new)/phiCur
     * Step 3.2 - Calculate the infinity norms of the Es
     *            E_abs_Linf = max(abs(E_abs))
     *            E_rel_Linf = max(abs(E_rel))
     * Step 4   - Check whether E_abs_Linf > atol
     *            if yes
     *              Check whether E_rel_Linf > rtol
     *              if yes, set phiCur = phi_new and start from step 1
     *              if no stop (phiCur, vort_cur and vortD_cur is returned)
     *            if no stop (phiCur, vort_cur and vortD_cur is returned)
     * Step 5   - If number of iteration is above maxit, trhow exception
     *
     *
     *
     * Input:
     * atol    - absolute tolerance
     * rtol    - relative tolerance
     * ln_n    - logarithm of the density
     * n       - the density
     * vortD   - the density-vorticity
     * phi     - the potential
     * vort    - the vorticity
     *
     * Output:
     * phi     - the potential
     * vort    - the vorticity
     */
    TRACE("Halt in NaulinSolver::NSolver");

    // A new object is created, these objects do not share memory
    Field3D phiCur   = phi;
    int curCount     = 0;   // Make the counter

    // Communicate the current phi
    mesh->communicate(phiCur);

    // Infinite loop
    for(;;){
        vort = (vortD/n) - own_Grad_perp(ln_n)*own_Grad_perp(phiCur);
        // Solve for phi (solve takes care of the perp boundaries)
        /* Explanation of solve
         * The second argument has one purpose in the fourier solver
         * 1. Sets the boundary condition if inner/outer flag 16 - "INVERT_SET is set"
         */
        phi = phiSolver->solve(vort, phi);

        // Calculate the errors (the true in the end means the max will be
        // taken over all processors)
        E_abs_Linf = (abs( phiCur - phi)        ).max(true);
        E_rel_Linf = (abs((phiCur - phi)/phiCur)).max(true);

        if(E_abs_Linf > atol){
            if(E_rel_Linf > rtol){
                // Force the boundary on outer and inner rho
                phi.applyBoundary();
                innerRhoCylinder(phi);
                phiCur = phi;            // Update phiCur
                // Communication is needed as we take derivatives
                mesh->communicate(phi, phiCur);
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
                                curCount, atol, rtol, E_abs_Linf, E_rel_Linf);
        }
    }

    if(monitor){
        output << "NaulinSolver success:"<<
                  " Iterations=" << curCount <<
                  " abserr=" << E_abs_Linf <<
                  " relerr=" << E_rel_Linf <<
                  std::endl;
    }

    // Force the boundary on outer and inner rho
    phi.applyBoundary();
    innerRhoCylinder(phi);
    phiCur = phi;            // Update phiCur

    return phiCur;
 }
// ############################################################################

// Own operators
// ############################################################################
// Own Grad_perp operator
// ****************************************************************************
const Vector3D NaulinSolver::own_Grad_perp(const Field3D &f) {
    /* Info:
     * Equivalent to Grad_perp in vecops.cxx, but without the non-diagonal
     * elements in the metric tensor
     */
    TRACE("Halt in NaulinSolver::own_Grad_perp");
    Vector3D result;

    result.x = DDX(f);
    result.y = 0.0;
    result.z = DDZ(f);

    result.covariant = true;

    return result;
}
// ****************************************************************************
// ############################################################################
