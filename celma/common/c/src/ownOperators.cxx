#ifndef __OWNOPERATORS_CXX__
#define __OWNOPERATORS_CXX__

#include "../include/ownOperators.hxx"
#include <fft.hxx>                     //Includes the FFT
#include <interpolation.hxx>           //Includes the interpolation

/*!
 * \brief Constructor
 *
 * Constructor which sets the private member data
 *
 * \param[in] section What section to get bndryFuncGen from.
 *                    Default is "phi"
 */
OwnOperators::OwnOperators(const string &section) :
    incXBndry(false)
{
    TRACE("Halt in OwnOperators::OwnOperators");

    // Get the function and lastX value (used in D3DX3)
    // ************************************************************************
    // Get the function
    Options *varOptions = Options::getRoot()->getSection(section);
    string bndryFuncString;
    // Last argument in get is the default
    varOptions->get("bndry_xout", bndryFuncString, "");
    if (bndryFuncString == ""){
        output << "\n\n\n\n" << "!!!WARNING!!!: No bndry_xout found in section "
               << section << " D3DX3 will not work"
               << "\n\n\n\n" << std::endl;

    }
    // Strip the function name
    int pos = bndryFuncString.find('(');
    int pos2 = bndryFuncString.rfind(')');
    // Recycle bndryFuncString
    // This time the string is stripped
    bndryFuncString = trim(bndryFuncString.substr(pos+1, pos2-pos-1));
    bndryFuncGen = FieldFactory::get()->parse(bndryFuncString);

    // Get the last x value
    if(mesh->lastX()){
        /* NOTE:
         * For a local index, globalX returns the global value between 0
         * and 1 corresponding to that index.
         * globalLastXVal will in this case return 1.0 When evaluating
         * the function given in bndryFuncGen, the x variable need to be
         * in range 0-1 (even if it contains expressions like geom:xl)
         * in order to be consistent with the rest of the code.
         */
        globalLastXVal = 0.5 * (mesh->GlobalX(mesh->xend+1) + mesh->GlobalX(mesh->xend));
    }
    // ************************************************************************

    // Calculate the powers of the Jacobian
    // ************************************************************************
    J  = mesh->J;
    J2 = mesh->J^(2.0);
    J3 = mesh->J^(3.0);
    J4 = mesh->J^(4.0);
    // ************************************************************************

}

// Operators
/*!
 * 3rd derivative in the \f$\rho\f$ direction using cylindrical geometry. As
 * this is a wide stencil, extra care needs to be taken at the boundaries of a
 * processor.  A 3rd order scheme (using one guard cell) is used of the points
 * closest to the processor edges, and a 3rd order scheme **not** using the
 * last guard cell is used for xout. This avoid the need of having two ghost
 * points
 *
 * \param[in] f The original field
 * \param[in] t The time
 *
 * \return result The result of the operation
 *
 * \note t is an input as the boundary may be time dependent
 * \note A 3rd order stencil around the processor boundaries are chosen so that
 *       the error doesn't dominate there
 * \warning This is not made for staggered grids
 */
Field3D OwnOperators::D3DX3(const Field3D &f,
                            const BoutReal &t)
{
    TRACE("Halt in OwnOperators::D3DX3");
    // Need to initialize Field3D
    Field3D result = 0.0;
    BoutReal y, valBndry;
    int xend;

    /* The stencil for the inner points
     * mesh->xstart would need a second ghost point, so we start at
     * mesh->xstart+1
     * In similar manner mesh->xend would need a second ghost point, so we
     * iterate up to mesh->xend-1
     * If we are at the last x-processor, we will use a special case for the
     * two last point in order to avoid 1st order convergence, so in this case,
     * we only iterate up mesh->xend-2
     */
    if (mesh->lastX()){
        xend = mesh->xend-2;
    }
    else{
        xend = mesh->xend-1;
    }
    for(int xInd = mesh->xstart+1; xInd <= xend; xInd++){
        for(int yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
            for(int zInd = 0; zInd < mesh->ngz -1; zInd ++){
                result(xInd, yInd, zInd) =
                    (
                    - 0.5*f(xInd-2, yInd, zInd)
                    +     f(xInd-1, yInd, zInd)
                    -     f(xInd+1, yInd, zInd)
                    + 0.5*f(xInd+2, yInd, zInd)
                    )/pow(mesh->dx(xInd, yInd), 3.0)
                    ;
            }
        }
    }

    // Loop for the innermost xpoint at this processor (3rd order)
    int xInd = mesh->xstart;
    for(int yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
        for(int zInd = 0; zInd < mesh->ngz -1; zInd ++){
            // 3rd order for first point
            result(xInd, yInd, zInd) =
                (
                -  7.0*f(xInd-1, yInd, zInd)
                + 25.0*f(xInd  , yInd, zInd)
                - 34.0*f(xInd+1, yInd, zInd)
                + 22.0*f(xInd+2, yInd, zInd)
                -  7.0*f(xInd+3, yInd, zInd)
                +      f(xInd+4, yInd, zInd)
                )/(4.0*pow(mesh->dx(xInd, yInd), 3.0))
                ;
        }
    }

    // If we are on the last x-processor, we have a special case for the two
    // last x-points
    if(mesh->lastX()){
        // Last point
        int xInd = mesh->xend;
        // Calculate the x half-way between guard cell and grid cell
        for(int yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
            // Calculate the y for the current yIndex
            /* NOTE:
             * For a local index, globalY returns the global value between 0
             * and 1 corresponding to that index.
             * When evaluating the function given in bndryFuncGen, the y
             * variable need to be in range 0-2*pi (even if it contains
             * expressions like geom:yl) in order to be consistent with the
             * rest of the code.
             */
            y = mesh->GlobalY(yInd)*TWOPI;
            for(int zInd = 0; zInd < mesh->ngz -1; zInd ++){
                // Getting the on the boundary
                valBndry = bndryFuncGen->generate(globalLastXVal,
                                             y,
                                             TWOPI*zInd/(mesh->ngz-1),
                                             t);

                // 3rd order for second last point
                result(xInd-1, yInd, zInd) =
                    (
                    +  14.0*f(xInd-4, yInd, zInd)
                    -  99.0*f(xInd-3, yInd, zInd)
                    + 189.0*f(xInd-2, yInd, zInd)
                    - 105.0*f(xInd-1, yInd, zInd)
                    -  63.0*f(xInd  , yInd, zInd)
                    +  64.0*valBndry
                    )/(63.0*pow(mesh->dx(xInd, yInd), 3.0))
                    ;

                // 3rd order for last point
                result(xInd, yInd, zInd) =
                    (
                    -   4.0*f(xInd-4, yInd, zInd)
                    +  27.0*f(xInd-3, yInd, zInd)
                    -  81.0*f(xInd-2, yInd, zInd)
                    + 129.0*f(xInd-1, yInd, zInd)
                    - 135.0*f(xInd  , yInd, zInd)
                    +  64.0*valBndry
                    )/(9.0*pow(mesh->dx(xInd, yInd), 3.0))
                    ;
            }
        }
    }
    else{
        // We are not on the last x-processor, so we only need to take care of
        // the last point
        int xInd = mesh->xend;
        for(int yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
            for(int zInd = 0; zInd < mesh->ngz -1; zInd ++){
                // 3rd order for first point
                result(xInd, yInd, zInd) =
                    (
                    -      f(xInd-4, yInd, zInd)
                    +  7.0*f(xInd-3, yInd, zInd)
                    - 22.0*f(xInd-2, yInd, zInd)
                    + 34.0*f(xInd-1, yInd, zInd)
                    - 25.0*f(xInd  , yInd, zInd)
                    +  7.0*f(xInd+1, yInd, zInd)
                    )/(4.0*pow(mesh->dx(xInd, yInd), 3.0))
                    ;
            }
        }
    }

    return result;
}

/*!
 * 3rd derivative in the \f$\rho\f$ direction using cylindrical geometry.
 * This is a natural extension of D2DZ2 in src/sys/derivs.cxx, with the only
 * execption that incXBndry is a private member data set by the
 * constructor.
 *
 * \param[in] f The original field
 *
 * \return result The result of the operation
 */
Field3D OwnOperators::D3DZ3(const Field3D &f)
{
    TRACE("Halt in OwnOperators::D3DZ3");
    CELL_LOC inloc = f.getLocation(); // Input location
    CELL_LOC diffloc = inloc; // Location of differential result
    CELL_LOC outloc    = inloc;

    Field3D result;

    // Use FFT

    BoutReal shift = 0.; // Shifting result in Z?

    result.allocate(); // Make sure data allocated

    int ncz = mesh->ngz-1;

    #ifndef _OPENMP
    static dcomplex *cv = (dcomplex*) NULL;
    #else
    static dcomplex *globalcv;
    static int nthreads = 0;
    #endif

    #pragma omp parallel
    {
        #ifndef _OPENMP
        // Serial, so can have a single static array
        if(cv == (dcomplex*) NULL)
            cv = new dcomplex[ncz/2 + 1];
        #else
        // Parallel, so allocate a separate array for each thread

        int th_id = omp_get_thread_num(); // thread ID
        int n_th = omp_get_num_threads();
        if(th_id == 0) {
            if(nthreads < n_th) {
                // Allocate memory in thread zero
                if(nthreads > 0)
                    delete[] globalcv;
                globalcv = new dcomplex[n_th*(ncz/2 + 1)];
                nthreads = n_th;
            }
        }
        // Wait for memory to be allocated
        #pragma omp barrier

        dcomplex *cv = globalcv + th_id*(ncz/2 + 1); // Separate array for each thread
        #endif
        int xs = mesh->xstart;
        int xe = mesh->xend;
        int ys = mesh->ystart;
        int ye = mesh->yend;
        if(incXBndry) { // Include x boundary region (for mixed XZ derivatives)
            xs = 0;
            xe = mesh->ngx-1;
        }
        if (mesh->freeboundary_xin && mesh->firstX() && !mesh->periodicX)
            xs = 0;
        if (mesh->freeboundary_xout && mesh->lastX() && !mesh->periodicX)
            xe = mesh->ngx-1;
        if (mesh->freeboundary_ydown)
            ys = 0;
        if (mesh->freeboundary_yup)
            ye = mesh->ngy-1;
        #pragma omp for
        for(int jx=xs;jx<=xe;jx++) {
            for(int jy=ys;jy<=ye;jy++) {
                rfft(f[jx][jy], ncz, cv); // Forward FFT

                for(int jz=0;jz<=ncz/2;jz++) {
                    BoutReal kwave=jz*2.0*PI/mesh->zlength(); // wave number is 1/[rad]

                    BoutReal flt;
                    if (jz>0.4*ncz) flt=1e-10; else flt=1.0;
                    cv[jz] *= dcomplex(0.0, - pow(kwave, 3.0)) * flt;
                            -SQ(kwave) * flt;
                    if(mesh->StaggerGrids)
                        cv[jz] *= exp(Im * (shift * kwave * mesh->dz));
                }

                irfft(cv, ncz, result[jx][jy]); // Reverse FFT

                result[jx][jy][ncz] = result[jx][jy][0];
            }
        }
    }
    // End of parallel section

    #ifdef CHECK
    // Mark boundaries as invalid
    if (mesh->freeboundary_xin) result.bndry_xin = true;
    else result.bndry_xin = false;
    if (mesh->freeboundary_xout) result.bndry_xout = true;
    else result.bndry_xout = false;
    if (mesh->freeboundary_yup) result.bndry_yup = true;
    else result.bndry_yup = false;
    if (mesh->freeboundary_ydown) result.bndry_ydown = true;
    else result.bndry_ydown = false;
    #endif

    result.setLocation(diffloc);

    return interp_to(result, outloc);
}

/*!
 * Operator for \f$\nabla\cdot_(f \nabla_\perp g)\f$ in cylindrical geometry.
 * We have that
 *
 * \f{eqnarray}{
 *   \nabla\cdot(f\nabla_\perp g) = f\nabla_\perp^2g + \nabla
 *   f\cdot \nabla_\perp g = f\nabla_\perp^2g + \nabla_\perp f\cdot
 *   \nabla_\perp g
 * \f}
 *
 * The expression for the perpendicular Laplacian can be found in the
 * coordinates manual. Note that in cylinder coordinates
 *
 * \f{eqnarray}{
 *   G^x &=& \frac{1}{J}\\
 *   G^y &=& 0\\
 *   G^z &=& 0\\
 *   g^{zz} &=& \frac{1}{\rho^2}\\
 *   g^{yy}\partial_y^2
 *   - \frac{1}{J}\partial_y\left(\frac{J}{g^{yy}}\partial_y\right)
 *   &=& 0
 * \f}
 *
 * \param[in] f The f field
 * \param[in] g The g field
 *
 * \return result The result of the operation
 */
Field3D OwnOperators::div_f_GradPerp_g(const Field3D &f,
                                       const Field3D &g)
{
    TRACE("Halt in OwnOperators::div_f_GradPerp_g");

    Field3D result;

    result =   f*D2DX2(g)
             + (f/J)*DDX(g)
             + (f/J2)*D2DZ2(g)
             + DDX(f)*DDX(g)
             + (1/J2)*DDZ(f)*DDZ(g)
             ;

    return result;
}

/*!
 * Operator for \f$\nabla\cdot_(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$ using
 * cylindrical geometry. The derivation can be found in the derivation folder.
 *
 * \param[in] n The density
 * \param[in] phi The potential
 *
 * \return result The result of the operation
 */
Field3D OwnOperators::div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                                     const Field3D &phi)
{
    TRACE("Halt in OwnOperators::div_uE_dot_grad_n_GradPerp_phi");

    Field3D result;

    phi_x   = DDX      (phi)       ;
    phi_xx  = D2DX2    (phi)       ;
    phi_xxx = D3DX3    (phi)       ;
    phi_z   = DDZ      (phi)       ;
    phi_zz  = D2DZ2    (phi)       ;
    phi_zzz = D3DZ3    (phi)       ;
    phi_xz  = DDX(DDZ  (phi, true));
    phi_xxz = D2DX2(DDZ(phi, true));
    phi_xzz = DDX(D2DZ2(phi, true));
    n_x     = DDX      (n)         ;
    n_xx    = D2DX2    (n)         ;
    n_z     = DDZ      (n)         ;
    n_zz    = D2DZ2    (n)         ;
    n_xz    = DDX(DDZ  (n, true))  ;

    result = (
                  n*(
                    - phi_x*(
                          phi_xxz*J3
                        + phi_xz*J2
                        + phi_z*J
                        + phi_zzz*J)
                     + phi_z*(
                           phi_xx*J2
                         + phi_xxx*J3
                         + phi_xzz*J
                         - 2.0*phi_zz)
                     )
                - (phi_x^(2.0))*(
                      n_xz*J3
                    + n_z*J2
                                  )
                + (phi_z^(2.0))*(
                      n_xz*J
                    - n_z
                                  )
                + phi_x*(
                      phi_z*(
                          n_x*J2
                        + n_xx*J3
                        - n_zz*J
                            )
                    - 2.0*n_z*(
                          phi_xx*J3
                        + phi_zz*J
                              )
                        )
                + 2.0*n_x*phi_z*(
                          phi_xx*J3
                        + phi_zz*J
                                )
             )/J4;

    return result;
}

/*!
 * \f$\nabla_\perp f\f$, equivalent to Grad_perp in vecops.cxx, but in
 * cylindrical geometry.  This means that there the y-component of the vector
 * is set to 0 (as there are no off-diagonal elements). The advantage of
 * introducing this operator is:
 *
 *      1. No need for setting the y-boundaries on the operand
 *      2. Reduced calculation time
 *
 * \param[in] f The original field
 *
 * \return result The result of the operation written in a covariant basis
 */
Vector3D OwnOperators::Grad_perp(const Field3D &f)
{
    TRACE("Halt in OwnOperators::own_Grad_perp");
    Vector3D result;

    result.x = DDX(f);
    result.y = 0.0;
    result.z = DDZ(f);

    result.covariant = true;

    return result;
}

// Auxiliary functions
/*!
 * Getter for incXBndry
 *
 * \return incXBndry Used to determine if derivatives of the ghost cells are
 *                   going to be taken
 */
bool OwnOperators::getIncXbndry(){
    return incXBndry;
}

/*!
 * Setter for incXBndry
 *
 * \param[in] option The state of incXBndry
 * \param[out] incXBndry Used to determine if derivatives of the ghost cells are
 *                       going to be taken
 */
void OwnOperators::setIncXbndry(bool option){
    incXBndry = option;
}
#endif
