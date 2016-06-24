#ifndef __OWNOPERATORS_CXX__
#define __OWNOPERATORS_CXX__

#include "../include/ownOperators.hxx"
#include <fft.hxx>                     //Includes the FFT
#include <interpolation.hxx>           //Includes the interpolation

// OwnOperators

/*!
 * \brief Constructor
 *
 * Constructor which sets the private member data
 *
 * \param[in] section What section to get bndryFuncGen from.
 *                    Default is "phi"
 */
OwnOperators::OwnOperators(Options *option) :
    incXBndry(false)
{
    TRACE("Halt in OwnOperators::OwnOperators");

    // Calculate the powers of the Jacobian
    // ************************************************************************
    J  = mesh->J;
    J2 = mesh->J^(2.0);
    J3 = mesh->J^(3.0);

    invJ  = 1.0/(mesh->J      );
    invJ2 = 1.0/(mesh->J^(2.0));
    invJ3 = 1.0/(mesh->J^(3.0));
    invJ4 = 1.0/(mesh->J^(4.0));
    // ************************************************************************

    /* Check that there are enough points
     * ngx is the size of local mesh including guard cells
     */
    Options *switchOptions = Options::getRoot()->getSection("switch");
    switchOptions->get("warnPoints", warnPoints, false);
    if (mesh->ngx - 2*mesh->xstart < 4){

        // Create a stream which we cast to a string
        std::ostringstream stream;
        stream << "Not enough inner points i the x-direction\n"
               << "D3DX3 needs 4 inner points  x\n"
               << "Currently the number of inner points is "
               << mesh->ngx- 2*mesh->xstart;

        if (warnPoints){
            output << "\n\n!!! WARNING" << stream << "\n\n" << std::endl;
        }
        else{
            std::string str =  stream.str();
            // Cast the stream to a const char in order to use it in BoutException
            const char* message = str.c_str();
            throw BoutException(message);
        }
    }
}

/*!
 * This function now works as a constructor of the child-classes of OwnFilters
 */
OwnOperators* OwnOperators::createOperators(Options *options)
{
    TRACE("Halt in OwnOperators::createOperators");

    // The filter option is by defualt found in the ownFilter section
    if(options == NULL){
        options = Options::getRoot()->getSection("ownOperators");
    }

    string type;
    options->get("type", type, "simpleStupid");

    if(lowercase(type) == lowercase("simpleStupid")){
        string phiBndrySec;
        options->get("phiBndrySec", phiBndrySec, "phi");
        output << "OwnOperators type set to 'simpleStupid'" << std::endl;
        output << "Will look for phi boundaries in section '" << phiBndrySec
               << "'"<< std::endl;
        return new OwnOpSimpleStupid(options, phiBndrySec);
    }
    else if(lowercase(type) == lowercase("onlyBracket")){
        output << "OwnOperators type set to 'onlyBracket'" << std::endl;
        return new OwnOpOnlyBracket(options);
    }
    else if(lowercase(type) == lowercase("2Brackets")){
        output << "OwnOperators type set to '2Brackets'" << std::endl;
        return new OwnOp2Brackets(options);
    }
    else if(lowercase(type) == lowercase("3Brackets")){
        output << "OwnOperators type set to '3Brackets'" << std::endl;
        return new OwnOp3Brackets(options);
    }
    else if(lowercase(type) == lowercase("3BasicBrackets")){
        output << "OwnOperators type set to '3BasicBrackets'" << std::endl;
        return new OwnOp3BasicBrackets(options);
    }
    else {
        // Create a stream which we cast to a string
        std::ostringstream stream;
        stream << "OwnOperators '"<< type << "' not implemented\n"
               << "Available operators:\n"
               << "simpleStupid   - Simple stupid implementation of "
                                    "div_uE_dot_grad_n_GradPerp_phi\n"
               << "onlyBracket    - Only {phi, vortD} will be used. "
                                    "NOTE: Not consistent\n"
               << "2Brackets      - Consistent implementation using 2 brackets.\n"
               << "3Brackets      - Consistent implementation using 3 brackets "
                                    "with modified Arakawa.\n"
               << "3BasicBrackets - Consistent implementation using 3 brackets "
                                    "with basic Arakawa.\n"
               ;
        std::string str =  stream.str();
        // Cast the stream to a const char in order to use it in BoutException
        const char* message = str.c_str();

        throw BoutException(message);
    }
}

/*!
 * Only implemented in the simpleStupid class
 *
 * \param[in] f The field
 * \param[in] t The time to evaluate the boundary function
 *
 * \returns f (never reached)
 */
Field3D OwnOperators::D3DX3(const Field3D &f, const BoutReal &t)
{
    TRACE("Halt in OwnOperators::D3DX3");

    throw BoutException("D3DX3 is only implemented in the simpleStupid class");

    return f;
}

/*!
 * Only implemented in the 3Brackets class
 *
 * \param[in] phi The potential (not (DDX(phi)^2.0))
 * \param[in] n   The density
 *
 * \returns phi (never reached)
 */
Field3D OwnOperators::ArakawaOfDDXPhi2N(const Field3D &phi, const Field3D &n)
{
    TRACE("Halt in OwnOperators::ArakawaOfDDXf2g");

    throw BoutException("D3DX3 is only implemented in the 3Brackets class");

    return phi;
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

/*!
 * Getter for incXBndry
 *
 * \return incXBndry Used to determine if derivatives of the ghost cells are
 *                   going to be taken
 */
bool OwnOperators::getIncXbndry()
{
    return incXBndry;
}

/*!
 * Setter for incXBndry
 *
 * \param[in] option The state of incXBndry
 * \param[out] incXBndry Used to determine if derivatives of the ghost cells are
 *                       going to be taken
 */
void OwnOperators::setIncXbndry(bool option)
{
    incXBndry = option;
}

// OwnOpSimpleStupid

/*!
 * \brief Constructor
 *
 * Constructor which sets the private member data
 *
 * \param[in] section What section to get bndryFuncGen from.
 *                    Default is "phi"
 *
 * \warning Only dirichlet boundary condition is implemented.
 */
OwnOpSimpleStupid::OwnOpSimpleStupid(Options *options, string &phiBndrySec) :
    OwnOperators(options)
{
    TRACE("Halt in OwnOpSimpleStupid::OwnOpSimpleStupid");

    // Get the function and lastX value (used in D3DX3)
    // ************************************************************************
    // Get the function
    Options *varOptions = Options::getRoot()->getSection(phiBndrySec);
    string bndryFuncString;
    // Last argument in get is the default
    varOptions->get("bndry_xout", bndryFuncString, "");
    if (bndryFuncString == ""){
        output << "\n\n\n\n" << "!!!WARNING!!!: No bndry_xout found in section "
               << phiBndrySec << " D3DX3 will not work"
               << "\n\n\n\n" << std::endl;

    }
    else if (!(lowercase(bndryFuncString).find("dirichlet") != string::npos)){
        // Create a stream which we cast to a string
        std::ostringstream stream;
        stream << "neumann boundary condition not implemented in D3DX3\n"
               << "Found:\n"
               << bndryFuncString
               << "\nin\n"
               << phiBndrySec
               << "\nHowever, current implementation only includes dirichlet BC"
               ;
        std::string str =  stream.str();
        // Cast the stream to a const char in order to use it in BoutException
        const char* message = str.c_str();

        throw BoutException(message);
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
}

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
Field3D OwnOpSimpleStupid::D3DX3(const Field3D &f,
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
 * Operator for \f$\nabla\cdot_(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$ using
 * cylindrical geometry. The derivation can be found in the derivation folder.
 *
 * \param[in] n The density
 * \param[in] phi The potential
 * \param[in] vortD The modified vorticity (not used here)
 *
 * \return result The result of the operation
 */
Field3D OwnOpSimpleStupid::div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                                          const Field3D &phi)
{
    TRACE("Halt in OwnOpSimpleStupid::div_uE_dot_grad_n_GradPerp_phi");

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
             )*invJ4;

    return result;
}

/*!
 * Not implemented in this child class
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns phi (never reached)
 */
Field3D OwnOpSimpleStupid::vortDAdv(const Field3D &phi, const Field3D &n)
{
    TRACE("Halt in OwnOpSimpleStupid::vortDAdv");

    throw BoutException("vortDAdv not used in the OwnOpSimpleStupid "
                        "implementation");

    return phi;
}

/*!
 * Not implemented in this child class
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns phi (never reached)
 */
Field3D OwnOpSimpleStupid::kinEnAdvN(const Field3D &phi, const Field3D &n)
{
    TRACE("Halt in OwnOpSimpleStupid::kinEnAdvN");

    throw BoutException("kinEnAdvN not used in the OwnOpSimpleStupid "
                        "implementation");
}

// OwnOpOnlyBracket

/*!
 * \brief Constructor
 *
 * Constructor which calls parent constructor and sets the bracket method
 */
OwnOpOnlyBracket::OwnOpOnlyBracket(Options *options) :
    OwnOperators(options)
{
    TRACE("Halt in OwnOpOnlyBracket::OwnOpOnlyBracket");

    bm = BRACKET_ARAKAWA;
}

/*!
 * Not implemented in this child class
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns phi (never reached)
 */
Field3D OwnOpOnlyBracket::div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                                         const Field3D &phi)
{
    TRACE("Halt in OwnOpSimpleStupid::div_uE_dot_grad_n_GradPerp_phi");

    throw BoutException("div_uE_dot_grad_n_GradPerp_phi not used in the "
                        "OnlyBracket implementation");

    return phi;
}

/*!
 * Calculates \f$\{\phi, Omega^D\}\f$
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns result The result of the operation
 */
Field3D OwnOpOnlyBracket::vortDAdv(const Field3D &phi, const Field3D &vortD)
{
    TRACE("Halt in OwnOpOnlyBracket::vortDAdv");

    return invJ*bracket(phi, vortD, bm);
}

/*!
 * Not implemented in this child class
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns phi (never reached)
 */
Field3D OwnOpOnlyBracket::kinEnAdvN(const Field3D &phi, const Field3D &n)
{
    TRACE("Halt in OwnOpOnlyBracket::kinEnAdvN");

    throw BoutException("kinEnAdvN not used in the OnlyBracket "
                        "implementation");
}

// OwnOp2Brackets

/*!
 * \brief Constructor
 *
 * Constructor which calls parent constructor and sets the bracket method
 */
OwnOp2Brackets::OwnOp2Brackets(Options *options) :
    OwnOperators(options)
{
    TRACE("Halt in OwnOp2Brackets::OwnOp2Brackets");

    bm = BRACKET_ARAKAWA;
}

/*!
 * Not implemented in this child class
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns phi (never reached)
 */
Field3D OwnOp2Brackets::div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                                       const Field3D &phi)
{
    TRACE("Halt in OwnOp2Brackets::div_uE_dot_grad_n_GradPerp_phi");

    throw BoutException("div_uE_dot_grad_n_GradPerp_phi not used in the "
                        "OwnOp2Brackets implementation");

    return phi;
}

/*!
 * Calculates \f$\{\phi, \Omega^D\}\f$
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns result The result of the operation
 */
Field3D OwnOp2Brackets::vortDAdv(const Field3D &phi, const Field3D &vortD)
{
    TRACE("Halt in OwnOp2Brackets::vortDAdv");

    return invJ*bracket(phi, vortD, bm);
}

/*!
 * Calculates \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns result The result of the operation
 */
Field3D OwnOp2Brackets::kinEnAdvN(const Field3D &phi, const Field3D &n)
{
    TRACE("Halt in OwnOp2Brackets::kinEnAdvN");

    Field3D result;

    result =   invJ *DDX(phi)*(D2DXDZ(phi)*DDX(n) - D2DX2(phi)*DDZ(n))
             + invJ3*DDZ(phi)*bracket(DDZ(phi, true), n, bm)
             + invJ4*DDZ(n)  *((DDZ(phi))^(2.0))
           ;

    return result;
}

// OwnOp3Brackets

/*!
 * \brief Constructor
 *
 * Constructor which calls parent constructor and sets the bracket method
 */
OwnOp3Brackets::OwnOp3Brackets(Options *options) :
    OwnOperators(options)
{
    TRACE("Halt in OwnOp3Brackets::OwnOp3Brackets");

    bm = BRACKET_ARAKAWA;
}

/*!
 * Not implemented in this child class
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns phi (never reached)
 */
Field3D OwnOp3Brackets::div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                                       const Field3D &phi)
{
    TRACE("Halt in OwnOp3Brackets::div_uE_dot_grad_n_GradPerp_phi");

    throw BoutException("div_uE_dot_grad_n_GradPerp_phi not used in the "
                        "OwnOp3Brackets implementation");

    return phi;
}

/*!
 * Calculates \f$\{\phi, \Omega^D\}\f$
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns result The result of the operation
 */
Field3D OwnOp3Brackets::vortDAdv(const Field3D &phi, const Field3D &vortD)
{
    TRACE("Halt in OwnOp3Brackets::vortDAdv");

    return invJ*bracket(phi, vortD, bm);
}

/*!
 * Calculates \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns result The result of the operation
 */
Field3D OwnOp3Brackets::kinEnAdvN(const Field3D &phi, const Field3D &n)
{
    TRACE("Halt in OwnOp3Brackets::kinEnAdvN");

    Field3D result;

    result =   ArakawaOfDDXPhi2N(phi, n)
             + bracket(((invJ*DDZ(phi, true)^(2.0))), n, bm)
           ;

    // Multiply with B/2
    return 0.5*invJ*result;
}

/*!
 * Calculates \f$\{(\partial_\rho\phi)^2, n\}\f$ using a centered difference
 * stencil for \f$(\partial_\rho\phi)^2\f$ in the center, and a third order
 * one-sided stencils close to processor or domain boundaries.
 *
 * \param[in] phi The potential (not (DDX(phi)^2.0))
 * \param[in] n   The density
 *
 * \returns result The result of the operation
 */
Field3D OwnOp3Brackets::ArakawaOfDDXPhi2N(const Field3D &phi, const Field3D &n)
{
    TRACE("Halt in OwnOp3Brackets::ArakawaOfDDXf2g");

    Field3D result;

    result.allocate();

    ncz = mesh->ngz - 1;
    /* Loop over all inner points:
     * The Arakawa bracket requires evaluation of the field in x-1 and x+1. In
     * our case, we are at the same time taking the derivative of the fields,
     * which means that we need to calculate the derivatives in x-1 and x+1
     * (i.e. the ghost points/guard cells). Hence, special treatment is needed
     * for the points closest to the domain/processor boundary
     */
    xstart = mesh->xstart + 1;
    xend   = mesh->xend   - 1;
    for(xInd=xstart; xInd<=xend; xInd++){
        xIndP1 = xInd + 1;
        xIndM1 = xInd - 1;
        for(int yInd=mesh->ystart; yInd<=mesh->yend; yInd++){
            for(int zInd=0;zInd<ncz;zInd++) {
                zIndP1 = (zInd + 1) % ncz;
                zIndM1 = (zInd - 1 + ncz) % ncz;

                /* Boilerplate
                 * phi_xyz = pow( (-phi(xInd  -1, yInd, zInd  ) + phi(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                 */
                phi_xyzP1   = pow( (-phi(xInd  -1, yInd, zIndP1) + phi(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                phi_xyzM1   = pow( (-phi(xInd  -1, yInd, zIndM1) + phi(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                phi_xP1yz   = pow( (-phi(xIndP1-1, yInd, zInd  ) + phi(xIndP1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                phi_xM1yz   = pow( (-phi(xIndM1-1, yInd, zInd  ) + phi(xIndM1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                phi_xP1yzP1 = pow( (-phi(xIndP1-1, yInd, zIndP1) + phi(xIndP1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                phi_xP1yzM1 = pow( (-phi(xIndP1-1, yInd, zIndM1) + phi(xIndP1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                phi_xM1yzP1 = pow( (-phi(xIndM1-1, yInd, zIndP1) + phi(xIndM1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                phi_xM1yzM1 = pow( (-phi(xIndM1-1, yInd, zIndM1) + phi(xIndM1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);

                // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
                Jpp = 0.25*( (phi_xyzP1               - phi_xyzM1               )*
                             (n(xInd+1, yInd, zInd  ) - n(xInd-1, yInd, zInd  ) )
                             -
                             (phi_xP1yz               - phi_xM1yz                )*
                             (n(xInd  , yInd, zIndP1) - n(xInd  , yInd, zIndM1)) )
                  / (mesh->dx(xInd,yInd) * mesh->dz);

                // J+x
                Jpx = 0.25*( n(xInd+1,yInd, zInd  )*(phi_xP1yzP1 - phi_xP1yzM1) -
                             n(xInd-1,yInd, zInd  )*(phi_xM1yzP1 - phi_xM1yzM1) -
                             n(xInd  ,yInd, zIndP1)*(phi_xP1yzP1 - phi_xM1yzP1) +
                             n(xInd  ,yInd, zIndM1)*(phi_xP1yzM1 - phi_xM1yzM1))
                  / (mesh->dx(xInd,yInd) * mesh->dz);

                // Jx+
                Jxp = 0.25*( n(xInd+1, yInd, zIndP1)*(phi_xyzP1 - phi_xP1yz) -
                             n(xInd-1, yInd, zIndM1)*(phi_xM1yz - phi_xyzM1) -
                             n(xInd-1, yInd, zIndP1)*(phi_xyzP1 - phi_xM1yz) +
                             n(xInd+1, yInd, zIndM1)*(phi_xP1yz - phi_xyzM1))
                  / (mesh->dx(xInd,yInd) * mesh->dz);


                result(xInd,yInd,zInd) = (Jpp + Jpx + Jxp) / 3.;
            }
        }
    }

    // Treatment of the first point
    // That means that xM1 is a ghost/guard point
    xInd = mesh->xstart;
    xIndP1 = xInd + 1;
    xIndM1 = xInd - 1;
    for(int yInd=mesh->ystart; yInd<=mesh->yend; yInd++){
        for(int zInd=0; zInd<ncz; zInd++) {
            zIndP1 = (zInd + 1) % ncz;
            zIndM1 = (zInd - 1 + ncz) % ncz;

            // Unsure why this needs to be one-sided, but is needed for convergence
            phi_xyzP1   = pow(
                            (
                             - (2.0)*phi(xInd  -1, yInd, zIndP1)
                             - (3.0)*phi(xInd    , yInd, zIndP1)
                             + (6.0)*phi(xInd  +1, yInd, zIndP1)
                             -       phi(xInd  +2, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xyzM1   = pow(
                            (
                             - (2.0)*phi(xInd  -1, yInd, zIndM1)
                             - (3.0)*phi(xInd    , yInd, zIndM1)
                             + (6.0)*phi(xInd  +1, yInd, zIndM1)
                             -       phi(xInd  +2, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xP1yz   = pow(
                            (
                             - (2.0)*phi(xIndP1-1, yInd, zInd  )
                             - (3.0)*phi(xIndP1  , yInd, zInd  )
                             + (6.0)*phi(xIndP1+1, yInd, zInd  )
                             -       phi(xIndP1+2, yInd, zInd  )
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xP1yzP1 = pow(
                            (
                             - (2.0)*phi(xIndP1-1, yInd, zIndP1)
                             - (3.0)*phi(xIndP1  , yInd, zIndP1)
                             + (6.0)*phi(xIndP1+1, yInd, zIndP1)
                             -       phi(xIndP1+2, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xP1yzM1 = pow(
                            (
                             - (2.0)*phi(xIndP1-1, yInd, zIndM1)
                             - (3.0)*phi(xIndP1  , yInd, zIndM1)
                             + (6.0)*phi(xIndP1+1, yInd, zIndM1)
                             -       phi(xIndP1+2, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            // Different one sided stencil as M1 is already a ghost point
            phi_xM1yz   = pow(
                            (
                             - (11.0)*phi(xIndM1  , yInd, zInd  )
                             + (18.0)*phi(xIndM1+1, yInd, zInd  )
                             - ( 9.0)*phi(xIndM1+2, yInd, zInd  )
                             + ( 2.0)*phi(xIndM1+3, yInd, zInd  )
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xM1yzP1 = pow(
                            (
                             - (11.0)*phi(xIndM1  , yInd, zIndP1)
                             + (18.0)*phi(xIndM1+1, yInd, zIndP1)
                             - ( 9.0)*phi(xIndM1+2, yInd, zIndP1)
                             + ( 2.0)*phi(xIndM1+3, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xM1yzM1 = pow(
                            (
                             - (11.0)*phi(xIndM1  , yInd, zIndM1)
                             + (18.0)*phi(xIndM1+1, yInd, zIndM1)
                             - ( 9.0)*phi(xIndM1+2, yInd, zIndM1)
                             + ( 2.0)*phi(xIndM1+3, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);

            // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
            Jpp = 0.25*( (phi_xyzP1               - phi_xyzM1               )*
                         (n(xInd+1, yInd, zInd  ) - n(xInd-1, yInd, zInd  ) )
                         -
                         (phi_xP1yz               - phi_xM1yz                )*
                         (n(xInd  , yInd, zIndP1) - n(xInd  , yInd, zIndM1)) )
              / (mesh->dx(xInd,yInd) * mesh->dz);

            // J+x
            Jpx = 0.25*( n(xInd+1,yInd, zInd  )*(phi_xP1yzP1 - phi_xP1yzM1) -
                         n(xInd-1,yInd, zInd  )*(phi_xM1yzP1 - phi_xM1yzM1) -
                         n(xInd  ,yInd, zIndP1)*(phi_xP1yzP1 - phi_xM1yzP1) +
                         n(xInd  ,yInd, zIndM1)*(phi_xP1yzM1 - phi_xM1yzM1))
              / (mesh->dx(xInd,yInd) * mesh->dz);

            // Jx+
            Jxp = 0.25*( n(xInd+1, yInd, zIndP1)*(phi_xyzP1 - phi_xP1yz) -
                         n(xInd-1, yInd, zIndM1)*(phi_xM1yz - phi_xyzM1) -
                         n(xInd-1, yInd, zIndP1)*(phi_xyzP1 - phi_xM1yz) +
                         n(xInd+1, yInd, zIndM1)*(phi_xP1yz - phi_xyzM1))
              / (mesh->dx(xInd,yInd) * mesh->dz);


            result(xInd,yInd,zInd) = (Jpp + Jpx + Jxp) / 3.;
        }
    }

    // Treatment of the last point
    // That means that xP1 is a ghost/guard point
    xInd = mesh->xend;
    xIndP1 = xInd + 1;
    xIndM1 = xInd - 1;
    for(int yInd=mesh->ystart; yInd<=mesh->yend; yInd++){
        for(int zInd=0;zInd<ncz;zInd++) {
            zIndP1 = (zInd + 1) % ncz;
            zIndM1 = (zInd - 1 + ncz) % ncz;

            // Unsure why this needs to be one-sided, but is needed for convergence
            phi_xyzP1   = pow(
                            (
                                     phi(xInd  -2, yInd, zIndP1)
                             - (6.0)*phi(xInd  -1, yInd, zIndP1)
                             + (3.0)*phi(xInd    , yInd, zIndP1)
                             + (2.0)*phi(xInd  +1, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xyzM1   = pow(
                            (
                                     phi(xInd  -2, yInd, zIndM1)
                             - (6.0)*phi(xInd  -1, yInd, zIndM1)
                             + (3.0)*phi(xInd    , yInd, zIndM1)
                             + (2.0)*phi(xInd  +1, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xM1yz   = pow(
                            (
                                     phi(xIndM1-2, yInd, zInd  )
                             - (6.0)*phi(xIndM1-1, yInd, zInd  )
                             + (3.0)*phi(xIndM1  , yInd, zInd  )
                             + (2.0)*phi(xIndM1+1, yInd, zInd  )
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xM1yzP1 = pow(
                            (
                                     phi(xIndM1-2, yInd, zIndP1)
                             - (6.0)*phi(xIndM1-1, yInd, zIndP1)
                             + (3.0)*phi(xIndM1  , yInd, zIndP1)
                             + (2.0)*phi(xIndM1+1, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xM1yzM1 = pow(
                            (
                                     phi(xIndM1-2, yInd, zIndM1)
                             - (6.0)*phi(xIndM1-1, yInd, zIndM1)
                             + (3.0)*phi(xIndM1  , yInd, zIndM1)
                             + (2.0)*phi(xIndM1+1, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);

            // Different one sided stencil as P1 is already a ghost point
            phi_xP1yz   = pow(
                            (
                             - ( 2.0)*phi(xIndP1-3, yInd, zInd  )
                             + ( 9.0)*phi(xIndP1-2, yInd, zInd  )
                             - (18.0)*phi(xIndP1-1, yInd, zInd  )
                             + (11.0)*phi(xIndP1  , yInd, zInd  )
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xP1yzP1 = pow(
                            (
                             - ( 2.0)*phi(xIndP1-3, yInd, zIndP1)
                             + ( 9.0)*phi(xIndP1-2, yInd, zIndP1)
                             - (18.0)*phi(xIndP1-1, yInd, zIndP1)
                             + (11.0)*phi(xIndP1  , yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            phi_xP1yzM1 = pow(
                            (
                             - ( 2.0)*phi(xIndP1-3, yInd, zIndM1)
                             + ( 9.0)*phi(xIndP1-2, yInd, zIndM1)
                             - (18.0)*phi(xIndP1-1, yInd, zIndM1)
                             + (11.0)*phi(xIndP1  , yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);

            // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
            Jpp = 0.25*( (phi_xyzP1                 - phi_xyzM1               )*
                         (n(xInd+1, yInd, zInd  ) - n(xInd-1, yInd, zInd  ) )
                         -
                         (phi_xP1yz                 - phi_xM1yz                )*
                         (n(xInd  , yInd, zIndP1) - n(xInd  , yInd, zIndM1)) )
              / (mesh->dx(xInd,yInd) * mesh->dz);

            // J+x
            Jpx = 0.25*( n(xInd+1,yInd, zInd  )*(phi_xP1yzP1 - phi_xP1yzM1) -
                         n(xInd-1,yInd, zInd  )*(phi_xM1yzP1 - phi_xM1yzM1) -
                         n(xInd  ,yInd, zIndP1)*(phi_xP1yzP1 - phi_xM1yzP1) +
                         n(xInd  ,yInd, zIndM1)*(phi_xP1yzM1 - phi_xM1yzM1))
              / (mesh->dx(xInd,yInd) * mesh->dz);

            // Jx+
            Jxp = 0.25*( n(xInd+1, yInd, zIndP1)*(phi_xyzP1 - phi_xP1yz) -
                         n(xInd-1, yInd, zIndM1)*(phi_xM1yz - phi_xyzM1) -
                         n(xInd-1, yInd, zIndP1)*(phi_xyzP1 - phi_xM1yz) +
                         n(xInd+1, yInd, zIndM1)*(phi_xP1yz - phi_xyzM1))
              / (mesh->dx(xInd,yInd) * mesh->dz);


            result(xInd,yInd,zInd) = (Jpp + Jpx + Jxp) / 3.;
        }
    }

    return result;
}

// OwnOp3BasicBrackets

/*!
 * \brief Constructor
 *
 * Constructor which calls parent constructor and sets the bracket method
 */
OwnOp3BasicBrackets::OwnOp3BasicBrackets(Options *options) :
    OwnOperators(options)
{
    TRACE("Halt in OwnOp3BasicBrackets::OwnOp3BasicBrackets");

    bm = BRACKET_ARAKAWA;
}

/*!
 * Not implemented in this child class
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns phi (never reached)
 */
Field3D OwnOp3BasicBrackets::div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                                       const Field3D &phi)
{
    TRACE("Halt in OwnOp3BasicBrackets::div_uE_dot_grad_n_GradPerp_phi");

    throw BoutException("div_uE_dot_grad_n_GradPerp_phi not used in the "
                        "OwnOp3BasicBrackets implementation");

    return phi;
}

/*!
 * Calculates \f$\{\phi, \Omega^D\}\f$
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns result The result of the operation
 */
Field3D OwnOp3BasicBrackets::vortDAdv(const Field3D &phi, const Field3D &vortD)
{
    TRACE("Halt in OwnOp3BasicBrackets::vortDAdv");

    return invJ*bracket(phi, vortD, bm);
}

/*!
 * Calculates \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns result The result of the operation
 */
Field3D OwnOp3BasicBrackets::kinEnAdvN(const Field3D &phi, const Field3D &n)
{
    TRACE("Halt in OwnOp3BasicBrackets::kinEnAdvN");

    Field3D result;

    // Calculate the derivative of phi
    DDXPhi = DDX(phi);
    // Reset inner boundary
    ownBC.innerRhoCylinder(DDXPhi);
    // Reset outer boundary
    if (mesh->lastX()){
        /* NOTE: xend
         *       xend = index value of last inner point on this processor
         *       xend+1 = first guard point
         */
        ghostIndX = mesh->xend + 1;
        // Newton polynomial of fourth order (including boundary) evaluated at ghost
        for(int yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
            for(int zInd = 0; zInd < mesh->ngz -1; zInd ++){
                DDXPhi(ghostIndX, yInd, zInd) =
                    -      DDXPhi(ghostIndX-4, yInd, zInd)
                    +  4.0*DDXPhi(ghostIndX-3, yInd, zInd)
                    -  6.0*DDXPhi(ghostIndX-2, yInd, zInd)
                    +  4.0*DDXPhi(ghostIndX-1, yInd, zInd)
                      ;
            }
        }
    }

    // Communicate before taking new derivative
    mesh->communicate(DDXPhi);
    DDXDDXPhi = DDX(DDXPhi);
    // Reset inner boundary
    ownBC.innerRhoCylinder(DDXDDXPhi);

    if (mesh->lastX()){
        ghostIndX = mesh->xend + 1;
        // Newton polynomial of fourth order (including boundary) evaluated at ghost
        for(int yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
            for(int zInd = 0; zInd < mesh->ngz -1; zInd ++){
                DDXDDXPhi(ghostIndX, yInd, zInd) =
                    -      DDXDDXPhi(ghostIndX-4, yInd, zInd)
                    +  4.0*DDXDDXPhi(ghostIndX-3, yInd, zInd)
                    -  6.0*DDXDDXPhi(ghostIndX-2, yInd, zInd)
                    +  4.0*DDXDDXPhi(ghostIndX-1, yInd, zInd)
                    ;
            }
        }
    }

    // Communicate before taking new derivative
    mesh->communicate(DDXDDXPhi);

    result =   bracket(((DDXDDXPhi)^(2.0)), n, bm)
             + bracket(((invJ*DDZ(phi, true))^(2.0)), n, bm)
           ;

    // Multiply with B/2
    return 0.5*invJ*result;
}

#endif
