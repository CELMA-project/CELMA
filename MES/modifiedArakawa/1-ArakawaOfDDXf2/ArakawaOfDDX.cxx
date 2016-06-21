// *************** Simulation of ArakawaOfDDX *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "ArakawaOfDDX.hxx"

Field3D ArakawaOfDDXf2(Field3D const &f, Field3D const &g);

inline BoutReal arakawaDifferencing(
                                    BoutReal const & f_xyzP1  ,
                                    BoutReal const & f_xyzM1  ,
                                    BoutReal const & f_xP1yz  ,
                                    BoutReal const & f_xM1yz  ,
                                    BoutReal const & f_xP1yzP1,
                                    BoutReal const & f_xP1yzM1,
                                    BoutReal const & f_xM1yzP1,
                                    BoutReal const & f_xM1yzM1,
                                    Field3D  const & g        ,
                                    int      const & xInd     ,
                                    int      const & yInd     ,
                                    int      const & zInd     ,
                                    int      const & zIndP1   ,
                                    int      const & zIndM1
                                   )
                                   ;

// Initialization of the physics
// ############################################################################
int ArakawaOfDDX::init(bool restarting) {
    TRACE("Halt in ArakawaOfDDXf2::init");

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Load from the geometry
    // ************************************************************************
    Options *geom = options->getSection("geom");
    geom->get("Lx", Lx, 0.0);
    // ************************************************************************

    // Obtain the fields
    // ************************************************************************
    // phi
    phi = FieldFactory::get()
        ->create3D("phi:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // n
    n = FieldFactory::get()
        ->create3D("n:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // S
    S = FieldFactory::get()
        ->create3D("S:solution", Options::getRoot(), mesh, CELL_CENTRE, 0);
    // ************************************************************************

    // Add a FieldGroup to communicate
    // ************************************************************************
    // Only these fields will be taken derivatives of
    com_group.add(phi);
    com_group.add(n);
    // ************************************************************************

    // Set boundaries manually
    // ************************************************************************
    phi.setBoundary("phi");
    phi.applyBoundary();
    ownBC.innerRhoCylinder(phi);
    n.setBoundary("n");
    n.applyBoundary();
    ownBC.innerRhoCylinder(n);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(com_group);

    output << "\n\n\n\n\n\n\nNow running test" << std::endl;

    // Calculate
    S_num = ArakawaOfDDXf2(phi, n);

    // Error in S
    e = S_num - S;

    // Save the variables
    SAVE_ONCE (Lx);
    SAVE_ONCE4(phi, n, S, S_num);
    SAVE_ONCE (e);

    // Finalize
    dump.write();
    dump.close();

    output << "\nFinished running test, now quitting\n\n\n\n\n\n" << std::endl;

    // Wait for all processors to write data
    MPI_Barrier(BoutComm::get());

    return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int ArakawaOfDDX::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(ArakawaOfDDX);

// Destructor
ArakawaOfDDX::~ArakawaOfDDX(){
    TRACE("ArakawaOfDDXf2::~ArakawaOfDDXf2");
}

// New implementation
// NOTE: When calculating the last point, we need either 2 ghost points, or a
// one-sided stencil.
Field3D ArakawaOfDDXf2(Field3D const &f, Field3D const &g)
{
    TRACE("ArakawaOfDDXf2");

    Field3D result;
    int xInd;
    int xIndP1;
    int xIndM1;
    int xstart;
    int xend  ;
    int zIndP1 ;
    int zIndM1 ;
    BoutReal f_xyzP1;
    BoutReal f_xyzM1;
    BoutReal f_xP1yz;
    BoutReal f_xM1yz;
    BoutReal f_xP1yzP1;
    BoutReal f_xP1yzM1;
    BoutReal f_xM1yzP1;
    BoutReal f_xM1yzM1;

    result.allocate();

    int ncz = mesh->ngz - 1;
    /* Loop over all inner points:
     * The Arakawa bracket requires evaluation of the field in x-1 and x+1. In
     * our case, we are at the same time taking the derivative of the fields,
     * which means that we need to calculate the derivatives in x-1 and x+1
     * (i.e. the ghost points/guard cells). Hence, special treatment is needed
     * for two points closest to the domain/processor boundary
     */
    xstart = mesh->xstart + 2;
    xend   = mesh->xend   - 2;
    for(xInd=xstart; xInd<=xend; xInd++){
        xIndP1 = xInd + 1;
        xIndM1 = xInd - 1;
        for(int yInd=mesh->ystart; yInd<=mesh->yend; yInd++){
            for(int zInd=0;zInd<ncz;zInd++) {
                int zIndP1 = (zInd + 1) % ncz;
                int zIndM1 = (zInd - 1 + ncz) % ncz;

                /* Boilerplate
                 * f_xyz = pow( (-f(xInd  -1, yInd, zInd  ) + f(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                 */
                f_xyzP1   = pow( (-f(xInd  -1, yInd, zIndP1) + f(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                f_xyzM1   = pow( (-f(xInd  -1, yInd, zIndM1) + f(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                f_xP1yz   = pow( (-f(xIndP1-1, yInd, zInd  ) + f(xIndP1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                f_xM1yz   = pow( (-f(xIndM1-1, yInd, zInd  ) + f(xIndM1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                f_xP1yzP1 = pow( (-f(xIndP1-1, yInd, zIndP1) + f(xIndP1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                f_xP1yzM1 = pow( (-f(xIndP1-1, yInd, zIndM1) + f(xIndP1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                f_xM1yzP1 = pow( (-f(xIndM1-1, yInd, zIndP1) + f(xIndM1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                f_xM1yzM1 = pow( (-f(xIndM1-1, yInd, zIndM1) + f(xIndM1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);


                result(xInd,yInd,zInd) = arakawaDifferencing(f_xyzP1  ,
                                                             f_xyzM1  ,
                                                             f_xP1yz  ,
                                                             f_xM1yz  ,
                                                             f_xP1yzP1,
                                                             f_xP1yzM1,
                                                             f_xM1yzP1,
                                                             f_xM1yzM1,
                                                             g        ,
                                                             xInd     ,
                                                             yInd     ,
                                                             zInd     ,
                                                             zIndP1   ,
                                                             zIndM1
                                                            );
            }
        }
    }

    // Treatment of the second first point (no special treatment of the innermost rho, as we will only use one ghost point)
    // That means that xM1 - 1 is a ghost/guard point
    xInd = mesh->xstart + 1;
    xIndP1 = xInd + 1;
    xIndM1 = xInd - 1;
    for(int yInd=mesh->ystart; yInd<=mesh->yend; yInd++){
        for(int zInd=0; zInd<ncz; zInd++) {
            zIndP1 = (zInd + 1) % ncz;
            zIndM1 = (zInd - 1 + ncz) % ncz;

            /* Boilerplate
             * BoutReal f_xyz = pow( (-f(xInd  -1, yInd, zInd  ) + f(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
             */
            f_xyzP1   = pow( (-f(xInd  -1, yInd, zIndP1) + f(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            f_xyzM1   = pow( (-f(xInd  -1, yInd, zIndM1) + f(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            f_xP1yz   = pow( (-f(xIndP1-1, yInd, zInd  ) + f(xIndP1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            f_xP1yzP1 = pow( (-f(xIndP1-1, yInd, zIndP1) + f(xIndP1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            f_xP1yzM1 = pow( (-f(xIndP1-1, yInd, zIndM1) + f(xIndP1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            // One-sides stencils for these guys
            f_xM1yz   = pow(
                             (
                              - (2.0)*f(xIndM1-1, yInd, zInd  )
                              - (3.0)*f(xIndM1  , yInd, zInd  )
                              + (6.0)*f(xIndM1+1, yInd, zInd  )
                              -       f(xIndM1+2, yInd, zInd  )
                             )/(6.0*mesh->dx(xInd,yInd)),
                           2.0);
            f_xM1yzP1 = pow(
                             (
                              - (2.0)*f(xIndM1-1, yInd, zIndP1)
                              - (3.0)*f(xIndM1  , yInd, zIndP1)
                              + (6.0)*f(xIndM1+1, yInd, zIndP1)
                              -       f(xIndM1+2, yInd, zIndP1)
                             )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xM1yzM1 = pow(
                             (
                              - (2.0)*f(xIndM1-1, yInd, zIndM1)
                              - (3.0)*f(xIndM1  , yInd, zIndM1)
                              + (6.0)*f(xIndM1+1, yInd, zIndP1)
                              -       f(xIndM1+2, yInd, zIndP1)
                             )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);

            result(xInd,yInd,zInd) = arakawaDifferencing(f_xyzP1  ,
                                                         f_xyzM1  ,
                                                         f_xP1yz  ,
                                                         f_xM1yz  ,
                                                         f_xP1yzP1,
                                                         f_xP1yzM1,
                                                         f_xM1yzP1,
                                                         f_xM1yzM1,
                                                         g        ,
                                                         xInd     ,
                                                         yInd     ,
                                                         zInd     ,
                                                         zIndP1   ,
                                                         zIndM1
                                                        );
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

            f_xyzP1   = pow(
                            (
                             - (2.0)*f(xInd  -1, yInd, zIndP1)
                             - (3.0)*f(xInd    , yInd, zIndP1)
                             + (6.0)*f(xInd  +1, yInd, zIndP1)
                             -       f(xInd  +2, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xyzM1   = pow(
                            (
                             - (2.0)*f(xInd  -1, yInd, zIndM1)
                             - (3.0)*f(xInd    , yInd, zIndM1)
                             + (6.0)*f(xInd  +1, yInd, zIndM1)
                             -       f(xInd  +2, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xP1yz   = pow(
                            (
                             - (2.0)*f(xIndP1-1, yInd, zInd  )
                             - (3.0)*f(xIndP1  , yInd, zInd  )
                             + (6.0)*f(xIndP1+1, yInd, zInd  )
                             -       f(xIndP1+2, yInd, zInd  )
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xP1yzP1 = pow(
                            (
                             - (2.0)*f(xIndP1-1, yInd, zIndP1)
                             - (3.0)*f(xIndP1  , yInd, zIndP1)
                             + (6.0)*f(xIndP1+1, yInd, zIndP1)
                             -       f(xIndP1+2, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xP1yzM1 = pow(
                            (
                             - (2.0)*f(xIndP1-1, yInd, zIndM1)
                             - (3.0)*f(xIndP1  , yInd, zIndM1)
                             + (6.0)*f(xIndP1+1, yInd, zIndM1)
                             -       f(xIndP1+2, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            // Different one sided stencil as M1 is already a ghost point
            f_xM1yz   = pow(
                            (
                             - (11.0)*f(xIndM1  , yInd, zInd  )
                             + (18.0)*f(xIndM1+1, yInd, zInd  )
                             - ( 9.0)*f(xIndM1+2, yInd, zInd  )
                             + ( 2.0)*f(xIndM1+3, yInd, zInd  )
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xM1yzP1 = pow(
                            (
                             - (11.0)*f(xIndM1  , yInd, zIndP1)
                             + (18.0)*f(xIndM1+1, yInd, zIndP1)
                             - ( 9.0)*f(xIndM1+2, yInd, zIndP1)
                             + ( 2.0)*f(xIndM1+3, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xM1yzM1 = pow(
                            (
                             - (11.0)*f(xIndM1  , yInd, zIndM1)
                             + (18.0)*f(xIndM1+1, yInd, zIndM1)
                             - ( 9.0)*f(xIndM1+2, yInd, zIndM1)
                             + ( 2.0)*f(xIndM1+3, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);

            result(xInd,yInd,zInd) = arakawaDifferencing(f_xyzP1  ,
                                                         f_xyzM1  ,
                                                         f_xP1yz  ,
                                                         f_xM1yz  ,
                                                         f_xP1yzP1,
                                                         f_xP1yzM1,
                                                         f_xM1yzP1,
                                                         f_xM1yzM1,
                                                         g        ,
                                                         xInd     ,
                                                         yInd     ,
                                                         zInd     ,
                                                         zIndP1   ,
                                                         zIndM1
                                                        );
        }
    }

    // Treatment of the second last point (no special treatment of the innermost rho, as we will only use one ghost point)
    // That means that xP1 + 1 is a ghost/guard point
    xInd = mesh->xend - 1;
    xIndP1 = xInd + 1;
    xIndM1 = xInd - 1;
    for(int yInd=mesh->ystart; yInd<=mesh->yend; yInd++){
        for(int zInd=0; zInd<ncz; zInd++) {
            zIndP1 = (zInd + 1) % ncz;
            zIndM1 = (zInd - 1 + ncz) % ncz;

            /* Boilerplate
             * BoutReal f_xyz = pow( (-f(xInd  -1, yInd, zInd  ) + f(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
             */
            f_xyzP1   = pow( (-f(xInd  -1, yInd, zIndP1) + f(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            f_xyzM1   = pow( (-f(xInd  -1, yInd, zIndM1) + f(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            f_xM1yz   = pow( (-f(xIndM1-1, yInd, zInd  ) + f(xIndM1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            f_xM1yzP1 = pow( (-f(xIndM1-1, yInd, zIndP1) + f(xIndM1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            f_xM1yzM1 = pow( (-f(xIndM1-1, yInd, zIndM1) + f(xIndM1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
            // One-sides stencils for these guys
            f_xP1yz   = pow(
                            (
                                     f(xIndP1-2, yInd, zInd  )
                             - (6.0)*f(xIndP1-1, yInd, zInd  )
                             + (3.0)*f(xIndP1  , yInd, zInd  )
                             + (2.0)*f(xIndP1+1, yInd, zInd  )
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xP1yzP1 = pow(
                            (
                                      f(xIndP1-2, yInd, zIndP1)
                              - (6.0)*f(xIndP1-1, yInd, zIndP1)
                              + (3.0)*f(xIndP1  , yInd, zIndP1)
                              + (2.0)*f(xIndP1+1, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xP1yzM1 = pow(
                            (
                                     f(xIndP1-2, yInd, zIndM1)
                             - (6.0)*f(xIndP1-1, yInd, zIndM1)
                             + (3.0)*f(xIndP1  , yInd, zIndM1)
                             + (2.0)*f(xIndP1+1, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);

            result(xInd,yInd,zInd) = arakawaDifferencing(f_xyzP1  ,
                                                         f_xyzM1  ,
                                                         f_xP1yz  ,
                                                         f_xM1yz  ,
                                                         f_xP1yzP1,
                                                         f_xP1yzM1,
                                                         f_xM1yzP1,
                                                         f_xM1yzM1,
                                                         g        ,
                                                         xInd     ,
                                                         yInd     ,
                                                         zInd     ,
                                                         zIndP1   ,
                                                         zIndM1
                                                        );
        }
    }

    // Treatment of the first point
    // That means that xP1 is a ghost/guard point
    xInd = mesh->xend;
    xIndP1 = xInd + 1;
    xIndM1 = xInd - 1;
    for(int yInd=mesh->ystart; yInd<=mesh->yend; yInd++){
        for(int zInd=0;zInd<ncz;zInd++) {
            zIndP1 = (zInd + 1) % ncz;
            zIndM1 = (zInd - 1 + ncz) % ncz;

            f_xyzP1   = pow(
                            (
                                     f(xInd  -2, yInd, zIndP1)
                             - (6.0)*f(xInd  -1, yInd, zIndP1)
                             + (3.0)*f(xInd    , yInd, zIndP1)
                             + (2.0)*f(xInd  +1, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xyzM1   = pow(
                            (
                                     f(xInd  -2, yInd, zIndM1)
                             - (6.0)*f(xInd  -1, yInd, zIndM1)
                             + (3.0)*f(xInd    , yInd, zIndM1)
                             + (2.0)*f(xInd  +1, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xM1yz   = pow(
                            (
                                     f(xIndM1-2, yInd, zInd  )
                             - (6.0)*f(xIndM1-1, yInd, zInd  )
                             + (3.0)*f(xIndM1  , yInd, zInd  )
                             + (2.0)*f(xIndM1+1, yInd, zInd  )
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xM1yzP1 = pow(
                            (
                                     f(xIndM1-2, yInd, zIndP1)
                             - (6.0)*f(xIndM1-1, yInd, zIndP1)
                             + (3.0)*f(xIndM1  , yInd, zIndP1)
                             + (2.0)*f(xIndM1+1, yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xM1yzM1 = pow(
                            (
                                     f(xIndM1-2, yInd, zIndM1)
                             - (6.0)*f(xIndM1-1, yInd, zIndM1)
                             + (3.0)*f(xIndM1  , yInd, zIndM1)
                             + (2.0)*f(xIndM1+1, yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);

            // Different one sided stencil as P1 is already a ghost point
            f_xP1yz   = pow(
                            (
                             - ( 2.0)*f(xIndP1-3, yInd, zInd  )
                             + ( 9.0)*f(xIndP1-2, yInd, zInd  )
                             - (18.0)*f(xIndP1-1, yInd, zInd  )
                             + (11.0)*f(xIndP1  , yInd, zInd  )
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xP1yzP1 = pow(
                            (
                             - ( 2.0)*f(xIndP1-3, yInd, zIndP1)
                             + ( 9.0)*f(xIndP1-2, yInd, zIndP1)
                             - (18.0)*f(xIndP1-1, yInd, zIndP1)
                             + (11.0)*f(xIndP1  , yInd, zIndP1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);
            f_xP1yzM1 = pow(
                            (
                             - ( 2.0)*f(xIndP1-3, yInd, zIndM1)
                             + ( 9.0)*f(xIndP1-2, yInd, zIndM1)
                             - (18.0)*f(xIndP1-1, yInd, zIndM1)
                             + (11.0)*f(xIndP1  , yInd, zIndM1)
                            )/(6.0*mesh->dx(xInd,yInd)),
                            2.0);

            result(xInd,yInd,zInd) = arakawaDifferencing(f_xyzP1  ,
                                                         f_xyzM1  ,
                                                         f_xP1yz  ,
                                                         f_xM1yz  ,
                                                         f_xP1yzP1,
                                                         f_xP1yzM1,
                                                         f_xM1yzP1,
                                                         f_xM1yzM1,
                                                         g        ,
                                                         xInd     ,
                                                         yInd     ,
                                                         zInd     ,
                                                         zIndP1   ,
                                                         zIndM1
                                                        );
        }
    }

    return result;
}

// The arakawa diffrencing
BoutReal arakawaDifferencing(
                              BoutReal const & f_xyzP1  ,
                              BoutReal const & f_xyzM1  ,
                              BoutReal const & f_xP1yz  ,
                              BoutReal const & f_xM1yz  ,
                              BoutReal const & f_xP1yzP1,
                              BoutReal const & f_xP1yzM1,
                              BoutReal const & f_xM1yzP1,
                              BoutReal const & f_xM1yzM1,
                              Field3D  const & g        ,
                              int      const & xInd     ,
                              int      const & yInd     ,
                              int      const & zInd     ,
                              int      const & zIndP1   ,
                              int      const & zIndM1
                             )
{
            BoutReal Jpp;
            BoutReal Jpx;
            BoutReal Jxp;

            // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
            Jpp = 0.25*( (f_xyzP1                 - f_xyzM1               )*
                         (g(xInd+1, yInd, zInd  ) - g(xInd-1, yInd, zInd  ) )
                         -
                         (f_xP1yz                 - f_xM1yz                )*
                         (g(xInd  , yInd, zIndP1) - g(xInd  , yInd, zIndM1)) )
              / (mesh->dx(xInd,yInd) * mesh->dz);

            // J+x
            Jpx = 0.25*( g(xInd+1,yInd, zInd  )*(f_xP1yzP1 - f_xP1yzM1) -
                         g(xInd-1,yInd, zInd  )*(f_xM1yzP1 - f_xM1yzM1) -
                         g(xInd  ,yInd, zIndP1)*(f_xP1yzP1 - f_xM1yzP1) +
                         g(xInd  ,yInd, zIndM1)*(f_xP1yzM1 - f_xM1yzM1))
              / (mesh->dx(xInd,yInd) * mesh->dz);

            // Jx+
            Jxp = 0.25*( g(xInd+1, yInd, zIndP1)*(f_xyzP1 - f_xP1yz) -
                         g(xInd-1, yInd, zIndM1)*(f_xM1yz - f_xyzM1) -
                         g(xInd-1, yInd, zIndP1)*(f_xyzP1 - f_xM1yz) +
                         g(xInd+1, yInd, zIndM1)*(f_xP1yz - f_xyzM1))
              / (mesh->dx(xInd,yInd) * mesh->dz);

            return (Jpp + Jpx + Jxp) / 3.;
}
