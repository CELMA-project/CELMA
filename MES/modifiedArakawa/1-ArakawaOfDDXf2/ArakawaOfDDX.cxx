// *************** Simulation of ArakawaOfDDX *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "ArakawaOfDDX.hxx"

Field3D ArakawaOfDDXf2(Field3D const &f, Field3D const &g);

// Initialization of the physics
// ############################################################################
int ArakawaOfDDXf2::init(bool restarting) {
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
    // f
    f = FieldFactory::get()
        ->create3D("f:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // S
    S = FieldFactory::get()
        ->create3D("S:solution", Options::getRoot(), mesh, CELL_CENTRE, 0);
    // ************************************************************************

    // Add a FieldGroup to communicate
    // ************************************************************************
    // Only these fields will be taken derivatives of
    com_group.add(f);
    // ************************************************************************

    // Set boundaries manually
    // ************************************************************************
    f.setBoundary("f");
    f.applyBoundary();
    ownBC.innerRhoCylinder(f);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(com_group);

    output << "\n\n\n\n\n\n\nNow running test" << std::endl;

    // Calculate
    S_num = DDXCylHighResCenter(f);

    // Error in S
    e = S_num - S;

    // Save the variables
    SAVE_ONCE (Lx);
    SAVE_ONCE3(f, S, S_num);
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
int ArakawaOfDDXf2::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(ArakawaOfDDXf2);

// Destructor
ArakawaOfDDXf2::~ArakawaOfDDXf2(){
    TRACE("ArakawaOfDDXf2::~ArakawaOfDDXf2");
}

// New implementation
// NOTE: When calculating the last point, we need either 2 ghost points, or a
// one-sided stencil. We will use a one-sided stencil, and we will not have to
// worry about the precision of the ghost point if we are using a neumann BC,
// as this is known to machine precision when boundaries are half between grid
// points
Field3D ArakawaOfDDXf2(Field3D const &f, Field3D const &g)
{
    TRACE("ArakawaOfDDXf2");

    Field3D result;

    result.allocate();

    int ncz = mesh->ngz - 1;
    if(!mesh->firstX() && !mesh->lastX()){
        for(int xInd=mesh->xstart;xInd<=mesh->xend;xInd++){
            for(int yInd=mesh->ystart;yInd<=mesh->yend;yInd++){
                for(int zInd=0;zInd<ncz;zInd++) {
                    int zIndP1 = (zInd + 1) % ncz;
                    int zIndM1 = (zInd - 1 + ncz) % ncz;

                    /* Boilerplate
                     * BoutReal f_xyz = pow( (-f(xInd  -1, yInd, zInd  ) + f(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                     */
                    int xIndP1 = xInd + 1;
                    int xIndM1 = xInd - 1;
                    BoutReal f_xyzP1   = pow( (-f(xInd  -1, yInd, zIndP1) + f(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xyzM1   = pow( (-f(xInd  -1, yInd, zIndM1) + f(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xP1yz   = pow( (-f(xIndP1-1, yInd, zInd  ) + f(xIndP1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xM1yz   = pow( (-f(xIndM1-1, yInd, zInd  ) + f(xIndM1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xP1yzP1 = pow( (-f(xIndP1-1, yInd, zIndP1) + f(xIndP1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xP1yzM1 = pow( (-f(xIndP1-1, yInd, zIndM1) + f(xIndP1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xM1yzP1 = pow( (-f(xIndM1-1, yInd, zIndP1) + f(xIndM1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xM1yzM1 = pow( (-f(xIndM1-1, yInd, zIndM1) + f(xIndM1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);

                    // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
                    BoutReal Jpp = 0.25*( (f_xyzP1                 - f_xyzM1               )*
                                          (g(xInd+1, yInd, zInd  ) - g(xInd-1, yInd, zInd  ) )
                                          -
                                          (f_xP1yz                 - f_xM1yz                )*
                                          (g(xInd  , yInd, zIndP1) - g(xInd  , yInd, zIndM1)) )
                      / (mesh->dx(xInd,yInd) * mesh->dz);

                    // J+x
                    BoutReal Jpx = 0.25*( g(xInd+1,yInd, zInd  )*(f_xP1yzP1 - f_xP1yzM1) -
                                          g(xInd-1,yInd, zInd  )*(f_xM1yzP1 - f_xM1yzM1) -
                                          g(xInd  ,yInd, zIndP1)*(f_xP1yzP1 - f_xM1yzP1) +
                                          g(xInd  ,yInd, zIndM1)*(f_xP1yzM1 - f_xM1yzM1))
                      / (mesh->dx(xInd,yInd) * mesh->dz);
                    // Jx+
                    BoutReal Jxp = 0.25*( g(xInd+1, yInd, zIndP1)*(f_xyzP1 - f_xP1yz) -
                                          g(xInd-1, yInd, zIndM1)*(f_xM1yz - f_xyzM1) -
                                          g(xInd-1, yInd, zIndP1)*(f_xyzP1 - f_xM1yz) +
                                          g(xInd+1, yInd, zIndM1)*(f_xP1yz - f_xyzM1))
                      / (mesh->dx(xInd,yInd) * mesh->dz);

                    result(xInd,yInd,zInd) = (Jpp + Jpx + Jxp) / 3.;
                }
            }
        }
    }
    else if(mesh->firstX()){
        // Care must be taken at the inner boundaries (2 ghost points are needed), therefore, we start one from the innermost xInd
        for(int xInd=mesh->xstart+1;xInd<=mesh->xend;xInd++){
            for(int yInd=mesh->ystart;yInd<=mesh->yend;yInd++){
                for(int zInd=0;zInd<ncz;zInd++) {
                    int zIndP1 = (zInd + 1) % ncz;
                    int zIndM1 = (zInd - 1 + ncz) % ncz;

                    /* Boilerplate
                     * BoutReal f_xyz = pow( (-f(xInd  -1, yInd, zInd  ) + f(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                     */
                    int xIndP1 = xInd + 1;
                    int xIndM1 = xInd - 1;
                    BoutReal f_xyzP1   = pow( (-f(xInd  -1, yInd, zIndP1) + f(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xyzM1   = pow( (-f(xInd  -1, yInd, zIndM1) + f(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xP1yz   = pow( (-f(xIndP1-1, yInd, zInd  ) + f(xIndP1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xM1yz   = pow( (-f(xIndM1-1, yInd, zInd  ) + f(xIndM1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xP1yzP1 = pow( (-f(xIndP1-1, yInd, zIndP1) + f(xIndP1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xP1yzM1 = pow( (-f(xIndP1-1, yInd, zIndM1) + f(xIndP1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xM1yzP1 = pow( (-f(xIndM1-1, yInd, zIndP1) + f(xIndM1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xM1yzM1 = pow( (-f(xIndM1-1, yInd, zIndM1) + f(xIndM1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);

                    // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
                    BoutReal Jpp = 0.25*( (f_xyzP1                 - f_xyzM1               )*
                                          (g(xInd+1, yInd, zInd  ) - g(xInd-1, yInd, zInd  ) )
                                          -
                                          (f_xP1yz                 - f_xM1yz                )*
                                          (g(xInd  , yInd, zIndP1) - g(xInd  , yInd, zIndM1)) )
                      / (mesh->dx(xInd,yInd) * mesh->dz);

                    // J+x
                    BoutReal Jpx = 0.25*( g(xInd+1,yInd, zInd  )*(f_xP1yzP1 - f_xP1yzM1) -
                                          g(xInd-1,yInd, zInd  )*(f_xM1yzP1 - f_xM1yzM1) -
                                          g(xInd  ,yInd, zIndP1)*(f_xP1yzP1 - f_xM1yzP1) +
                                          g(xInd  ,yInd, zIndM1)*(f_xP1yzM1 - f_xM1yzM1))
                      / (mesh->dx(xInd,yInd) * mesh->dz);
                    // Jx+
                    BoutReal Jxp = 0.25*( g(xInd+1, yInd, zIndP1)*(f_xyzP1 - f_xP1yz) -
                                          g(xInd-1, yInd, zIndM1)*(f_xM1yz - f_xyzM1) -
                                          g(xInd-1, yInd, zIndP1)*(f_xyzP1 - f_xM1yz) +
                                          g(xInd+1, yInd, zIndM1)*(f_xP1yz - f_xyzM1))
                      / (mesh->dx(xInd,yInd) * mesh->dz);

                    result(xInd,yInd,zInd) = (Jpp + Jpx + Jxp) / 3.;
                }
            }
        }

        int piIndex = ncz/2;
        // Take derivative of first point only
        // NOTE: First ghost-point is without any problem
        // Loop over first half of cylinder
        int xInd=mesh->xstart;
        for(int yInd=mesh->ystart;yInd<=mesh->yend;yInd++){
            // Loop until pi index
            for(int zInd=0; zInd<piIndex; zInd++) {
                int zIndP1 = (zInd + 1) % ncz;
                int zIndM1 = (zInd - 1 + ncz) % ncz;

                /* Boilerplate
                 * BoutReal f_xyz = pow( (-f(xInd  -1, yInd, zInd  ) + f(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                 */
                int xIndP1 = xInd + 1;
                int xIndM1 = xInd - 1;
                BoutReal f_xyzP1   = pow( (-f(xInd  -1, yInd, zIndP1) + f(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xyzM1   = pow( (-f(xInd  -1, yInd, zIndM1) + f(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xP1yz   = pow( (-f(xIndP1-1, yInd, zInd  ) + f(xIndP1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xP1yzP1 = pow( (-f(xIndP1-1, yInd, zIndP1) + f(xIndP1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xP1yzM1 = pow( (-f(xIndP1-1, yInd, zIndM1) + f(xIndP1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                // Two ghost points are needed. Will use same tricks as in innerRhoBoundary
                BoutReal f_xM1yz   = pow( (-f(xInd+2, yInd, zInd   + piIndex) + f(xIndM1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                // Note: If zInd = piIndex will remain witin scope due to the addition of the piIndex
                BoutReal f_xM1yzM1 = pow( (-f(xInd+2, yInd, zInd-1 + piIndex) + f(xIndM1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                // Note: If zInd = piIndex will overflow due to the piIndex (unless we wrap around)
                int zIndP1PPi = (zInd + 1 + piIndex) % ncz;
                BoutReal f_xM1yzP1 = pow( (-f(xInd+2, yInd, zIndP1PPi) + f(xIndM1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);

                // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
                BoutReal Jpp = 0.25*( (f_xyzP1                 - f_xyzM1               )*
                                      (g(xInd+1, yInd, zInd  ) - g(xInd-1, yInd, zInd  ) )
                                      -
                                      (f_xP1yz                 - f_xM1yz                )*
                                      (g(xInd  , yInd, zIndP1) - g(xInd  , yInd, zIndM1)) )
                  / (mesh->dx(xInd,yInd) * mesh->dz);

                // J+x
                BoutReal Jpx = 0.25*( g(xInd+1,yInd, zInd  )*(f_xP1yzP1 - f_xP1yzM1) -
                                      g(xInd-1,yInd, zInd  )*(f_xM1yzP1 - f_xM1yzM1) -
                                      g(xInd  ,yInd, zIndP1)*(f_xP1yzP1 - f_xM1yzP1) +
                                      g(xInd  ,yInd, zIndM1)*(f_xP1yzM1 - f_xM1yzM1))
                  / (mesh->dx(xInd,yInd) * mesh->dz);
                // Jx+
                BoutReal Jxp = 0.25*( g(xInd+1, yInd, zIndP1)*(f_xyzP1 - f_xP1yz) -
                                      g(xInd-1, yInd, zIndM1)*(f_xM1yz - f_xyzM1) -
                                      g(xInd-1, yInd, zIndP1)*(f_xyzP1 - f_xM1yz) +
                                      g(xInd+1, yInd, zIndM1)*(f_xP1yz - f_xyzM1))
                  / (mesh->dx(xInd,yInd) * mesh->dz);

                result(xInd,yInd,zInd) = (Jpp + Jpx + Jxp) / 3.;
            }
        }

        // Loop over second half of cylinder
        for(int yInd=mesh->ystart;yInd<=mesh->yend;yInd++){
            // Loop until pi index
            for(int zInd=piIndex; zInd<ncz; zInd++) {
                int zIndP1 = (zInd + 1) % ncz;
                int zIndM1 = (zInd - 1 + ncz) % ncz;

                /* Boilerplate
                 * BoutReal f_xyz = pow( (-f(xInd  -1, yInd, zInd  ) + f(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                 */
                int xIndP1 = xInd + 1;
                int xIndM1 = xInd - 1;
                BoutReal f_xyzP1   = pow( (-f(xInd  -1, yInd, zIndP1) + f(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xyzM1   = pow( (-f(xInd  -1, yInd, zIndM1) + f(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xP1yz   = pow( (-f(xIndP1-1, yInd, zInd  ) + f(xIndP1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xP1yzP1 = pow( (-f(xIndP1-1, yInd, zIndP1) + f(xIndP1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xP1yzM1 = pow( (-f(xIndP1-1, yInd, zIndM1) + f(xIndP1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                // Two ghost points are needed. Will use same tricks as in innerRhoBoundary
                BoutReal f_xM1yz   = pow( (-f(xInd+2, yInd, zInd   - piIndex) + f(xIndM1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                // Note: If zInd = piIndex will remain witin scope due to the addition of the piIndex
                BoutReal f_xM1yzP1 = pow( (-f(xInd+2, yInd, zInd+1 - piIndex) + f(xIndM1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                // Note: If zInd = piIndex will underflow due to the piIndex (unless we wrap around)
                int zIndM1MPi = (zInd - 1 - piIndex + ncz) % ncz;
                BoutReal f_xM1yzM1 = pow( (-f(xInd+2, yInd, zIndM1MPi) + f(xIndM1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);

                // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
                BoutReal Jpp = 0.25*( (f_xyzP1                 - f_xyzM1               )*
                                      (g(xInd+1, yInd, zInd  ) - g(xInd-1, yInd, zInd  ) )
                                      -
                                      (f_xP1yz                 - f_xM1yz                )*
                                      (g(xInd  , yInd, zIndP1) - g(xInd  , yInd, zIndM1)) )
                  / (mesh->dx(xInd,yInd) * mesh->dz);

                // J+x
                BoutReal Jpx = 0.25*( g(xInd+1,yInd, zInd  )*(f_xP1yzP1 - f_xP1yzM1) -
                                      g(xInd-1,yInd, zInd  )*(f_xM1yzP1 - f_xM1yzM1) -
                                      g(xInd  ,yInd, zIndP1)*(f_xP1yzP1 - f_xM1yzP1) +
                                      g(xInd  ,yInd, zIndM1)*(f_xP1yzM1 - f_xM1yzM1))
                  / (mesh->dx(xInd,yInd) * mesh->dz);
                // Jx+
                BoutReal Jxp = 0.25*( g(xInd+1, yInd, zIndP1)*(f_xyzP1 - f_xP1yz) -
                                      g(xInd-1, yInd, zIndM1)*(f_xM1yz - f_xyzM1) -
                                      g(xInd-1, yInd, zIndP1)*(f_xyzP1 - f_xM1yz) +
                                      g(xInd+1, yInd, zIndM1)*(f_xP1yz - f_xyzM1))
                  / (mesh->dx(xInd,yInd) * mesh->dz);

                result(xInd,yInd,zInd) = (Jpp + Jpx + Jxp) / 3.;
            }
        }
    }
    else if(mesh->lastX()){
        // Loop until point before last inner point
        // NOTE: First ghost-point is without any problem
        for(int xInd=mesh->xstart;xInd<=mesh->xend-1;xInd++){
            for(int yInd=mesh->ystart;yInd<=mesh->yend;yInd++){
                for(int zInd=0;zInd<ncz;zInd++) {
                    int zIndP1 = (zInd + 1) % ncz;
                    int zIndM1 = (zInd - 1 + ncz) % ncz;

                    /* Boilerplate
                     * BoutReal f_xyz = pow( (-f(xInd  -1, yInd, zInd  ) + f(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                     */
                    int xIndP1 = xInd + 1;
                    int xIndM1 = xInd - 1;
                    BoutReal f_xyzP1   = pow( (-f(xInd  -1, yInd, zIndP1) + f(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xyzM1   = pow( (-f(xInd  -1, yInd, zIndM1) + f(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xP1yz   = pow( (-f(xIndP1-1, yInd, zInd  ) + f(xIndP1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xM1yz   = pow( (-f(xIndM1-1, yInd, zInd  ) + f(xIndM1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xP1yzP1 = pow( (-f(xIndP1-1, yInd, zIndP1) + f(xIndP1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xP1yzM1 = pow( (-f(xIndP1-1, yInd, zIndM1) + f(xIndP1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xM1yzP1 = pow( (-f(xIndM1-1, yInd, zIndP1) + f(xIndM1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                    BoutReal f_xM1yzM1 = pow( (-f(xIndM1-1, yInd, zIndM1) + f(xIndM1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);

                    // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
                    BoutReal Jpp = 0.25*( (f_xyzP1                 - f_xyzM1               )*
                                          (g(xInd+1, yInd, zInd  ) - g(xInd-1, yInd, zInd  ) )
                                          -
                                          (f_xP1yz                 - f_xM1yz                )*
                                          (g(xInd  , yInd, zIndP1) - g(xInd  , yInd, zIndM1)) )
                      / (mesh->dx(xInd,yInd) * mesh->dz);

                    // J+x
                    BoutReal Jpx = 0.25*( g(xInd+1,yInd, zInd  )*(f_xP1yzP1 - f_xP1yzM1) -
                                          g(xInd-1,yInd, zInd  )*(f_xM1yzP1 - f_xM1yzM1) -
                                          g(xInd  ,yInd, zIndP1)*(f_xP1yzP1 - f_xM1yzP1) +
                                          g(xInd  ,yInd, zIndM1)*(f_xP1yzM1 - f_xM1yzM1))
                      / (mesh->dx(xInd,yInd) * mesh->dz);
                    // Jx+
                    BoutReal Jxp = 0.25*( g(xInd+1, yInd, zIndP1)*(f_xyzP1 - f_xP1yz) -
                                          g(xInd-1, yInd, zIndM1)*(f_xM1yz - f_xyzM1) -
                                          g(xInd-1, yInd, zIndP1)*(f_xyzP1 - f_xM1yz) +
                                          g(xInd+1, yInd, zIndM1)*(f_xP1yz - f_xyzM1))
                      / (mesh->dx(xInd,yInd) * mesh->dz);

                    result(xInd,yInd,zInd) = (Jpp + Jpx + Jxp) / 3.;
                }
            }
        }
        // Special care with last outer point
        int xInd = mesh->xend;
        for(int yInd=mesh->ystart;yInd<=mesh->yend;yInd++){
            for(int zInd=0;zInd<ncz;zInd++) {
                int zIndP1 = (zInd + 1) % ncz;
                int zIndM1 = (zInd - 1 + ncz) % ncz;

                /* Boilerplate
                 * BoutReal f_xyz = pow( (-f(xInd  -1, yInd, zInd  ) + f(xInd  +1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                 */
                int xIndP1 = xInd + 1;
                int xIndM1 = xInd - 1;
                BoutReal f_xyzP1   = pow( (-f(xInd  -1, yInd, zIndP1) + f(xInd  +1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xyzM1   = pow( (-f(xInd  -1, yInd, zIndM1) + f(xInd  +1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xM1yz   = pow( (-f(xIndM1-1, yInd, zInd  ) + f(xIndM1+1, yInd, zInd  ))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xM1yzP1 = pow( (-f(xIndM1-1, yInd, zIndP1) + f(xIndM1+1, yInd, zIndP1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                BoutReal f_xM1yzM1 = pow( (-f(xIndM1-1, yInd, zIndM1) + f(xIndM1+1, yInd, zIndM1))/(2.0*mesh->dx(xInd,yInd)) , 2.0);
                // xIndP1 + 1 does not exist, therefore, we will use a one-sided stencil
                // We use a second order one sided FD stencil (1/2, -2, 3/2)
                YOU ARE HERE: USE DERIVATION D3DX3 TO FIND A ONE-SIDED STENCIL
                    ALSO: MAKE MES ARAKAWA TEST
                BoutReal f_xP1yzP1 = pow( (
                                                  f(xIndP1-2, yInd, zIndP1)
                                           - 4.0* f(xIndP1-1, yInd, zIndP1)
                                           + 3.0* f(xIndP1  , yInd, zIndP1)
                                          )
                                          /(2.0*mesh->dx(xInd,yInd))
                                        , 2.0);
                BoutReal f_xP1yzM1 = pow( (
                                                  f(xIndP1-2, yInd, zIndM1)
                                           - 4.0* f(xIndP1-1, yInd, zIndM1)
                                           + 3.0* f(xIndP1  , yInd, zIndM1)
                                          )/(2.0*mesh->dx(xInd,yInd))
                                        , 2.0);
                BoutReal f_xP1yz   = pow( (
                                                   f(xIndP1-1, yInd, zInd  )
                                            - 4.0* f(xIndP1-2, yInd, zInd  )
                                            + 3.0* f(xIndP1  , yInd, zInd  )
                                           )/(2.0*mesh->dx(xInd,yInd))
                                        , 2.0);

                // J++ = DDZ(f)*DDX(g) - DDX(f)*DDZ(g)
                BoutReal Jpp = 0.25*( (f_xyzP1                 - f_xyzM1               )*
                                      (g(xInd+1, yInd, zInd  ) - g(xInd-1, yInd, zInd  ) )
                                      -
                                      (f_xP1yz                 - f_xM1yz                )*
                                      (g(xInd  , yInd, zIndP1) - g(xInd  , yInd, zIndM1)) )
                  / (mesh->dx(xInd,yInd) * mesh->dz);

                // J+x
                BoutReal Jpx = 0.25*( g(xInd+1,yInd, zInd  )*(f_xP1yzP1 - f_xP1yzM1) -
                                      g(xInd-1,yInd, zInd  )*(f_xM1yzP1 - f_xM1yzM1) -
                                      g(xInd  ,yInd, zIndP1)*(f_xP1yzP1 - f_xM1yzP1) +
                                      g(xInd  ,yInd, zIndM1)*(f_xP1yzM1 - f_xM1yzM1))
                  / (mesh->dx(xInd,yInd) * mesh->dz);
                // Jx+
                BoutReal Jxp = 0.25*( g(xInd+1, yInd, zIndP1)*(f_xyzP1 - f_xP1yz) -
                                      g(xInd-1, yInd, zIndM1)*(f_xM1yz - f_xyzM1) -
                                      g(xInd-1, yInd, zIndP1)*(f_xyzP1 - f_xM1yz) +
                                      g(xInd+1, yInd, zIndM1)*(f_xP1yz - f_xyzM1))
                  / (mesh->dx(xInd,yInd) * mesh->dz);

                result(xInd,yInd,zInd) = (Jpp + Jpx + Jxp) / 3.;
            }
        }
    }

    return result;
}
