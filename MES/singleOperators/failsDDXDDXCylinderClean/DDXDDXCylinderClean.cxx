// *************** Simulation of DDXDDXCylinderClean *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "DDXDDXCylinderClean.hxx"

// Set xout ghost point of the resulting field of DDX
// *****************************************************************************
void xout_after_DDX(Field3D &out, Field3D const &in){
    /* Info:
     * Set xout ghost point of the resulting field of DDX
     *
     * Input:
     * in    - The field that was used to find out
     * out   - DDX(in), but without the ghost point in xout
     *
     * Output:
     * out- DDi(in), with the ghost point of xout set
     */
    TRACE("Halt in xout_after_DDX");

    if (mesh->lastX()){
        /* NOTE: xend
         *       xend = index value of last inner point on this processor
         *       xend+1 = first guard point
         */
        int x_ind = mesh->xend + 1;

        // Newton polynomial of fourth order (including boundary) evaluated at ghost
        for(int y_ind = mesh->ystart; y_ind <= mesh->yend; y_ind++){
            for(int z_ind = 0; z_ind < mesh->ngz -1; z_ind ++){
                out(x_ind, y_ind, z_ind) =
                    - (1.0/5.0)*out(x_ind-3, y_ind, z_ind)
                    +           out(x_ind-2, y_ind, z_ind)
                    - 3.0*      out(x_ind-1, y_ind, z_ind)
                    + (16.0/5.0)*
                      // Calculation of boundary value
                      ((in(x_ind, y_ind, z_ind)-in(x_ind-1, y_ind, z_ind))/
                       mesh->dx(x_ind, y_ind));
            }
        }
    }
}
// *****************************************************************************

// Initialization of the physics
// ############################################################################
int DDXDDXCylinderClean::init(bool restarting) {
    TRACE("Halt in DDXDDXCylinderClean::init");

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
    DDXf = DDX(f);
    ownBC.innerRhoCylinder(DDXf); // Set inner rho boundary
    xout_after_DDX(DDXf, f);      // Set outer rho boundary
    mesh->communicate(DDXf);
    S_num = DDX(DDXf);

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
int DDXDDXCylinderClean::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(DDXDDXCylinderClean);

// Destructor
DDXDDXCylinderClean::~DDXDDXCylinderClean(){
    TRACE("DDXDDXCylinderClean::~DDXDDXCylinderClean");
}
