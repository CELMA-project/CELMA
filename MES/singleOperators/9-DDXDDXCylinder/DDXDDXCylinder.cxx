// *************** Simulation of DDXDDXCylinder *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "DDXDDXCylinder.hxx"

// Set xout ghost point of the resulting field of DDX
// *****************************************************************************
void xOutAfterDDX(Field3D &out, Field3D const &in){
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
    TRACE("Halt in xOutAfterDDX");

    if (mesh->lastX()){
        /* NOTE: xend
         *       xend = index value of last inner point on this processor
         *       xend+1 = first guard point
         */
        int xInd = mesh->xend + 1;

        // Newton polynomial of fourth order evaluated at ghost
        for(int yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
            for(int zInd = 0; zInd < mesh->ngz -1; zInd ++){
                out(xInd, yInd, zInd) =
                    -      out(xInd-4, yInd, zInd)
                    +  4.0*out(xInd-3, yInd, zInd)
                    -  6.0*out(xInd-2, yInd, zInd)
                    +  4.0*out(xInd-1, yInd, zInd)
                      ;
            }
        }
    }
}
// *****************************************************************************

// Initialization of the physics
// ############################################################################
int DDXDDXCylinder::init(bool restarting) {
    TRACE("Halt in DDXDDXCylinder::init");

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
    xOutAfterDDX(DDXf, f);      // Set outer rho boundary
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
int DDXDDXCylinder::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(DDXDDXCylinder);

// Destructor
DDXDDXCylinder::~DDXDDXCylinder(){
    TRACE("DDXDDXCylinder::~DDXDDXCylinder");
}
