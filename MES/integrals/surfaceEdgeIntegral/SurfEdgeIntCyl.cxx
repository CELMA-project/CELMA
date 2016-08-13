// *************** Simulation of SurfEdgeIntCyl *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "SurfEdgeIntCyl.hxx"

// Initialization of the physics
// ############################################################################
int SurfEdgeIntCyl::init(bool restarting) {
    TRACE("Halt in SurfEdgeIntCyl::init");

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
    Options *SOpt = options->getSection("S");
    SOpt->get("solution", S, 0.0);
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

    // Calculate the integral
    volumeIntegral(f, S_num);

    // Error in S
    e = S_num - S;

    // Save the variables
    SAVE_ONCE2(Lx, Ly);
    SAVE_ONCE3(e, S, S_num);
    SAVE_ONCE (f);

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
int SurfEdgeIntCyl::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(SurfEdgeIntCyl);

// Destructor
SurfEdgeIntCyl::~SurfEdgeIntCyl(){
    TRACE("SurfEdgeIntCyl::~SurfEdgeIntCyl");
}
