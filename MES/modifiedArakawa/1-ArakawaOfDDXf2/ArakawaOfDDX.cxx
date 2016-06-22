// *************** Simulation of ArakawaOfDDX *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "ArakawaOfDDX.hxx"

// Initialization of the physics
// ############################################################################
int ArakawaOfDDX::init(bool restarting) {
    TRACE("Halt in ArakawaOfDDXf2::init");

    /* NOTE: Calls createOperators without making an object of OwnOperators.
     *       The child is typecasted to the parent
     */
    ownOp = OwnOperators::createOperators();

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
    S_num = ownOp->ArakawaOfDDXPhi2N(phi, n);

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
