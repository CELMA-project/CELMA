// *************** Simulation of DivExBAdv *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "divExBAdv.hxx"

// Initialization of the physics
// ############################################################################
int DivExBAdv::init(bool restarting) {
    TRACE("Halt in DivExBAdv::init");

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Load from the geometry
    // ************************************************************************
    Options *geom = options->getSection("geom");
    geom->get("Lx", Lx, 0.0);
    // ************************************************************************

    // Obtain the fields
    // ************************************************************************
    // n
    n = FieldFactory::get()
             ->create3D("n:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // phi
    phi = FieldFactory::get()
             ->create3D("phi:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // The source
    S = FieldFactory::get()
        ->create3D("S:solution", Options::getRoot(), mesh, CELL_CENTRE, 0);
    // ************************************************************************

    // Add a FieldGroup to communicate
    // ************************************************************************
    // Only these fields will be taken derivatives of
    com_group.add(n);
    com_group.add(phi);
    // ************************************************************************

    // Set boundaries manually
    // ************************************************************************
    phi.setBoundary("phi");
    n.setBoundary("n");

    phi.applyBoundary();
    n.applyBoundary();

    ownBC.innerRhoCylinder(phi);
    ownBC.innerRhoCylinder(n);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(com_group);

    output << "\n\n\n\n\n\n\nNow running test" << std::endl;

    // Calculate
    S_num = ownOp.div_uE_dot_grad_n_GradPerp_phi(n, phi);

    // Error in phi
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
int DivExBAdv::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(DivExBAdv);

// Destructor
DivExBAdv::~DivExBAdv(){
    TRACE("DivExBAdv::~DivExBAdv");
}
