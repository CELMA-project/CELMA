// *************** Simulation of JTimesDivSource *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "JTimesDivSource.hxx"

// Initialization of the physics
// ############################################################################
int JTimesDivSource::init(bool restarting) {
    TRACE("Halt in JTimesDivSource::init");

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Load from the geometry
    // ************************************************************************
    Options *geom = options->getSection("geom");
    geom->get("Lx", Lx, 0.0);
    // ************************************************************************

    // Obtain the fields
    // ************************************************************************
    // S_n
    S_n = FieldFactory::get()
             ->create3D("S_n:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

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
    com_group.add(S_n);
    com_group.add(phi);
    // ************************************************************************

    // Set boundaries manually
    // ************************************************************************
    S_n.setBoundary("S_n");
    phi.setBoundary("phi");

    S_n.applyBoundary();
    phi.applyBoundary();

    ownBC.innerRhoCylinder(S_n);
    ownBC.innerRhoCylinder(phi);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(com_group);

    output << "\n\n\n\n\n\n\nNow running test" << std::endl;

    // Calculate
    S_num = mesh->J*ownOp->div_f_GradPerp_g(S_n, phi);

    // Error in phi
    e = S_num - S;

    // Save the variables
    SAVE_ONCE (Lx);
    SAVE_ONCE4(phi, S_n, S, S_num);
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
int JTimesDivSource::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(JTimesDivSource);

// Destructor
JTimesDivSource::~JTimesDivSource(){
    TRACE("JTimesDivSource::~JTimesDivSource");
}
