// *************** Simulation of NaulinSolver *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "NaulinSolver.hxx"

// Initialization of the physics
// ############################################################################
int NaulinSolver::init(bool restarting) {
    TRACE("Halt in NaulinSolver::init");

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Load from the geometry
    // ************************************************************************
    Options *geom = options->getSection("geom");
    geom->get("Lx", Lx, 0.0);
    // ************************************************************************

    // Create the solver
    // ************************************************************************
    ownLapl.create(ownOp, ownBC);
    // ************************************************************************

    // Obtain the fields
    // ************************************************************************
    // n
    n = FieldFactory::get()
             ->create3D("n:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // phi
    phi = FieldFactory::get()
             ->create3D("phi:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // vortD
    vortD = FieldFactory::get()
            ->create3D("vortD:function", Options::getRoot(), mesh, CELL_CENTRE, 0);
    // ************************************************************************

    // Set boundaries manually
    // ************************************************************************
    phi.setBoundary("phi");
    n.  setBoundary("n");

    phi.applyBoundary();
    n.  applyBoundary();

    ownBC.innerRhoCylinder(n);
    ownBC.innerRhoCylinder(phi);    // Used to set BC in lapalace inversion
    // ************************************************************************

    // Preparations
    // ************************************************************************
    ln_n = log(n);
    mesh->communicate(ln_n);
    gradPerp_ln_n = ownOp.Grad_perp(ln_n);
    // ************************************************************************

    // Add a FieldGroup to communicate
    // ************************************************************************
    // Only these fields will be taken derivatives of
    com_group.add(ln_n);
    com_group.add(phi);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(com_group);

    output << "\n\n\n\n\n\n\nNow running test" << std::endl;

    // Call the solver
    // phi_num and vort are output
    phi_num = ownLapl.NaulinSolver(gradPerp_ln_n, n, vortD, phi, vort);

    // Error in phi
    e = phi - phi_num;

    // Save the variables
    SAVE_ONCE (Lx);
    SAVE_ONCE2(vortD, vort);
    SAVE_ONCE (n);
    SAVE_ONCE2(phi, phi_num);
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
int NaulinSolver::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(NaulinSolver);
