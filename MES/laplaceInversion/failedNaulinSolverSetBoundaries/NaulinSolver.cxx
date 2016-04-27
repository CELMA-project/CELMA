// *************** Simulation of NaulinSolver *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "NaulinSolver.hxx"
#include "cxxRoutines/ownImplimentations.cxx"

// Initialization of the physics
// ############################################################################
int NaulinSolver::init(bool restarting) {
    TRACE("Halt in NaulinSolver::init");

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Create the laplace object
    // ************************************************************************
    // The laplace object will look in the section phiSolver in the BOUT.inp
    // file
    Options *phiSol_opt = options->getSection("phiSolver");
    phiSolver = Laplacian::create(phiSol_opt);
    // Set the coefficients manually (should also be set to this by default)
    phiSolver->setCoefD(1.0);
    phiSolver->setCoefC(1.0);
    phiSolver->setCoefA(0.0);
    // Get the tolerances
    phiSol_opt->get("atol",  atol,  1.0e-10);
    phiSol_opt->get("rtol",  rtol,  1.0e-5);
    phiSol_opt->get("maxit", maxit, 300);
    // Check if the Naulin iterator should be monitored
    phiSol_opt->get("NaulinMonitor", NaulinMonitor, true);
    // ************************************************************************

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

    innerRhoCylinder(n);
    innerRhoCylinder(phi);    // Used to set BC in lapalace inversion
    // ************************************************************************

    // Calculate ln_n
    ln_n = log(n);

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
    phi_num = NSolver(atol, rtol, maxit,
                      ln_n, n, vortD,
                      phi, vort,
                      NaulinMonitor);

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

// Destructor
NaulinSolver::~NaulinSolver(){
    TRACE("NaulinSolver::~NaulinSolver");

    delete phiSolver;
}
