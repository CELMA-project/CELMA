// *************** Simulation of DivPerp *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "divPerp.hxx"

// Initialization of the physics
// ############################################################################
int DivPerp::init(bool restarting) {
    TRACE("Halt in DivPerp::init");

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
    f.x = FieldFactory::get()
             ->create3D("fx:function", Options::getRoot(), mesh, CELL_CENTRE, 0);
    f.y = FieldFactory::get()
             ->create3D("fy:function", Options::getRoot(), mesh, CELL_CENTRE, 0);
    f.z = FieldFactory::get()
             ->create3D("fz:function", Options::getRoot(), mesh, CELL_CENTRE, 0);
    f.covariant = false;

    // The source
    S = FieldFactory::get()
        ->create3D("S:Solution", Options::getRoot(), mesh, CELL_CENTRE, 0);
    // ************************************************************************

    // Add a FieldGroup to communicate
    // ************************************************************************
    // Only these fields will be taken derivatives of
    com_group.add(f);
    // ************************************************************************

    // Set boundaries manually
    // ************************************************************************
    f.x.setBoundary("fx");
    f.y.setBoundary("fy");
    f.z.setBoundary("fz");

    f.x.applyBoundary();
    f.y.applyBoundary();
    f.z.applyBoundary();

    ownBC.innerRhoCylinder(f.x);
    ownBC.innerRhoCylinder(f.y);
    ownBC.innerRhoCylinder(f.z);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(com_group);

    output << "\n\n\n\n\n\n\nNow running test" << std::endl;

    // Calculate
    S_num = Div(f);

    // Error in phi
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
int DivPerp::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(DivPerp);

// Destructor
DivPerp::~DivPerp(){
    TRACE("DivPerp::~DivPerp");
}
