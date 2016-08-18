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

    // Load from the geometry
    // ************************************************************************
    Options *switches = options->getSection("switches");
    switches->get("saveFields", saveFields, false);
    // ************************************************************************

    // Obtain the fields
    // ************************************************************************
    // See initialprofiles.cxx for creation of vectors
    // v
    v.covariant = false; // Referring to the components
    initial_profile("v", v);

    // S
    Options *SOpt = options->getSection("S");
    SOpt->get("S_Xout" , S_Xout , 0.0);
    SOpt->get("S_Yup"  , S_Yup  , 0.0);
    SOpt->get("S_Ydown", S_Ydown, 0.0);
    // Add them to the result vector
    // ************************************************************************

    // Add a FieldGroup to communicate
    // ************************************************************************
    // Only these fields will be taken derivatives of
    com_group.add(v);
    // ************************************************************************

    // Set boundaries manually
    // ************************************************************************
    v.x.setBoundary("vx");
    v.y.setBoundary("vy");
    v.z.setBoundary("vz");
    v.applyBoundary();
    ownBC.innerRhoCylinder(v.x);
    ownBC.innerRhoCylinder(v.y);
    ownBC.innerRhoCylinder(v.z);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(com_group);

    output << "\n\n\n\n\n\n\nNow running test" << std::endl;

    // Calculate the integral
    helper.surfaceEdgeIntegral(v, results);
    S_Xout_num  = results[1];
    S_Ydown_num = results[2];
    S_Yup_num   = results[3];

    // Error in S
    e_Xout  = S_Xout_num  - S_Xout;
    e_Ydown = S_Ydown_num - S_Ydown;
    e_Yup   = S_Yup_num   - S_Yup;

    // Save the variables
    SAVE_ONCE2(Lx, Ly);
    SAVE_ONCE3(S_Xout_num, S_Ydown_num, S_Yup_num);
    SAVE_ONCE3(S_Xout,     S_Ydown,     S_Yup);
    SAVE_ONCE3(e_Xout,     e_Ydown,     e_Yup);

    if(saveFields){
        SAVE_ONCE(v);
    }

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
