// *************** Simulation of D2DZ2Cylinder *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "D2DZ2Cylinder.hxx"

// Initialization of the physics
// ############################################################################
int D2DZ2Cylinder::init(bool restarting) {
  TRACE("Halt in D2DZ2Cylinder::init");

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
  f = FieldFactory::get()->create3D("f:function", Options::getRoot(), mesh,
                                    CELL_CENTRE, 0);

  // S
  S = FieldFactory::get()->create3D("S:solution", Options::getRoot(), mesh,
                                    CELL_CENTRE, 0);
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

  // Save the variables
  SAVE_ONCE(Lx);
  SAVE_ONCE3(f, S, S_num);
  SAVE_ONCE(e);

  return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int D2DZ2Cylinder::rhs(BoutReal t) {
  TRACE("D2DZ2Cylinder::~D2DZ2Cylinder");

  // Calculate
  S_num = D2DZ2(f);

  // Error in S
  e = S_num - S;
  return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(D2DZ2Cylinder);

// Destructor
D2DZ2Cylinder::~D2DZ2Cylinder() { TRACE("D2DZ2Cylinder::~D2DZ2Cylinder"); }
