// *************** Simulation of cauchyBC *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "cauchyBC.hxx"

// Initialization of the physics
// ############################################################################
int CauchyBC::init(bool restarting) {
  TRACE("Halt in CauchyBC::init");

  // Get the option (before any sections) in the BOUT.inp file
  Options *options = Options::getRoot();

  // Load from the geometry
  // ************************************************************************
  Options *geom = options->getSection("geom");
  geom->get("Lx", Lx, 0.0);
  geom->get("Ly", Ly, 0.0);
  // ************************************************************************

  // Obtain the fields
  // ************************************************************************
  // fOrigin
  fOrigin = FieldFactory::get()->create3D("f:function", Options::getRoot(),
                                          mesh, CELL_CENTRE, 0);
  // ************************************************************************

  // Add a FieldGroup to communicate
  // ************************************************************************
  // Only these fields will be taken derivatives of
  com_group.add(fOrigin);
  // ************************************************************************

  // Communicate before taking derivatives
  mesh->communicate(com_group);

  // Copy ensures that the two doesn't share memory
  fCauchy = copy(fOrigin);

  // Save the variables
  SAVE_ONCE2(Lx, Ly);
  SAVE_ONCE2(fCauchy, fOrigin);
  SAVE_ONCE(e);

  return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int CauchyBC::rhs(BoutReal t) {

  TRACE("CauchyBC::rhs");

  // Prepare cauchy
  ownBC.prepareCauchy("f");
  ownBC.cauchyYDown(fCauchy);

  // Error in S
  e = fCauchy - fOrigin;

  return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(CauchyBC);

// Destructor
CauchyBC::~CauchyBC() { TRACE("CauchyBC::~CauchyBC"); }
