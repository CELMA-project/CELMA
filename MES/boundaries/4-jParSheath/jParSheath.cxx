// *************** Simulation of jParSheath *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "jParSheath.hxx"

// Initialization of the physics
// ############################################################################
int JParSheath::init(bool restarting) {
  TRACE("Halt in JParSheath::init");

  // Get the option (before any sections) in the BOUT.inp file
  Options *options = Options::getRoot();

  // Load from the geometry
  // ************************************************************************
  Options *geom = options->getSection("geom");
  geom->get("Lx", Lx, 0.0);
  geom->get("Ly", Ly, 0.0);
  // ************************************************************************

  // Load from the constants
  // ************************************************************************
  Options *cst = options->getSection("cst");
  cst->get("Lambda", Lambda, 0.0);
  cst->get("phiRef", phiRef, 0.0);
  // ************************************************************************

  // Obtain the fields
  // ************************************************************************
  // uEPar
  uEPar = FieldFactory::get()->create3D("uEPar:function", Options::getRoot(),
                                        mesh, CELL_CENTRE, 0);
  // uIPar
  uIPar = FieldFactory::get()->create3D("uIPar:function", Options::getRoot(),
                                        mesh, CELL_CENTRE, 0);
  // phi
  phi = FieldFactory::get()->create3D("phi:function", Options::getRoot(), mesh,
                                      CELL_CENTRE, 0);
  // n
  n = FieldFactory::get()->create3D("n:function", Options::getRoot(), mesh,
                                    CELL_CENTRE, 0);
  // jParOrigin
  jParOrigin = FieldFactory::get()->create3D(
      "jPar:function", Options::getRoot(), mesh, CELL_CENTRE, 0);
  // ************************************************************************

  // Add a FieldGroup to communicate
  // ************************************************************************
  // Only these fields will be taken derivatives of
  com_group.add(uEPar);
  // ************************************************************************

  // Communicate before taking derivatives
  mesh->communicate(com_group);

  jParWBC = copy(jParOrigin);

  // Save the variables
  SAVE_ONCE2(Lx, Ly);
  SAVE_ONCE2(jParWBC, jParOrigin);
  SAVE_ONCE4(uIPar, uEPar, phi, n);
  SAVE_ONCE(e);

  return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int JParSheath::rhs(BoutReal t) {

  TRACE("JParSheath::rhs");

  // Extrapolate
  ownBC.jParSheath(jParWBC, uEPar, uIPar, phi, n, Lambda, phiRef);

  // Error in S
  e = jParWBC - jParOrigin;

  return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(JParSheath);

// Destructor
JParSheath::~JParSheath() { TRACE("JParSheath::~JParSheath"); }
