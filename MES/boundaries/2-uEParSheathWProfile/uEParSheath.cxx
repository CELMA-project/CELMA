// *************** Simulation of uEParSheath *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "uEParSheath.hxx"

// Initialization of the physics
// ############################################################################
int UeSheath::init(bool restarting) {
  TRACE("Halt in UeSheath::init");

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
  // uEParOrigin
  uEParOrigin = FieldFactory::get()->create3D(
      "uEPar:function", Options::getRoot(), mesh, CELL_CENTRE, 0);
  // phi
  phi = FieldFactory::get()->create3D("phi:function", Options::getRoot(), mesh,
                                      CELL_CENTRE, 0);
  // profile
  profile = FieldFactory::get()->create3D(
      "dampProf:profile", Options::getRoot(), mesh, CELL_CENTRE, 0);
  // ************************************************************************

  // Add a FieldGroup to communicate
  // ************************************************************************
  // Only these fields will be taken derivatives of
  com_group.add(uEParOrigin);
  // ************************************************************************

  // Communicate before taking derivatives
  mesh->communicate(com_group);

  output << "\n\n\n\n\n\n\nNow running test" << std::endl;

  // Copy
  uEParWBC = copy(uEParOrigin);

  // Extrapolate
  ownBC.uEParSheath(uEParWBC, phi, Lambda, phiRef, profile);

  // Error in S
  e = uEParWBC - uEParOrigin;

  // Save the variables
  SAVE_ONCE2(Lx, Ly);
  SAVE_ONCE2(uEParWBC, uEParOrigin);
  SAVE_ONCE2(phi, profile);
  SAVE_ONCE(e);

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
int UeSheath::rhs(BoutReal t) { return 0; }
// ############################################################################

// Create a simple main() function
BOUTMAIN(UeSheath);

// Destructor
UeSheath::~UeSheath() { TRACE("UeSheath::~UeSheath"); }
