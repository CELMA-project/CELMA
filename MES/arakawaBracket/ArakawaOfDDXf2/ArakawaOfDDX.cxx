// *************** Simulation of ArakawaOfDDX *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "ArakawaOfDDX.hxx"

// Initialization of the physics
// ############################################################################
int ArakawaOfDDX::init(bool restarting) {
  TRACE("Halt in ArakawaOfDDX::init");

  // Get the option (before any sections) in the BOUT.inp file
  Options *options = Options::getRoot();

  // Load from the geometry
  // ************************************************************************
  Options *geom = options->getSection("geom");
  geom->get("Lx", Lx, 0.0);
  // ************************************************************************

  // Obtain the fields
  // ************************************************************************
  // phi
  phi = FieldFactory::get()->create3D("phi:function", Options::getRoot(), mesh,
                                      CELL_CENTRE, 0);

  // n
  n = FieldFactory::get()->create3D("n:function", Options::getRoot(), mesh,
                                    CELL_CENTRE, 0);

  // S
  S = FieldFactory::get()->create3D("S:solution", Options::getRoot(), mesh,
                                    CELL_CENTRE, 0);
  // ************************************************************************

  // Add a FieldGroup to communicate
  // ************************************************************************
  // Only these fields will be taken derivatives of
  com_group.add(phi);
  com_group.add(n);
  // ************************************************************************

  // Set boundaries manually
  // ************************************************************************
  phi.setBoundary("phi");
  phi.applyBoundary();
  ownBC.innerRhoCylinder(phi);
  n.setBoundary("n");
  n.applyBoundary();
  ownBC.innerRhoCylinder(n);
  // ************************************************************************

  // Communicate before taking derivatives
  mesh->communicate(com_group);

  output << "\n\n\n\n\n\n\nNow running test" << std::endl;

  // Calculate the derivative of phi
  DDXPhi = DDX(phi);
  // Reset inner boundary
  ownBC.innerRhoCylinder(DDXPhi);
  // Reset outer boundary
  if (mesh->lastX()) {
    /* NOTE: xend
     *       xend = index value of last inner point on this processor
     *       xend+1 = first guard point
     */
    int ghostIndX = mesh->xend + 1;
    // Newton polynomial of fourth order (including boundary) evaluated at ghost
    for (int yInd = mesh->ystart; yInd <= mesh->yend; yInd++) {
      for (int zInd = 0; zInd < mesh->LocalNz; zInd++) {
        DDXPhi(ghostIndX, yInd, zInd) =
            -DDXPhi(ghostIndX - 4, yInd, zInd) +
            4.0 * DDXPhi(ghostIndX - 3, yInd, zInd) -
            6.0 * DDXPhi(ghostIndX - 2, yInd, zInd) +
            4.0 * DDXPhi(ghostIndX - 1, yInd, zInd);
      }
    }
  }

  // Communicate before taking new derivative
  mesh->communicate(DDXPhi);

  // Save the variables
  SAVE_ONCE(Lx);
  SAVE_ONCE4(phi, n, S, S_num);
  SAVE_ONCE(e);

  return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int ArakawaOfDDX::rhs(BoutReal t) {

TRACE("ArakawaOfDDX::rhs");

  // Calculate
  S_num = bracket(pow(DDXPhi, 2.0), n, bm);

  // Error in S
  e = S_num - S;

    return 0; }
// ############################################################################

// Create a simple main() function
BOUTMAIN(ArakawaOfDDX);

// Destructor
ArakawaOfDDX::~ArakawaOfDDX() { TRACE("ArakawaOfDDX::~ArakawaOfDDXf2"); }
