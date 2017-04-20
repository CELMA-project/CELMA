// *************** Simulation of DDXCylinderHighResCenter *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "DDXCylinderHighResCenter.hxx"

Field3D DDXCylHighResCenter(Field3D const &f);

// Initialization of the physics
// ############################################################################
int DDXCylinderHighResCenter::init(bool restarting) {
  TRACE("Halt in DDXCylinderHighResCenter::init");

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
int DDXCylinderHighResCenter::rhs(BoutReal t) {

  TRACE("DDXCylHighResCenter::rhs");

  // Calculate
  S_num = DDXCylHighResCenter(f);

  // Error in S
  e = S_num - S;
  return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(DDXCylinderHighResCenter);

// Destructor
DDXCylinderHighResCenter::~DDXCylinderHighResCenter() {
  TRACE("DDXCylinderHighResCenter::~DDXCylinderHighResCenter");
}

// New implementation
Field3D DDXCylHighResCenter(Field3D const &f) {
  TRACE("DDXCylHighResCenter");

  Field3D result;

  // Calculate points away from the center
  result = DDX(f);

  /* NOTE: The index corresponding to pi
   *       Since z in [0, 2 pi[, the z index corresponding to pi is
   *       mesh->LocalNz/ 2
   */
  int piIndex = mesh->LocalNz / 2;

  if (mesh->firstX()) {
    // For all z indices corresponding to a theta angle below pi
    int xInd = mesh->xstart;
    for (int yInd = mesh->ystart; yInd <= mesh->yend; yInd++) {
      for (int zInd = 0; zInd < piIndex; zInd++) {
        // First point
        // Use a seven point stencil
        result(xInd, yInd, zInd) =
            (
                // Third "ghost"
                -(1.0 / 60.0) * f(xInd + 2, yInd, zInd + piIndex)
                // Second "ghost"
                + (3.0 / 20.0) * f(xInd + 1, yInd, zInd + piIndex)
                // First ghost (set ownBC.innerRhoCylinder(f))
                - (3.0 / 4.0) * f(xInd - 1, yInd, zInd) +
                (3.0 / 4.0) * f(xInd + 1, yInd, zInd) -
                (3.0 / 20.0) * f(xInd + 2, yInd, zInd) +
                (1.0 / 60.0) * f(xInd + 3, yInd, zInd)) /
            mesh->coordinates()->dx(xInd, yInd);
        // Second point
        // Use a seven point stencil
        result(xInd + 1, yInd, zInd) =
            (
                // Second "ghost"
                -(1.0 / 60.0) * f(xInd + 1, yInd, zInd + piIndex)
                // First ghost (set ownBC.innerRhoCylinder(f))
                + (3.0 / 20.0) * f(xInd - 1, yInd, zInd) -
                (3.0 / 4.0) * f(xInd, yInd, zInd) +
                (3.0 / 4.0) * f(xInd + 2, yInd, zInd) -
                (3.0 / 20.0) * f(xInd + 3, yInd, zInd) +
                (1.0 / 60.0) * f(xInd + 4, yInd, zInd)) /
            mesh->coordinates()->dx(xInd + 1, yInd);
      }
    }
    // For all z indices corresponding to a theta value including and above
    // pi
    for (int yInd = mesh->ystart; yInd <= mesh->yend; yInd++) {
      for (int zInd = piIndex; zInd < mesh->LocalNz; zInd++) {
        // First point
        // Use a seven point stencil
        result(xInd, yInd, zInd) =
            (
                // Third "ghost"
                -(1.0 / 60.0) * f(xInd + 2, yInd, zInd - piIndex)
                // Second "ghost"
                + (3.0 / 20.0) * f(xInd + 1, yInd, zInd - piIndex)
                // First ghost (set ownBC.innerRhoCylinder(f))
                - (3.0 / 4.0) * f(xInd - 1, yInd, zInd) +
                (3.0 / 4.0) * f(xInd + 1, yInd, zInd) -
                (3.0 / 20.0) * f(xInd + 2, yInd, zInd) +
                (1.0 / 60.0) * f(xInd + 3, yInd, zInd)) /
            mesh->coordinates()->dx(xInd, yInd);
        // Second point
        // Use a seven point stencil
        result(xInd + 1, yInd, zInd) =
            (
                // Second "ghost"
                -(1.0 / 60.0) * f(xInd + 1, yInd, zInd - piIndex)
                // First ghost (set ownBC.innerRhoCylinder(f))
                + (3.0 / 20.0) * f(xInd - 1, yInd, zInd) -
                (3.0 / 4.0) * f(xInd, yInd, zInd) +
                (3.0 / 4.0) * f(xInd + 2, yInd, zInd) -
                (3.0 / 20.0) * f(xInd + 3, yInd, zInd) +
                (1.0 / 60.0) * f(xInd + 4, yInd, zInd)) /
            mesh->coordinates()->dx(xInd + 1, yInd);
      }
    }
  }

  return result;
}
