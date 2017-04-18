// *************** Simulation of RadialLowPass *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "radialLowPass.hxx"

// Initialization of the physics
// ############################################################################
int RadialLowPass::init(bool restarting) {
  TRACE("Halt in RadialLowPass::init");

  // Get the option (before any sections) in the BOUT.inp file
  Options *options = Options::getRoot();

  // Load from the geometry
  // ************************************************************************
  Options *geom = options->getSection("geom");
  geom->get("Lx", Lx, 0.0);
  // ************************************************************************

  // Create the filter
  // ************************************************************************
  ownFilter = OwnFilters::createFilter();
  // ************************************************************************

  // Allocate the fields
  // ************************************************************************
  unfiltered = 0.0;
  // ************************************************************************

  // Create the unfiltered
  // ************************************************************************
  /* NOTE: The Nyquist mode
   *       As the offset mode is the first mode on the positive frequencies,
   *       the Nyquist mode is the first mode on the negative frequencies
   */
  nyquistMode = int((mesh->LocalNz - 1) / 2.0);

  for (int mode = 1; mode <= nyquistMode; mode++) {
    for (int xInd = 0; xInd < mesh->LocalNx; xInd++) {
      for (int yInd = 0; yInd < mesh->LocalNy; yInd++) {
        for (int zInd = 0; zInd < mesh->LocalNz - 1; zInd++) {
          unfiltered(xInd, yInd, zInd) +=
              cos((TWOPI * mode / TWOPI) * zInd * (mesh->coordinates()->dz));
        }
      }
    }
  }
  // Add offset
  unfiltered += 2;
  // ************************************************************************

  output << "\n\n\n\n\n\n\nNow filtering" << std::endl;

  // Filter
  // ************************************************************************
  filtered = ownFilter->ownFilter(unfiltered);
  // ************************************************************************

  // Save the variables
  SAVE_ONCE(Lx);
  SAVE_ONCE2(unfiltered, filtered);

  // Finalize
  dump.write();
  dump.close();

  output << "\nFinished, now quitting\n\n\n\n\n\n" << std::endl;

  // Wait for all processors to write data
  MPI_Barrier(BoutComm::get());

  return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int RadialLowPass::rhs(BoutReal t) { return 0; }
// ############################################################################

// Create a simple main() function
BOUTMAIN(RadialLowPass);

// Destructor
RadialLowPass::~RadialLowPass() { TRACE("RadialLowPass::~RadialLowPass"); }
