// *************** Simulation of PolAvgTest *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __PolAvgTest_H__
#define __PolAvgTest_H__

#include <bout/constants.hxx> // Gives PI and TWOPI
#include <bout/physicsmodel.hxx>
#include <field_factory.hxx> // Gives field factory
// Gives own boundaries (doing so by setting ghost points)
#include "../../common/BOUTExtensions/include/ownBCs.hxx"
// Give the integrators
#include "../../common/BOUTExtensions/include/helpers.hxx"

class PolAvgTest : public PhysicsModel {
public:
  // Destructor
  ~PolAvgTest();

protected:
  int init(bool restarting);
  int rhs(BoutReal t);

private:
  // Global variable initialization
  // ############################################################################
  // Variables
  // *****************************************************************************
  Field3D f, S, S_num, e;
  // *****************************************************************************

  // Constants
  // *****************************************************************************
  BoutReal Lx, Ly; // The box dimensions
  // *****************************************************************************

  // Switches
  // *****************************************************************************
  bool saveFields;
  // *****************************************************************************

  // Make a field group to communicate
  // *****************************************************************************
  FieldGroup com_group;
  // *****************************************************************************

  // Other objects
  // *****************************************************************************
  VolumeIntegral volInt; // Helper class
  PolAvg avg;            // Helper class
  // *****************************************************************************
  // ############################################################################
};

#endif
