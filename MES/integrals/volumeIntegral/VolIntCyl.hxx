// *************** Simulation of VolIntCyl *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef MES_INTEGRALS_VOLUMEINTEGRAL_VOLINTCYL_HXX_
#define MES_INTEGRALS_VOLUMEINTEGRAL_VOLINTCYL_HXX_

#include <bout/constants.hxx> // Gives PI and TWOPI
#include <bout/physicsmodel.hxx>
#include <field_factory.hxx> // Gives field factory
// Gives own boundaries (doing so by setting ghost points)
#include "../../../common/BOUTExtensions/include/ownBCs.hxx"
// Give the integrators
#include "../../../common/BOUTExtensions/include/helpers.hxx"

class VolIntCyl : public PhysicsModel {
public:
  // Destructor
  ~VolIntCyl();

protected:
  int init(bool restarting);
  int rhs(BoutReal t);

private:
  // Global variable initialization
  // ############################################################################
  // Variables
  // *****************************************************************************
  Field3D f;
  BoutReal S, S_num;
  BoutReal e;
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
  OwnBCs ownBC;          // Class containing methods which sets the ghost points
  VolumeIntegral volInt; // Helper class
  // *****************************************************************************
  // ############################################################################
};

#endif // MES_INTEGRALS_VOLUMEINTEGRAL_VOLINTCYL_HXX_
