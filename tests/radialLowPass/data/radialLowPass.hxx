// *************** Simulation of RadialLowPass *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __RadialLowPass_H__
#define __RadialLowPass_H__

#include "../../common/c/include/ownFilters.hxx"
#include <bout/constants.hxx> // Gives PI and TWOPI
#include <bout/physicsmodel.hxx>
#include <field_factory.hxx> // Gives field factory

class RadialLowPass : public PhysicsModel {
public:
  // Destructor
  ~RadialLowPass();

protected:
  int init(bool restarting);
  int rhs(BoutReal t);

private:
  // Global variable initialization
  // ############################################################################
  // Variables
  // *****************************************************************************
  Field3D unfiltered, filtered;
  // *****************************************************************************

  // Constants
  // *****************************************************************************
  BoutReal Lx;     // The box dimensions
  int nyquistMode; // The nyquist mode
  // *****************************************************************************

  // Other objects
  // *****************************************************************************
  OwnFilters *ownFilter; // Pointer to the filter class
  // *****************************************************************************
  // ############################################################################
};

#endif
