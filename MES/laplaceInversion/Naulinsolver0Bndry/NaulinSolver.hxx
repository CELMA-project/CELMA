// *************** Simulation of NaulinSolver *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __NaulinSolver_H__
#define __NaulinSolver_H__

#include <bout/constants.hxx> // Gives PI and TWOPI
#include <bout/physicsmodel.hxx>
#include <derivs.hxx>         // Gives the bracket method
#include <difops.hxx>         // Gives the diff options
#include <field_factory.hxx>  // Gives field factory
#include <invert_laplace.hxx> // Gives invert laplace option
#include <vecops.hxx>         // Gives the vec diff options
// Gives own boundaries (doing so by setting ghost points)
#include "../../../common/c/include/ownBCs.hxx"
// Gives own operators
#include "../../common/c/include/ownOperators.hxx"
// Gives own laplacian inversions
#include "../../common/c/include/ownLaplacianInversions.hxx"

class NaulinSolver : public PhysicsModel {
protected:
  int init(bool restarting);
  int rhs(BoutReal t);

private:
  // Global variable initialization
  // ############################################################################
  // Variables
  // *****************************************************************************
  Field3D n, phi, vortD;
  Field3D phi_num, e;
  Field3D ln_n, vort;
  Vector3D gradPerp_ln_n;
  // *****************************************************************************

  // Constants
  // *****************************************************************************
  BoutReal Lx; // The box dimensions
  // *****************************************************************************

  // Make a field group to communicate
  // *****************************************************************************
  FieldGroup com_group;
  // *****************************************************************************

  // Other objects
  // *****************************************************************************
  OwnBCs ownBC;        // Class containing methods which sets the ghost points
  OwnOperators *ownOp; // Class containing own differential operators
  OwnLaplacianInversions ownLapl; // Class containing own laplacian invertors
  // *****************************************************************************
  // ############################################################################
};

#endif
