// *************** Simulation of cauchyBC *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __cauchyBC_H__
#define __cauchyBC_H__

#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>              // Gives field factory
#include <bout/constants.hxx>             // Gives PI and TWOPI
#include <derivs.hxx>                     // Gives the derivatives
#include <difops.hxx>                     // Gives the diff options
#include <vecops.hxx>                     // Gives the vec diff options
// Gives own boundaries (doing so by setting ghost points)
#include "../../common/c/include/ownBCs.hxx"
// Gives own operators
#include "../../common/c/include/ownOperators.hxx"

class UeSheath : public PhysicsModel {
public:
    // Destructor
    ~UeSheath();
protected:
    int init(bool restarting);
    int rhs(BoutReal t);
private:
    // Global variable initialization
    // ############################################################################
    // Variables
    // *****************************************************************************
    Field3D uEParOrigin, uEParWBC;
    Field3D phi;
    Field3D profile;
    Field3D e;
    // *****************************************************************************

    // Constants
    // *****************************************************************************
    BoutReal Lx, Ly;    // The box dimensions
    BoutReal Lambda;    // log(sqrt(mu/(2.0*pi)))
    BoutReal phiRef;    // Reference potential
    // *****************************************************************************

    // Make a field group to communicate
    // *****************************************************************************
    FieldGroup com_group;
    // *****************************************************************************

    // Other objects
    // *****************************************************************************
    OwnBCs ownBC;           // Class containing methods which sets the ghost points
    OwnOperators ownOp;     // Class containing own differential operators
    // *****************************************************************************
    // ############################################################################
};

#endif
