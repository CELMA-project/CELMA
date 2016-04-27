// *************** Simulation of DDXDDXCylinderClean *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __DDXDDXCylinderClean_H__
#define __DDXDDXCylinderClean_H__

#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>              // Gives field factory
#include <bout/constants.hxx>             // Gives PI and TWOPI
#include <derivs.hxx>                     // Gives the derivatives
#include <difops.hxx>                     // Gives the diff options
#include <vecops.hxx>                     // Gives the vec diff options
#include "../../common/c/OwnBc.hxx"          // Gives method to set own BC
#include "../../common/c/SetInnerRho.hxx"    // Gives method to set inner rho

class DDXDDXCylinderClean : public PhysicsModel {
public:
    // Destructor
    ~DDXDDXCylinderClean();
protected:
    int init(bool restarting);
    int rhs(BoutReal t);
private:
    // Global variable initialization
    // ############################################################################
    // Variables
    // *****************************************************************************
    Field3D f;
    Field3D S, S_num;
    Field3D e;
    // *****************************************************************************

    // Auxiliary variables
    // *****************************************************************************
    Field3D DDXf;
    // *****************************************************************************

    // Constants
    // *****************************************************************************
    BoutReal Lx;    // The box dimensions
    // *****************************************************************************

    // Make a field group to communicate
    // *****************************************************************************
    FieldGroup com_group;
    // *****************************************************************************
    // ############################################################################
};

#endif
