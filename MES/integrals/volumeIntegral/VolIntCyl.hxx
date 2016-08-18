// *************** Simulation of VolIntCyl *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __VolIntCyl_H__
#define __VolIntCyl_H__

#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>              // Gives field factory
#include <bout/constants.hxx>             // Gives PI and TWOPI
// Gives own boundaries (doing so by setting ghost points)
#include "../../common/c/include/ownBCs.hxx"
// Give the integrators
#include "../../common/c/include/helperFunctions.hxx"

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
    BoutReal Lx, Ly;    // The box dimensions
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
    OwnBCs ownBC;           // Class containing methods which sets the ghost points
    // *****************************************************************************
    // ############################################################################
};

#endif
