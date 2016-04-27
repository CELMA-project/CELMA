// *************** Simulation of DDXDDXCylinder *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __DDXDDXCylinder_H__
#define __DDXDDXCylinder_H__

#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>              // Gives field factory
#include <bout/constants.hxx>             // Gives PI and TWOPI
#include <derivs.hxx>                     // Gives the derivatives
#include <difops.hxx>                     // Gives the diff options
#include <vecops.hxx>                     // Gives the vec diff options
#include "../../common/c/OwnBc.hxx"          // Gives method to set own BC
#include "../../common/c/SetInnerRho.hxx"    // Gives method to set inner rho

class DDXDDXCylinder : public PhysicsModel {
public:
    // Destructor
    ~DDXDDXCylinder();
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

    // Own operators
    // *****************************************************************************
    Field3D calcPartialRhoPhi(const Field3D &phi);
    // *****************************************************************************

    // Auxiliary functions
    // *****************************************************************************
    void DDX_yup_and_ydown(const Field3D &in, Field3D &out);
    void DDX_one_xz_plane_without_BC(Field3D const &in,
                                     int const &y_ind,
                                     Field3D &out);
    // *****************************************************************************
    // ############################################################################
};

#endif
