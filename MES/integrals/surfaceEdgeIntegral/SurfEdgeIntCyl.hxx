// *************** Simulation of SurfEdgeIntCyl *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __SurfEdgeIntCyl_H__
#define __SurfEdgeIntCyl_H__

#include <bout/physicsmodel.hxx>
#include <field_factory.hxx>              // Gives field factory
#include <bout/constants.hxx>             // Gives PI and TWOPI
#include <initialprofiles.hxx>            // Give initial_profiles
// Gives own boundaries (doing so by setting ghost points)
#include "../../common/c/include/ownBCs.hxx"
// Give the integrators
#include "../../common/c/include/helpers.hxx"

class SurfEdgeIntCyl : public PhysicsModel {
public:
    // Constructor (initialize the results vector using list initialization)
    SurfEdgeIntCyl() : results(4, 0.0) {}
    // Destructor
    ~SurfEdgeIntCyl();
protected:
    int init(bool restarting);
    int rhs(BoutReal t);
private:
    // Global variable initialization
    // ############################################################################
    // Variables
    // *****************************************************************************
    Vector3D v;
    BoutReal S_Xout, S_Yup, S_Ydown;
    BoutReal S_Xout_num, S_Yup_num, S_Ydown_num;
    BoutReal e_Xout, e_Yup, e_Ydown;
    std::vector<BoutReal> results;
    // *****************************************************************************

    // Constants
    // *****************************************************************************
    BoutReal Lx, Ly;    // The box dimensions
    // *****************************************************************************

    // Switches
    // *****************************************************************************
    bool saveFields;
    // *****************************************************************************

    // Monitor
    // *****************************************************************************
    int xIndInner;
    int xIndOuter;
    int yIndLower;
    int yIndUpper;
    int MXG;
    int MYG;
    // *****************************************************************************


    // Make a field group to communicate
    // *****************************************************************************
    FieldGroup com_group;
    // *****************************************************************************

    // Other objects
    // *****************************************************************************
    OwnBCs ownBC;            // Class containing methods which sets the ghost points
    SurfaceIntegral surfInt; // Surface integration class
    // *****************************************************************************
    // ############################################################################
};

#endif
