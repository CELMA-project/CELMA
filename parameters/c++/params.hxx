// *************** Simulation of Params *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __PARAMS_H__
#define __PARAMS_H__

#include <bout/physicsmodel.hxx>
#include "../../celma/common/c/include/parameters.hxx"

class Params : public PhysicsModel {
public:
    // Destructor
    ~Params();
protected:
    int init(bool restarting);
    int rhs(BoutReal t);
private:
    // Global variable initialization
    // ############################################################################
    // Constants
    // *****************************************************************************
    BoutReal radius, Lx; // Plasma radius
    BoutReal len, Ly;    // Cylinder length
    // *****************************************************************************

    // Input parameters
    // *****************************************************************************
    BoutReal n0;
    BoutReal Ti0;
    BoutReal Te0;
    BoutReal B0;
    BoutReal S;
    // *****************************************************************************

    // Output parameters
    // *****************************************************************************
    BoutReal nuEI;
    BoutReal mu;
    BoutReal nuS;
    BoutReal beta;
    BoutReal omCI;
    BoutReal rhoS;
    // *****************************************************************************
    // ############################################################################
};

#endif
