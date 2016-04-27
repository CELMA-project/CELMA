// *************** Simulation of NaulinSolver *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __NaulinSolver_H__
#define __NaulinSolver_H__

#include <bout/physicsmodel.hxx>
#include <invert_laplace.hxx>                   // Gives invert laplace option
#include <field_factory.hxx>                    // Gives field factory
#include <bout/constants.hxx>                   // Gives PI and TWOPI
#include <float.h>                              // Includes DBL_EPSILON
#include <math.h>                               // Includes fabs
#include <derivs.hxx>                           // Gives the bracket method
#include <difops.hxx>                           // Gives the diff options
#include <vecops.hxx>                           // Gives the vec diff options
#include "../../common/c/OwnBc.hxx"          // Gives method to set own BC
#include "../../common/c/SetInnerRho.hxx"    // Gives method to set inner rho

class NaulinSolver : public PhysicsModel {
public:
    // Destructor
    ~NaulinSolver();
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
    // *****************************************************************************

    // Swithces
    // *****************************************************************************
    bool NaulinMonitor;
    // *****************************************************************************

    // Additional methods and solvers
    // *****************************************************************************
    Laplacian* phiSolver;            // The laplacian solver
    // *****************************************************************************

    // Constants
    // *****************************************************************************
    BoutReal Lx;                     // The box dimensions
    // *****************************************************************************

    // Make a field group to communicate
    // *****************************************************************************
    FieldGroup com_group;
    // *****************************************************************************

    // Own laplacian solver
    // *****************************************************************************
    Field3D NSolver(const BoutReal &atol,
                    const BoutReal &rtol,
                    const int &maxit,
                    const Field3D &ln_n,
                    const Field3D &n,
                    const Field3D &vortD,
                    Field3D &phi,
                    Field3D &vort,
                    const bool &monitor);

    // Tolerances
    BoutReal atol;
    BoutReal rtol;
    int maxit;

    // Auxiliary variables
    Field3D phiCur;
    BoutReal E_abs_Linf;
    BoutReal E_rel_Linf;
    // *****************************************************************************

    // Own operators
    // *****************************************************************************
    const Vector3D own_Grad_perp(const Field3D &f);
    // *****************************************************************************
    // ############################################################################
};

#endif
