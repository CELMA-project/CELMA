// *************** Simulation of CelmaCurMom *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __CELMACURMOM_H__
#define __CELMACURMOM_H__

#include <bout/physicsmodel.hxx>
#include <invert_laplace.hxx>                   // Gives invert laplace option
#include <field_factory.hxx>                    // Gives field factory
#include <derivs.hxx>                           // Gives the bracket method
#include <difops.hxx>                           // Gives the diff options
#include <vecops.hxx>                           // Gives the vec diff options
// Gives own boundaries (doing so by setting ghost points)
#include "../common/c/include/ownBCs.hxx"
// Gives own operators
#include "../common/c/include/ownOperators.hxx"
// Gives own laplacian inversions
#include "../common/c/include/ownLaplacianInversions.hxx"
// Gives the noise generator
#include "../common/c/include/noiseGenerator.hxx"
// Gives own lowPass filter
#include "../common/c/include/ownFilters.hxx"

class CelmaCurMom : public PhysicsModel
{
protected:
    int init(bool restarting);
    int convective(BoutReal t);
    int diffusive(BoutReal t, bool linear);
private:
    // Global variable initialization
    // ############################################################################
    // Evolved variables
    // *****************************************************************************
    Field3D jPar, momDensPar;     // Parallel velocities
    Field3D lnN;               // Logarithm of density
    Field3D vortD;              // Vorticity like quantity
    // *****************************************************************************

    // Non-evolved stored fields
    // *****************************************************************************
    Field3D phi, vort;      // Potential and vorticity
    // *****************************************************************************

    // Fields making up the evolved variables which are going to be stored
    // *****************************************************************************
    // lnN fields
    Field3D lnNAdv, lnNRes, gradUEPar;
    Field3D lnNUeAdv, srcN, lnNParArtVisc, lnNPerpArtVisc;
    // jPar fields
    Field3D jParAdv, uIParAdvSum, uEParDoubleAdv;
    Field3D jParRes, elField;
    Field3D elPressure, neutralERes, neutralIRes;
    Field3D jParParArtVisc, jParPerpArtVisc;
    // momDensPar fields
    Field3D momDensAdv, densDiffusion;
    Field3D momDensParArtVisc, momDensPerpArtVisc;
    // Vorticity fields
    Field3D vortNeutral, potNeutral;
    Field3D parDerDivUIParNGradPerpPhi;
    Field3D vortDAdv, kinEnAdvN;
    Field3D divParCur, vortDParArtVisc, vortDPerpArtVisc;
    Field3D vortDhyperVisc;
    // *****************************************************************************

    // Temporary fields
    // *****************************************************************************
    Field3D DivUIParNGradPerpPhi;
    // *****************************************************************************

    // Auxiliary fields
    // *****************************************************************************
    Field3D n;               // Density
    Field3D uIPar, uEPar;    // Parallel currents
    Field3D dampingProfile;  // Radial profile
    Field3D S;               // Particle source
    Field3D invJ;            // 1/J (used in front of the bracket operator)
    Vector3D gradPerpLnN;    // gradPerpLnN
    // *****************************************************************************

    // Constants
    // *****************************************************************************
    BoutReal mu, Lambda;            // i-e mass ratio and ln(sqrt(mu/(2*PI)))
    BoutReal nuEI, nuEN, nuIN;      // electron, ion, and neutral collision freqs
    BoutReal a, bRho, bZ, cRho, cZ; // Variables used in the source term
    BoutReal Lx, Ly;                // The box dimensions
    // Parallel aritifical dissipation
    BoutReal artViscParLnN, artViscParJpar, artViscPerpJPar;
    BoutReal artViscParVortD;
    // Perpendicular dissipation
    BoutReal artViscPerpLnN, artViscParMomDens, artViscPerpMomDens;
    BoutReal artViscPerpVortD;
    // Azimuthal hyperviscosities
    BoutReal artHyperAzVortD;
    // *****************************************************************************

    // Additional methods and solvers
    // *****************************************************************************
    BRACKET_METHOD bm;              // The bracket method
    OwnBCs ownBC;                   // Class containing methods which sets the ghost points
    OwnOperators *ownOp;            // Pointer to the chosen operators class
    OwnLaplacianInversions ownLapl; // Class containing own laplacian
    OwnFilters *ownFilter;          // Pointer to the chosen filter class
    // *****************************************************************************

    // Switches
    // *****************************************************************************
    string  ownOpType;           // Type of own operators
    bool saveDdt;             // If ddt's should be saved
    bool saveTerms;           // If terms should be saved
    bool includeNoise;        // Include noise
    bool forceAddNoise;       // Add noise on restart as well
    bool noiseAdded;          // A check whether the noise is added or not
    bool useHyperViscAzVortD; // If hyperviscosity should be used in the vorticity
    // *****************************************************************************

    // Make a field group to communicate
    // *****************************************************************************
    FieldGroup comGroup;
    // *****************************************************************************
    // ############################################################################
};

#endif
