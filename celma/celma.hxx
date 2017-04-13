// ********************************* Celma ************************************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __CELMA_H__
#define __CELMA_H__

#include <bout/physicsmodel.hxx>
#include <invert_laplace.hxx>     // Gives invert laplace option
#include <field_factory.hxx>      // Gives field factory
#include <derivs.hxx>             // Gives the bracket method
#include <difops.hxx>             // Gives the diff options
#include <vecops.hxx>             // Gives the vec diff options
#include <map>                    // Gives std::map
#include <float.h>                // Includes DBL_EPSILON
// Gives own boundaries (doing so by setting ghost points)
#include "../common/c/include/ownBCs.hxx"
// Gives own operators
#include "../common/c/include/ownOperators.hxx"
// Gives own laplacian inversions
#include "../common/c/include/ownLaplacianInversions.hxx"
// Gives own lowPass filter
#include "../common/c/include/ownFilters.hxx"
// Gives the monitors
#include "../common/c/include/ownMonitors.hxx"
// Gives the parameters
#include "../common/c/include/parameters.hxx"

class Celma : public PhysicsModel
{
protected:
    int init(bool restarting);
    int rhs(BoutReal t);
public:
    // Constructor
    // ############################################################################
    Celma();
    // ############################################################################

    // Variable initialization
    // ############################################################################
    // Evolved variables
    // *****************************************************************************
    Field3D jPar, momDensPar;  // Parallel velocities
    Field3D lnN;               // Logarithm of density
    Field3D vortD;             // Vorticity like quantity
    // *****************************************************************************

    // Non-evolved stored fields
    // *****************************************************************************
    Field3D phi, vort;      // Potential and vorticity
    // *****************************************************************************

    // Auxiliary fields
    // *****************************************************************************
    Field3D n;               // Density
    Field3D uIPar, uEPar;    // Parallel currents
    Field3D S;               // Particle source
    Field3D invJ;            // 1/J (used in front of the bracket operator)
    Vector3D gradPerpLnN;    // gradPerpLnN
    // *****************************************************************************

    // Additional methods and solvers
    // *****************************************************************************
    OwnOperators *ownOp;            // Pointer to the chosen operators class
    // *****************************************************************************
    // ############################################################################
private:
    // Variable initialization
    // ############################################################################
    // Fields making up the evolved variables which are going to be stored
    // *****************************************************************************
    // lnN fields
    Field3D lnNAdv, lnNRes, gradUEPar;
    Field3D lnNUeAdv, srcN, lnNParArtVisc, lnNPerpArtVisc;
    // jPar fields
    Field3D jParAdv, uEParAdv, uIParAdv, jParParAdv;
    Field3D jParRes, gradPhiLnN;
    Field3D neutralERes, neutralIRes;
    Field3D jParParArtVisc, jParPerpArtVisc;
    // momDensPar fields
    Field3D momDensAdv, uIFluxAdv, elPressure, densDiffusion;
    Field3D neutralEResMu, momDensParArtVisc, momDensPerpArtVisc;
    // Vorticity fields
    Field3D vortNeutral, potNeutral;
    Field3D parDerDivUIParNGradPerpPhi;
    Field3D vortDAdv, kinEnAdvN;
    Field3D divParCur, vortDParArtVisc, vortDPerpArtVisc;
    Field3D vortDHyperVisc;
    // *****************************************************************************

    // Artificial viscosities
    // *****************************************************************************
    // Parallel aritifical viscosities
    BoutReal artViscParLnN, artViscParJpar, artViscPerpJPar;
    BoutReal artViscParVortD;
    // Perpendicular viscosities
    BoutReal artViscPerpLnN, artViscParMomDens, artViscPerpMomDens;
    BoutReal artViscPerpVortD;
    // Azimuthal hyperviscosities
    BoutReal artHyperAzVortD;
    // *****************************************************************************

    // Temporary fields
    // *****************************************************************************
    Field3D divUIParNGradPerpPhi;
    Vector3D gradPerpPhi;
    // *****************************************************************************

    // Geometry parameters
    // *****************************************************************************
    BoutReal Lx, Ly; // The box dimensions
    // *****************************************************************************

    // Input parameters
    // *****************************************************************************
    BoutReal    radius; // Plasma radius
    BoutReal    length; // Cylinder length
    BoutReal    n0;
    BoutReal    Ti0,Te0;
    BoutReal    B0;
    BoutReal    Sn;    // Particle source amplitude
    BoutReal    nn;
    std::string gas;
    bool        warningForException;
    // *****************************************************************************

    // Output parameters
    // *****************************************************************************
    BoutReal LxParam;          // Normalized plasma radius
    BoutReal LyParam;          // Normalized cylinder length
    BoutReal nuEI, nuEN, nuIN; // Electron, ion, and neutral collision freqs
    BoutReal mu, Lambda;       // Mass ratio and Lambda
    BoutReal SNorm;
    BoutReal beta;
    BoutReal omCI, rhoS;
    // *****************************************************************************

    // Helping variables
    // *****************************************************************************
    BoutReal eta0INorm; // Viscous coefficient
    // *****************************************************************************

    // Additional methods and solvers
    // *****************************************************************************
    string ownOpType;               // Type of own operators
    BRACKET_METHOD bm;              // The bracket method
    OwnBCs ownBC;                   // Class containing methods which sets the ghost points
    OwnLaplacianInversions ownLapl; // Class containing own laplacian
    OwnFilters *ownFilter;          // Pointer to the filter class
    OwnMonitors ownMon;             // Own monitors
    int outputMonitor(BoutReal simtime, int iter, int NOUT);    // Monitors every output
    // *****************************************************************************

    // Input switches
    // *****************************************************************************
    bool saveDdt;             // If ddt's should be saved
    bool saveTerms;           // If terms should be saved
    bool forceAddNoise;       // Add noise on restart as well
    bool useHyperViscAzVortD; // If hyperviscosity should be used in the vorticity
    bool monitorEnergy;       // If energy should be monitored
    bool monitorParticleNumber;            // If total particle number should be monitored
    bool constViscPar;        // If the input par viscosity is the simulation viscosity
    bool constViscPerp;       // If the input perp viscosity is the simulation viscosity
    bool constViscHyper;      // If the input hyper viscosity is the simulation viscosity
    bool viscosityGuard;      // If a check should be performed for the artificial viscosity
    // *****************************************************************************

    // Make a field group to communicate
    // *****************************************************************************
    FieldGroup comGroup;
    // *****************************************************************************

    // Monitors
    // *****************************************************************************
    std::map<std::string, BoutReal> kinE;
    std::map<std::string, BoutReal> potE;
    std::map<std::string, BoutReal> particleNumber;
    // *****************************************************************************

    // Initialization functions
    // *****************************************************************************
    void initializeOwnObjects();
    void setAndSaveParameters();
    void printPointsPerRhoS();
    void setAndSaveSource();
    void setSwithces(bool &restarting);
    void setAndSaveViscosities();
    // *****************************************************************************

    // Timestep initialization
    // *****************************************************************************
    void timestepInitialization();
    // *****************************************************************************
    // ############################################################################
};

#endif
