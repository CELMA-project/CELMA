// *************** Simulation of CelmaCurMom *********************
/* Geometry
 *  x - The radial coordinate (rho)      [nx and dx set from the grid file]
 *  y - The height of the cylinder (z)   [ny and dy set from the grid file]
 *  z - The azimuthal coordinate (theta) [nz and dz set from the grid file and
 *                                        internally in the BOUT++ framework]
 */

#ifndef __CELMARC_H__
#define __CELMARC_H__

#include <bout/physicsmodel.hxx>
#include <invert_laplace.hxx>     // Gives invert laplace option
#include <field_factory.hxx>      // Gives field factory
#include <derivs.hxx>             // Gives the bracket method
#include <difops.hxx>             // Gives the diff options
#include <vecops.hxx>             // Gives the vec diff options
#include <map>                    // Gives std::map
// Gives own boundaries (doing so by setting ghost points)
#include "../common/c/include/ownBCs.hxx"
// Gives own operators
#include "../common/c/include/ownOperators.hxx"
// Gives the noise generator
#include "../common/c/include/noiseGenerator.hxx"
// Gives own lowPass filter
#include "../common/c/include/ownFilters.hxx"
// Gives the monitors
#include "../common/c/include/ownMonitors.hxx"
// Gives the parameters
#include "../common/c/include/parameters.hxx"

class CelmaCurMom : public PhysicsModel
{
protected:
    int init(bool restarting);
    int rhs(BoutReal t);
public:
    // Constructor
    // ############################################################################
    CelmaCurMom();
    // ############################################################################

    // Variable initialization
    // ############################################################################
    // Evolved variables
    // *****************************************************************************
    Field3D jPar, momDensPar;  // Parallel velocities
    Field3D lnN;               // Logarithm of density
    Field3D vort;              // Vorticity
    // *****************************************************************************

    // Non-evolved stored fields
    // *****************************************************************************
    Field3D phi;      // Potential
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
    Field3D DDYGradPerpPhiGradPerpUI;
    Field3D vortAdv, vortParAdv;
    Field3D divParCur, divSourcePhi;
    Field3D vortParArtVisc, vortPerpArtVisc, vortHyperVisc;
    // *****************************************************************************

    // Artificial viscosities
    // *****************************************************************************
    // Parallel aritifical viscosities
    BoutReal artViscParLnN, artViscParJpar, artViscPerpJPar;
    BoutReal artViscParVort;
    // Perpendicular viscosities
    BoutReal artViscPerpLnN, artViscParMomDens, artViscPerpMomDens;
    BoutReal artViscPerpVort;
    // Azimuthal hyperviscosities
    BoutReal artHyperAzVort;
    // *****************************************************************************

    // Temporary fields
    // *****************************************************************************
    Vector3D gradPerpPhi, sourcePhi;
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
    Laplacian *phiSolver;           // Solver object for the FFT solver
    OwnFilters *ownFilter;          // Pointer to the chosen filter class
    OwnMonitors ownMon;             // Own monitors
    int outputMonitor(BoutReal simtime, int iter, int NOUT);    // Monitors every output
    // *****************************************************************************

    // Input switches
    // *****************************************************************************
    bool saveDdt;             // If ddt's should be saved
    bool saveTerms;           // If terms should be saved
    bool includeNoise;        // Include noise
    bool forceAddNoise;       // Add noise on restart as well
    bool useHyperViscAzVort;  // If hyperviscosity should be used in the vorticity
    bool monitorEnergy;       // If energy should be monitored
    bool monitorParticleNumber;            // If total particle number should be monitored
    bool constViscPar;        // If the input par viscosity is the simulation viscosity
    bool constViscPerp;       // If the input perp viscosity is the simulation viscosity
    bool constViscHyper;      // If the input hyper viscosity is the simulation viscosity
    bool viscosityGuard;      // If a check should be performed for the artificial viscosity
    // *****************************************************************************

    // Runtime switches
    // *****************************************************************************
    bool noiseAdded;          // A check whether the noise is added or not
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
