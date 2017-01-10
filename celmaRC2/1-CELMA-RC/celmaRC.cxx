// *************** Simulation of CELMACURMOM *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "celmaRC.hxx"

// Initialization and solving of the physics
// ############################################################################
int CelmaCurMom::init(bool restarting)
{
    TRACE("Halt in CelmaCurMom::init");

    // Initialize non-standard BOUT++ objects
    initializeOwnObjects();
    // Set the input
    setAndSaveParameters();
    // Print points per rhoS
    printPointsPerRhoS();
    // Set the source
    setAndSaveSource();
    // Set the switches
    setSwithces(restarting);
    // Set and save the viscosities
    setAndSaveViscosities();

    // Additional fields
    // ************************************************************************
    // The initial potential (obtained from the input)
    phi = FieldFactory::get()
          ->create3D("phi:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // The metric coefficient (needed in front of the arakawa bracket)
    invJ = (1.0/mesh->J);
    // ************************************************************************

    // Specifying the brackets to the arakawa scheme
    // (see examples/MMS/hw for a elegant way to choose from the input file)
    // ************************************************************************
    bm = BRACKET_ARAKAWA;
    // ************************************************************************

    // Add a FieldGroup to communicate
    // ************************************************************************
    // NOTE: We only communicate variables we are taking derivatives of
    comGroup.add(lnN);
    comGroup.add(n);
    comGroup.add(momDensPar);
    comGroup.add(jPar);
    comGroup.add(uEPar);
    comGroup.add(uIPar);
    comGroup.add(phi);
    comGroup.add(vortD);
    // ************************************************************************

    // Specify BC for n and uIPar (used to set BC for jPar and momDensPar)
    // ************************************************************************
    n    .setBoundary("n"    );
    uIPar.setBoundary("uIPar");
    // ************************************************************************

    // Specify what values should be stored in the .dmp file
    // ************************************************************************
    // Variables to be saved repeatedly
    // vort and phi
    SAVE_REPEAT2(vort, phi);
    if(saveTerms){
        // lnN terms
        SAVE_REPEAT3(lnNAdv, lnNRes, gradUEPar);
        SAVE_REPEAT4(lnNUeAdv, srcN, lnNParArtVisc, lnNPerpArtVisc);
        // jPar terms
        SAVE_REPEAT4(jParAdv, uEParAdv, uIParAdv, jParParAdv);
        SAVE_REPEAT4(jParRes, gradPhiLnN, neutralERes, neutralIRes);
        SAVE_REPEAT2(jParParArtVisc, jParPerpArtVisc);
        // momDensPar terms
        SAVE_REPEAT4(momDensAdv, uIFluxAdv, elPressure, densDiffusion);
        SAVE_REPEAT3(neutralEResMu, momDensParArtVisc, momDensPerpArtVisc);
        // Vorticity terms
        SAVE_REPEAT2(vortNeutral, potNeutral);
        SAVE_REPEAT3(divParCur, vortDParArtVisc, vortDPerpArtVisc);
        SAVE_REPEAT (parDerDivUIParNGradPerpPhi);
        SAVE_REPEAT2(vortDAdv, kinEnAdvN);
        // Helping fields
        SAVE_REPEAT2(uIPar, uEPar);

        if(saveDdt){
            SAVE_REPEAT4(ddt(vortD), ddt(lnN), ddt(momDensPar), ddt(jPar));
        }
    }
    // Monitor variables to be solved for
    if(monitorEnergy){
        for (std::map<std::string, BoutReal>::iterator it=kinE.begin();
             it!=kinE.end();
             ++it){
            dump.add(it->second, it->first.c_str(), 1);
        }
        for (std::map<std::string, BoutReal>::iterator it=potE.begin();
             it!=potE.end();
             ++it){
            dump.add(it->second, it->first.c_str(), 1);
        }
    }
    if(monitorParticleNumber){
        for (std::map<std::string,BoutReal>::iterator it=particleNumber.begin();
             it!=particleNumber.end();
             ++it){
            dump.add(it->second, it->first.c_str(), 1);
        }
    }
    // Variables to be solved for
    SOLVE_FOR4(vortD, lnN, momDensPar, jPar);
    //*************************************************************************

    return 0;
}

int CelmaCurMom::rhs(BoutReal t)
{
    TRACE("Halt in CelmaCurMom::rhs");

    timestepInitialization();

    /* NOTE: Bracket operator
     * The bracket operator returns -(grad(phi) x b/B)*grad(f)
     * As B is defined by in field aligned coordinates, and we are working with
     * cylindrical coordinates (which happens to coincide with the metrics of
     * that of the field aligned coordinate system) we need to multiply with
     * 1/J
     */

    // Terms in lnNPar
    // ************************************************************************
    lnNAdv         = - invJ*bracket(phi, lnN, bm);
    gradUEPar      = - DDY(uEPar);
    lnNUeAdv       = - Vpar_Grad_par(uEPar, lnN);
    srcN           =   S/n;
    lnNRes         = (nuEI/mu)*(Laplace_perp(lnN)+gradPerpLnN*gradPerpLnN);
    lnNParArtVisc  = artViscParLnN*D2DY2(lnN);
    lnNPerpArtVisc = artViscPerpLnN*Laplace_perp(lnN);
    ddt(lnN) =
          lnNAdv
        + gradUEPar
        + lnNUeAdv
        + srcN
        + lnNRes
        + lnNParArtVisc
        + lnNPerpArtVisc
        ;
    // Filtering highest modes
    ddt(lnN) = ownFilter->ownFilter(ddt(lnN));
    // ************************************************************************


    // Terms in jPar
    // ************************************************************************
    jParAdv         = - invJ*bracket(phi, jPar, bm);
    uEParAdv        =   n*(Vpar_Grad_par(uEPar, uEPar));
    uIParAdv        = - n*(Vpar_Grad_par(uIPar, uIPar));
    jParParAdv      = - Vpar_Grad_par((jPar/n), uEPar);
    gradPhiLnN      = mu*n*DDY(lnN - phi);
    jParRes         = - 0.51*nuEI*jPar;
    neutralERes     = n*nuEN*uEPar;
    neutralIRes     = - n*nuIN*uIPar;
    jParParArtVisc  = (artViscParJpar)*D2DY2(jPar);
    jParPerpArtVisc = (artViscPerpJPar)*Laplace_perp(jPar);

    ddt(jPar) =
          jParAdv
        + uEParAdv
        + uIParAdv
        + jParParAdv
        + gradPhiLnN
        + jParRes
        + neutralERes
        + neutralIRes
        + jParParArtVisc
        + jParPerpArtVisc
        ;
    // Filtering highest modes
    ddt(jPar) = ownFilter->ownFilter(ddt(jPar));
    // ************************************************************************


    // Terms in momDensPar
    // ************************************************************************
    momDensAdv = - invJ*bracket(phi, momDensPar, bm);
 // UIParAdv calculated above
    uIFluxAdv     = - Vpar_Grad_par(uIPar, n*uEPar);
    elPressure    = - DDY(n);
    densDiffusion = nuEI*(uIPar/mu)*Laplace_perp(n);
 // neutralIRes calculated above
    neutralEResMu      = neutralERes/mu;
    momDensParArtVisc  = (artViscParMomDens)*D2DY2(momDensPar);
    momDensPerpArtVisc = (artViscPerpMomDens)*Laplace_perp(momDensPar);

    ddt(momDensPar) =
        + momDensAdv
        + uIFluxAdv
        + uIParAdv
        + elPressure
        + densDiffusion
        + neutralIRes
        + neutralEResMu
        + momDensParArtVisc
        + momDensPerpArtVisc
        ;
    // Filtering highest modes
    ddt(momDensPar) = ownFilter->ownFilter(ddt(momDensPar));
    // ************************************************************************


    // Preparation
    // ************************************************************************
    divUIParNGradPerpPhi = ownOp->div_f_GradPerp_g(uIPar*n, phi);
    // Set the ghost points in order to take DDY
    ownBC.extrapolateYGhost(divUIParNGradPerpPhi);
    // We must communicate as we will take DDY
    mesh->communicate(divUIParNGradPerpPhi);
    // Saving gradPerpPhi for use in monitors
    gradPerpPhi = ownOp->Grad_perp(phi);
    // ************************************************************************


    // Terms in vorticity
    // ************************************************************************
    vortNeutral = - nuIN*n*vort;
    potNeutral  = - nuIN*gradPerpPhi*ownOp->Grad_perp(n);
    vortDAdv    = - ownOp->vortDAdv(phi, vortD);
    kinEnAdvN   = - ownOp->kinEnAdvN(phi, n);
    parDerDivUIParNGradPerpPhi = - DDY(divUIParNGradPerpPhi);
    divParCur                  =   DDY(jPar);
    vortDParArtVisc  =   artViscParVortD*D2DY2(vortD);
    vortDPerpArtVisc =   artViscPerpVortD*Laplace_perp(vortD);
    vortDHyperVisc   = - artHyperAzVortD*D4DZ4(vortD);

    ddt(vortD) =
          vortNeutral
        + potNeutral
        + parDerDivUIParNGradPerpPhi
        + divParCur
        + vortDAdv
        + kinEnAdvN
        + vortDParArtVisc
        + vortDPerpArtVisc
        + vortDHyperVisc
        ;

    // Filtering highest modes
    ddt(vortD) = ownFilter->ownFilter(ddt(vortD));
    // ************************************************************************
    return 0;
}
// ############################################################################

// Constructor
// ############################################################################
CelmaCurMom::CelmaCurMom()
/* FIXME: c++11 is unsupported on jess
 * :
 *     kinE ({{"perpKinEE", 0.0}, {"parKinEE", 0.0}, {"sumKinEE", 0.0},
 *            {"perpKinEI", 0.0}, {"parKinEI", 0.0}, {"sumKinEI", 0.0},
 *            {"polAvgPerpKinEE", 0.0}, {"polAvgParKinEE", 0.0}, {"polAvgSumKinEE", 0.0},
 *            {"polAvgPerpKinEI", 0.0}, {"polAvgParKinEI", 0.0}, {"polAvgSumKinEI", 0.0}}
 *            ),
 *     potE ({{"potEE", 0.0}, {"polAvgPotEE", 0.0}}),
 *     particleNumber ({{"particleNumber", 0.0}})
 */
{
    TRACE("Halt in CelmaCurMom::CelmaCurMom");

    // Non c++11 initialization
    kinE["perpKinEE"]                = 0.0;
    kinE["parKinEE"]                 = 0.0;
    kinE["sumKinEE"]                 = 0.0;
    kinE["perpKinEI"]                = 0.0;
    kinE["parKinEI"]                 = 0.0;
    kinE["sumKinEI"]                 = 0.0;
    kinE["polAvgPerpKinEE"]          = 0.0;
    kinE["polAvgParKinEE"]           = 0.0;
    kinE["polAvgSumKinEE"]           = 0.0;
    kinE["polAvgPerpKinEI"]          = 0.0;
    kinE["polAvgParKinEI"]           = 0.0;
    kinE["polAvgSumKinEI"]           = 0.0;
    potE["potEE"]                    = 0.0;
    potE["polAvgPotEE"]              = 0.0;
    particleNumber["particleNumber"] = 0.0;
}
// ############################################################################

// Monitor
// ############################################################################
int CelmaCurMom::outputMonitor(BoutReal simtime, int iter, int NOUT)
{
    TRACE("Halt in CelmaCurMom::outputMonitor");

    if(monitorEnergy || monitorParticleNumber){
        ownMon.calcPolAvgN(n);
    }
    if(monitorEnergy){
        ownMon.kinEnergy(n, gradPerpPhi, uEPar, uIPar, &kinE);
        ownMon.potEnergy(n, &potE);
    }
    if(monitorParticleNumber){
        ownMon.numberOfParticles(n, &particleNumber);
    }

    return 0;
}
// ############################################################################

// Initialization helpers
// ############################################################################
void CelmaCurMom::initializeOwnObjects()
{
    TRACE("CelmaCurMom::initializeOwnObjects");

    // Create OwnOperators
    // ************************************************************************
    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    /* NOTE: Calls createOperators without making an object of OwnOperators.
     *       The child is typecasted to the parent
     */
    ownOp = OwnOperators::createOperators();
    ownLapl.create(ownOp, ownBC);
    // As a small hack we find the own operator type to figure out how to
    // treat the div(ue*grad(n grad_perp phi))
    Options *ownOp = options->getSection("ownOperators");
    ownOp->get("type", ownOpType, "BasicBrackets");
    if (
        (lowercase(ownOpType) == lowercase("simpleStupid")) ||
        (lowercase(ownOpType) == lowercase("onlyBracket"))
       ){
        throw BoutException(
                "simpleStupid and onlyBracket not available in this model"
                );
    }
    // ************************************************************************

    // Create the filter
    // ************************************************************************
    /* NOTE: Calls createFilter without making an object of ownFilter.
     *       The child is typecasted to the parent
     */
    ownFilter = OwnFilters::createFilter();
    // ************************************************************************
}

void CelmaCurMom::setAndSaveParameters()
{
    TRACE("Halt in CelmaCurMom::setAndSaveParameters");

    // Load from the geometry
    // ************************************************************************
    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Get the section of the variables from [geom] specified in BOUT.inp
    // or in the command-line arguments
    Options *geom = options->getSection("geom");
    geom->get("Lx", Lx, 0.0);
    geom->get("Ly", Ly, 0.0);
    // ************************************************************************

    // Load the input and print the output
    // ************************************************************************
    Options *input = options->getSection("input");
    input->get("radius", radius, 0.0);
    input->get("length", length, 0.0);
    input->get("n0"    , n0    , 0.0);
    input->get("Te0"   , Te0   , 0.0);
    input->get("Ti0"   , Ti0   , 0.0);
    input->get("B0"    , B0    , 0.0);
    input->get("Sn"    , Sn    , 0.0);
    input->get("nn"    , nn    , 0.0);
    input->get("gas"   , gas   , "H");
    // Cast to upper
    gas[0] = toupper(gas[0]);

    input->get("warningForException",
                warningForException , false);

    Parameters params(radius, length, n0, Te0, Ti0, B0, Sn, nn, gas,
                      warningForException);
    // ************************************************************************

    // Get the variables
    // ************************************************************************
    LxParam   = params.getLx();
    LyParam   = params.getLy();
    nuEI      = params.getNuEINorm();
    nuEN      = params.getNuENNorm();
    nuIN      = params.getNuINNorm();
    SNorm     = params.getSNorm();
    mu        = params.getMu();
    Lambda    = params.getLambda();
    beta      = params.getBeta();
    omCI      = params.getOmCI();
    rhoS      = params.getRhoS();
    eta0INorm = params.getEta0INorm();
    // ************************************************************************

    // Check that Lx and LxParams is the same up until the fourt decimal point
    // ************************************************************************
    std::ostringstream stream;
    int precision = 4;
    bool throwError = false;
    if ( ( fabs(round(LxParam*1e4)/1e4 - round(Lx*1e4)/1e4))>DBL_EPSILON ){
        stream << "Mismatch between 'Lx' calculated from 'radius' "
               << "and input 'Lx'\n"
               << "Calculated = "
               << std::fixed
               << std::setprecision(precision) << round(LxParam*1e4)/1e4
               << "\n"
               << "Input      = "
               << std::fixed
               << std::setprecision(precision) << round(Lx*1e4)/1e4
               << "\n\n";
        throwError = true;
    }
    if ( ( fabs(round(LyParam*1e4)/1e4 - round(Ly*1e4)/1e4))>DBL_EPSILON ){
        stream << "Mismatch between 'Ly' calculated from 'length' "
               << "and input 'Ly'\n"
               << "Calculated = "
               << std::fixed
               << std::setprecision(precision) << round(LyParam*1e4)/1e4
               << "\n"
               << "Input      = "
               << std::fixed
               << std::setprecision(precision) << round(Ly*1e4)/1e4
               << "\n\n";
        throwError = true;
    }
    if(throwError){
        std::string str = stream.str();
        // Cast the stream to a const char in order to use it in BoutException
        const char* message = str.c_str();
        throw BoutException(message);
    }
    // ************************************************************************

    // Save the input
    // ************************************************************************
    SAVE_ONCE5(n0, Te0, Ti0, B0, Sn);
    // Save the variables
    SAVE_ONCE2(Lx, Ly);
    SAVE_ONCE4(nuEI, nuEN, nuIN, SNorm);
    SAVE_ONCE2(mu, Lambda);
    SAVE_ONCE (beta);
    SAVE_ONCE2(omCI, rhoS);
    SAVE_ONCE (eta0INorm);
    // ************************************************************************
}

void CelmaCurMom::printPointsPerRhoS()
{
    TRACE("Halt in CelmaCurMom::printPointsPerRhoS");

    // Get the option (before any sections) in the BOUT.inp file
    int MXG;
    BoutReal pointsPerRhoSRadially;
    BoutReal pointsPerRhoSParallely;
    BoutReal pointsPerRhoSAzimuthally;
    BoutReal minPointsPerRhoSXZ;
    BoutReal minPointsPerRhoSY;

    Options *root = Options::getRoot();
    root->get("MXG", MXG, 0);

    if (MXG == 0){
        throw BoutException("No MXG found, please specify.");
    }

    // dx = Lx/(nx-2*MXG) => nx = (Lx/dx) + 2*MXG
    pointsPerRhoSRadially = ((Lx/mesh->dx(0, 0)) + 2*MXG)/Lx;
    // dy = Ly/ny => ny = Ly/dy => ny/Ly = 1/dy
    pointsPerRhoSParallely = 1.0/mesh->dy(0 ,0);
    // O=2*pi*r, so on edge nz/rho_s = nz/(2*pi*Lx)
    pointsPerRhoSAzimuthally = (mesh->ngz - 1)/(2.0*PI*Lx);

    root->getSection("geom")->get("minPointsPerRhoSXZ",minPointsPerRhoSXZ,3.0);
    root->getSection("geom")->get("minPointsPerRhoSY",minPointsPerRhoSY,1.0e-1);

    std::ostringstream stream;
    bool throwError = false;
    if (pointsPerRhoSRadially < minPointsPerRhoSXZ){
        stream << "Minimum points per rhoS not fulfilled in x.\n"
               << "Limit is         " << minPointsPerRhoSXZ << "\n"
               << "Current value is " << pointsPerRhoSRadially << "\n\n";
        throwError = true;
    }
    if (mesh->ngz > 3 && pointsPerRhoSAzimuthally < minPointsPerRhoSXZ){
        stream << "Minimum points per rhoS not fulfilled on outer circumference.\n"
               << "Limit is         " << minPointsPerRhoSXZ << "\n"
               << "Current value is " << pointsPerRhoSAzimuthally << "\n\n";
        throwError = true;
    }
    if (pointsPerRhoSParallely < minPointsPerRhoSY){
        stream << "Minimum points per rhoS not fulfilled in y.\n"
               << "Limit is         " << minPointsPerRhoSY << "\n"
               << "Current value is " << pointsPerRhoSParallely << "\n\n";
        throwError = true;
    }
    if(throwError){
        std::string str =  stream.str();
        // Cast the stream to a const char in order to use it in BoutException
        const char* message = str.c_str();
        throw BoutException(message);
    }

    output << '\n' << std::string(51, '*') << std::endl;
    output << "Points per rhoS in x                   = "
           << pointsPerRhoSRadially << std::endl;
    output << "Points per rhoS on outer circumference = "
           << pointsPerRhoSAzimuthally << std::endl;
    output << "Points per rhoS in y                   = "
           << pointsPerRhoSParallely << std::endl;
    output << std::string(51, '*') << '\n' << std::endl;
}

void CelmaCurMom::setAndSaveSource()
{
    TRACE("Halt in CelmaCurMom::setAndSaveSource");

    BoutReal radialWidth, radialCentre, radialSteepness;
    BoutReal parallelWidth, parallelCentre, parallelSteepness;

    // Get the source constants, save them and create the source
    // ************************************************************************
    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    Options *thesource = options->getSection("thesource");
    thesource->get("radialWidth"      , radialWidth      , 0.0);
    thesource->get("radialCentre"     , radialCentre     , 0.0);
    thesource->get("radialSteepness"  , radialSteepness  , 0.0);
    thesource->get("parallelWidth"    , parallelWidth    , 0.0);
    thesource->get("parallelCentre"   , parallelCentre   , 0.0);
    thesource->get("parallelSteepness", parallelSteepness, 0.0);

    SAVE_ONCE3(radialWidth  , radialCentre  , radialSteepness  );
    SAVE_ONCE3(parallelWidth, parallelCentre, parallelSteepness);

    // The source (obtained from the input, and multiplied with SNorm)
    S = SNorm*
        FieldFactory::get()
      ->create3D("theSource:S", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // Save the source field
    SAVE_ONCE(S);
    // ************************************************************************
}

void CelmaCurMom::setSwithces(bool &restarting)
{
    TRACE("Halt in CelmaCurMom::setSwithces");

    // Get the switches
    // ************************************************************************
    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    Options *switches = options->getSection("switch");
    switches->get("useHyperViscAzVortD"   , useHyperViscAzVortD   , false);
    switches->get("includeNoise"          , includeNoise          , false);
    switches->get("forceAddNoise"         , forceAddNoise         , false);
    switches->get("saveDdt"               , saveDdt               , false);
    switches->get("constViscPar"          , constViscPar          , false);
    switches->get("constViscPerp"         , constViscPerp         , false);
    switches->get("constViscHyper"        , constViscHyper        , false);
    switches->get("saveTerms"             , saveTerms             , true );
    switches->get("monitorEnergy"         , monitorEnergy         , true );
    switches->get("monitorParticleNumber" , monitorParticleNumber , true );
    switches->get("viscosityGuard"        , viscosityGuard        , true );
    noiseAdded = false;
    // Decide whether noise should be added upon restart
    if (restarting && includeNoise && !(forceAddNoise)){
        output << "\n\n!!!!Warning!!!\n"
               << "restarting = true, includeNoise = true, forceAddNoise = false\n"
               << "Since forceAddNoise = false => program reset includeNoise to false"
               << "\n\n"
               << std::endl;
        includeNoise = false;
        noiseAdded = true; // For extra safety measurements
    }
    // ************************************************************************
}

void CelmaCurMom::setAndSaveViscosities()
{
    TRACE("Halt in CelmaCurMom::setAndSaveViscosities");

    // Get and save the viscosities
    // ************************************************************************
    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    Options *visc = options->getSection("visc");
    visc->get("artViscParLnN"     , artViscParLnN     , 0.0);
    visc->get("artViscParJpar"    , artViscParJpar    , 0.0);
    visc->get("artViscParMomDens" , artViscParMomDens , 0.0);
    visc->get("artViscParVortD"   , artViscParVortD   , 0.0);
    visc->get("artViscPerpLnN"    , artViscPerpLnN    , 0.0);
    visc->get("artViscPerpJPar"   , artViscPerpJPar   , 0.0);
    visc->get("artViscPerpMomDens", artViscPerpMomDens, 0.0);
    visc->get("artViscPerpVortD"  , artViscPerpVortD  , 0.0);
    visc->get("artHyperAzVortD"   , artHyperAzVortD   , 0.0);

    // Calculate diffusion from grid size
    if (!constViscPar){
        // SQ is squaring the expression
        // dx and dy are Field2D (0th index is ghost, but gives no problems)
        artViscParLnN     *= SQ(mesh->dy(0,0));
        artViscParJpar    *= SQ(mesh->dy(0,0));
        artViscParMomDens *= SQ(mesh->dy(0,0));
        artViscParVortD   *= SQ(mesh->dy(0,0));
    }

    if (!constViscPerp){
        // The perpednicular diffusion is in our model
        /* NOTE: Chosen independent of dz
         *       This makes artVisc constant when expanding restarts
         */
        artViscPerpLnN     *= SQ(mesh->dx(0,0));
        artViscPerpJPar    *= SQ(mesh->dx(0,0));
        artViscPerpMomDens *= SQ(mesh->dx(0,0));
        artViscPerpVortD   *= SQ(mesh->dx(0,0));
    }

    // Set artificial viscosities to 0 if useHyperViscAzVortD is false
    if (!useHyperViscAzVortD){
        output << "Setting artHyperAzVortD = 0.0 as useHyperViscAzVortD = False"
               << std::endl;
        artHyperAzVortD = 0.0;
    }
    if (!constViscHyper){
        // Azimuthal hyperviscosities
        artHyperAzVortD *= SQ(SQ(mesh->dz));
    }

    // Print and store the variables
    output << "\nPrinting real artificial viscosity coefficients" << std::endl;
    output << "***********************************************"   << std::endl;
    output << "Perpendicular";
    if (!constViscPerp){
        output << " (SQ(mesh->dx(0,0)) = " << SQ(mesh->dx(0,0)) << "):";
    }
    output << std::endl;
    output << "    For ln(n)    : "  << artViscPerpLnN     << std::endl;
    output << "    For j_{\\|}   : " << artViscPerpJPar    << std::endl;
    output << "    For nu_{i,\\|}: " << artViscPerpMomDens << std::endl;
    output << "    For vortD    : "  << artViscPerpVortD   << std::endl;
    output << "Parallel";
    if (!constViscPar){
        output << " (SQ(mesh->dy(0,0)) = "
                  << SQ(mesh->dy(0,0)) << "):";
    }
    output << std::endl;
    output << "    For ln(n)    : "  << artViscParLnN     << std::endl;
    output << "    For j_{\\|}   : " << artViscParJpar    << std::endl;
    output << "    For nu_{i,\\|}: " << artViscParMomDens << std::endl;
    output << "    For vortD    : "  << artViscParVortD   << std::endl;
    output << "Azimuthal hyperviscosity";
    if (!constViscHyper){
        output << "Azimuthal hyperviscosity (SQ(SQ(mesh->dz)) = "
                  << SQ(SQ(mesh->dz)) << "):";
    }
    output << std::endl;
    output << "    For vortD   : " << artHyperAzVortD << std::endl;
    output << "***********************************************\n" << std::endl;

    // Guards
    std::vector<BoutReal> perpVisc (4, 0.0);
    perpVisc.push_back(artViscPerpLnN    );
    perpVisc.push_back(artViscPerpJPar   );
    perpVisc.push_back(artViscPerpMomDens);
    perpVisc.push_back(artViscPerpVortD  );

    for (std::vector<BoutReal>::iterator it = perpVisc.begin();
         it != perpVisc.end();
         ++it){
        if(*it > nuEI){
            throw BoutException("One of the perpendicular viscosities "
                                "is larger than nuEI");
        }
    }

    std::vector<BoutReal> parVisc (4, 0.0);
    perpVisc.push_back(artViscParLnN    );
    perpVisc.push_back(artViscParJpar   );
    perpVisc.push_back(artViscParMomDens);
    perpVisc.push_back(artViscParVortD  );

    if(viscosityGuard){
        for (std::vector<BoutReal>::iterator it = perpVisc.begin();
             it != perpVisc.end();
             ++it){
            if(*it > 100.0*eta0INorm){
                throw BoutException("One of the parallel viscosities "
                                    "is 100.0 times larger than eta0I");
            }
        }
    }

    SAVE_ONCE2(artViscParLnN, artViscParJpar);
    SAVE_ONCE2(artViscParMomDens, artViscParVortD);
    SAVE_ONCE2(artViscPerpLnN, artViscPerpJPar);
    SAVE_ONCE2(artViscPerpMomDens, artViscPerpVortD);
    SAVE_ONCE (artHyperAzVortD);
    // ************************************************************************
}
// ############################################################################

// Timestep initialization
// ############################################################################
void CelmaCurMom::timestepInitialization()
{
    TRACE("Halt in CelmaCurMom::timestepInitialization");

    if (includeNoise && !noiseAdded){
        /* NOTE: Positioning of includeNoise
         * Field is reloaded from restart files if noise is added in init
         */
        // Class containing the noise generators
        // Calls the constructor with default arguments
        NoiseGenerator noise("geom", 3);
        // Add noise
        //*********************************************************************
        noise.generateRandomPhases(lnN, 1.0e-3);
        //*********************************************************************

        // Declare that the noise has been added
        noiseAdded = true;
    }

    // Manually specifying rho inner ghost points
    // ************************************************************************
    ownBC.innerRhoCylinder(lnN);
    ownBC.innerRhoCylinder(jPar);
    ownBC.innerRhoCylinder(momDensPar);
    ownBC.innerRhoCylinder(phi);    // Used to set BC in lapalace inversion
    ownBC.innerRhoCylinder(vortD);  // Later taken derivative of
    // The inner boundaries of phi is set in the inversion procedure
    // ************************************************************************

    // Preparations
    // ************************************************************************
    /* NOTE: Preparations
     * 1. gradPerp(lnN) is input to the NaulinSolver , thus
     *    a) gradPerpLnN needs to be calculated...
     *    b) ...which requires that lnN has been communicated
     * 2. n is input to the NaulinSolver, so
     *    a) It needs to be calculated
     *    b) We will communicate as the derivative is needed in vortD
     *
     * Note also:
     * 1. No further derivatives is taken for GradPerp(lnN), so this is not
     *    needed to be communicated.
     * 2. We are taking the derivative of phi in the NaulinSolver, however,
     *    this happens in a loop, and it is communicated before taking the
     *    derivative
     */
    mesh->communicate(lnN);
    gradPerpLnN = ownOp->Grad_perp(lnN);
    n = exp(lnN);
    uIPar = momDensPar/n;           // momDensPar = n*uIPar
    uEPar = - (jPar)/n + uIPar;     // j = n(uIPar - uEPar)
    // ************************************************************************

    // Laplace inversion
    // ************************************************************************
    /* NOTE: Solve for phi
     * 1. phi and vort is output
     * 2. The solver takes care of perpendicular boundaries
     * 3. Filtering of higher modes is done internally with the filter flag in
     *    BOUT.inp
     */
    phi = ownLapl.NaulinSolver(gradPerpLnN, n, vortD, phi, vort);
    // Filter
    phi = ownFilter->ownFilter(phi);
    // ************************************************************************

    // Treating parallel boundary conditions
    // ************************************************************************
    /* NOTE: No boundary condition on phi
     * - We need to know the parallel ghost point of phi, as we are using the
     *   Grad_par(phi) to calculate uEPar.
     * - phi is not an evolved variable.
     * - Setting a boundary condition on this would be to constrain the result.
     * - Extrapolation is used instead.
     */
    ownBC.extrapolateYGhost(phi);
    // Set BC on n and uIPar
    n    .applyBoundary();
    uIPar.applyBoundary();
    // Use the sheath boundary condition with constant Te for jPar and
    // uEPar at SE. Note that we have parallel uEPar derivatives
    // Here we have phiRef = Lambda
    ownBC.jParSheath (jPar, uEPar, uIPar, phi, n, Lambda, Lambda);
    ownBC.uEParSheath(uEPar, phi, Lambda, Lambda);
    // Set the BC for momDensPar
    ownBC.parDensMomSheath(momDensPar, uIPar, n);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(comGroup);
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(CelmaCurMom);
