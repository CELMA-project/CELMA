// *************** Simulation of CELMACURMOM *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "celmaCurMom.hxx"

// Initialization of the physics
// ############################################################################
int CelmaCurMom::init(bool restarting)
{
    TRACE("Halt in CelmaCurMom::init");

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Create the solver
    // ************************************************************************
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

    // Create the monitors
    // ************************************************************************
    ownMon.create(this);
    solver->addMonitor(ownMon.energyIntMon);
    // ************************************************************************

    // Get the constants
    // ************************************************************************
    // Get the section of the variables from [cst] specified in BOUT.inp
    // or in the command-line arguments
    Options *cst = options->getSection("cst");
    // Storing the variables with the following syntax
    // sectionName->get("nameInInp", varNameInCxx, defaultVal)
    cst->get("mu"                , mu                , 0.0);
    cst->get("Lambda"            , Lambda            , 0.0);
    cst->get("nuEI"              , nuEI              , 0.0);
    cst->get("nuEN"              , nuEN              , 0.0);
    cst->get("nuIN"              , nuIN              , 0.0);
    cst->get("artViscParLnN"     , artViscParLnN     , 0.0);
    cst->get("artViscParJpar"    , artViscParJpar    , 0.0);
    cst->get("artViscParMomDens" , artViscParMomDens , 0.0);
    cst->get("artViscParVortD"   , artViscParVortD   , 0.0);
    cst->get("artViscPerpLnN"    , artViscPerpLnN    , 0.0);
    cst->get("artViscPerpJPar"   , artViscPerpJPar   , 0.0);
    cst->get("artViscPerpMomDens", artViscPerpMomDens, 0.0);
    cst->get("artViscPerpVortD"  , artViscPerpVortD  , 0.0);
    cst->get("artHyperAzVortD"   , artHyperAzVortD   , 0.0);
    // ************************************************************************

    // Get the source constants
    // ************************************************************************
    Options *src = options->getSection("src");
    src->get("a",    a,    0.0);
    src->get("bRho", bRho, 0.0);
    src->get("bZ",   bZ,   0.0);
    src->get("cRho", cRho, 0.0);
    src->get("cZ",   cZ,   0.0);
    // ************************************************************************

    // Get the switches
    // ************************************************************************
    Options *switches = options->getSection("switch");
    switches->get("useHyperViscAzVortD", useHyperViscAzVortD, false);
    switches->get("includeNoise"       , includeNoise       , false);
    switches->get("forceAddNoise"      , forceAddNoise      , false);
    switches->get("saveDdt"            , saveDdt            , false);
    switches->get("saveTerms"          , saveTerms          , true );
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
    // Set artificial dissipation to 0 if useHyperViscAzVortD is false
    if (!useHyperViscAzVortD){
        artHyperAzVortD = 0.0;
    }
    // ************************************************************************

    // Calculate diffusion from grid size
    // ************************************************************************
    // SQ is squaring the expression
    // dx and dy are Field2D (0th index is ghost, but gives no problems)
    artViscParLnN     *= SQ(mesh->dy(0,0));
    artViscParJpar    *= SQ(mesh->dy(0,0));
    artViscParMomDens *= SQ(mesh->dy(0,0));
    artViscParVortD   *= SQ(mesh->dy(0,0));

    // The perpednicular diffusion is in our model
    /* NOTE: Chosen independent of dz
     *       This makes artVisc constant when expanding restarts
     */
    artViscPerpLnN     *= SQ(mesh->dx(0,0));
    artViscPerpJPar    *= SQ(mesh->dx(0,0));
    artViscPerpMomDens *= SQ(mesh->dx(0,0));
    artViscPerpVortD   *= SQ(mesh->dx(0,0));

    // Azimuthal hyperviscosities
    artHyperAzVortD *= SQ(SQ(mesh->dz));

    output << "Printing real artificial viscosity coefficients" << std::endl;
    output << "***********************************************" << std::endl;
    output << "Perpendicular (SQ(mesh->dx(0,0)) = "
              << SQ(mesh->dx(0,0)) << "):" << std::endl;
    output << "    For ln(n)    : "  << artViscPerpLnN     << std::endl;
    output << "    For j_{\\|}   : " << artViscPerpJPar    << std::endl;
    output << "    For nu_{i,\\|}: " << artViscPerpMomDens << std::endl;
    output << "    For vortD    : "  << artViscPerpVortD   << std::endl;
    output << "Parallel (SQ(mesh->dy(0,0)) = "
              << SQ(mesh->dy(0,0)) << "):" << std::endl;
    output << "    For ln(n)    : "  << artViscParLnN     << std::endl;
    output << "    For j_{\\|}   : " << artViscParJpar    << std::endl;
    output << "    For nu_{i,\\|}: " << artViscParMomDens << std::endl;
    output << "    For vortD    : "  << artViscParVortD   << std::endl;
    output << "Azimuthal hyperviscosity (SQ(SQ(mesh->dz)) = "
              << SQ(SQ(mesh->dz)) << "):" << std::endl;
    output << "    For vortD   : " << artHyperAzVortD << std::endl;
    output << "***********************************************" << std::endl;
    // ************************************************************************

    // Load from the geometry
    // ************************************************************************
    Options *geom = options->getSection("geom");
    geom->get("Lx", Lx, 0.0);
    geom->get("Ly", Ly, 0.0);
    // ************************************************************************

    // Additional fields
    // ************************************************************************
    // The source (obtained from the input)
    S = FieldFactory::get()
      ->create3D("theSource:S", Options::getRoot(), mesh, CELL_CENTRE, 0);

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

    // Split operator
    // ************************************************************************
    /* NOTE: Only active if a split solver is active
     *
     *       Convective is called first (solver.cxx run_rhs)
     *
     *       For implicit-explicit schemes, convective() will typically
     *       be treated explicitly, whilst diffusive() will be treated implicitly.
     *       For unsplit methods, both convective and diffusive will be called
     *       and the sum used to evolve the system:
     *       rhs() = convective() + diffusive()
     *       (comment in physicsmodel.hxx)
     */
    setSplitOperator();
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
    // Variables to be saved once
    // Constants
    SAVE_ONCE2(mu, Lambda);
    SAVE_ONCE3(nuEI, nuEN, nuIN);
    SAVE_ONCE5(a, bRho, bZ, cRho, cZ);
    SAVE_ONCE2(Lx, Ly);
    SAVE_ONCE4(artViscParLnN, artViscParJpar, artViscParMomDens, artViscParVortD);
    SAVE_ONCE4(artViscPerpLnN, artViscPerpJPar, artViscPerpMomDens, artViscPerpVortD);
    // Fields
    SAVE_ONCE(S);
    // Variables to be saved repeatedly
    // vort and phi
    SAVE_REPEAT2(vort, phi);
    if(saveTerms){
        // lnN terms
        SAVE_REPEAT3(lnNAdv, lnNRes, gradUEPar);
        SAVE_REPEAT4(lnNUeAdv, srcN, lnNParArtVisc, lnNPerpArtVisc);
        // jPar terms
        SAVE_REPEAT3(jParAdv, uIParAdvSum, uEParDoubleAdv);
        SAVE_REPEAT2(jParRes, elField);
        SAVE_REPEAT3(muElPressure, neutralERes, neutralIRes);
        SAVE_REPEAT2(jParParArtVisc, jParPerpArtVisc);
        // momDensPar terms
        SAVE_REPEAT2(momDensAdv, elPressure);
        SAVE_REPEAT3(densDiffusion, momDensParArtVisc, momDensPerpArtVisc);
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
    // Variables to be solved for
    SOLVE_FOR4(vortD, lnN, momDensPar, jPar);
    //*************************************************************************

    return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
/* The convective part:
 * Should contain all slow dynamics.
 * Will be stepped forward explicitly
 */
// NOTE: The convective part is called first, then the diffusive part
int CelmaCurMom::convective(BoutReal t)
{
    TRACE("Halt in CelmaCurMom::convective");

    if (includeNoise && !noiseAdded){
        /* NOTE: Positioning of includeNoise
         *
         * As the convective part is calculated first, the addnoise routine
         * will only be needed here.
         *
         * If the add noise is going to be called when doing a restart, the
         * includeNoise must be added in the rhs
         *
         * int main(int argc, char **argv)
         *    # Solver *solver = Solver::create();
         *    # solver->setModel(model);                  # Call of model init
         *    # solver->addMonitor(bout_monitor, Solver::BACK);
         *    # solver->outputVars(dump);                 # Call of load restart files
         *    # solver->solve();
         *int Solver::solve(int NOUT, BoutReal TIMESTEP)
         *int PvodeSolver::init(bool restarting, int nout, BoutReal tstep)
         *    # Around line 83 there is the call
         *    # if(Solver::init(restarting, nout, tstep)){
         *    # This loads the restart files
         *void solver_f(integer N, BoutReal t, N_Vector u, N_Vector udot, void *f_data
         *PvodeSolver::rhs(int N, BoutReal t, BoutReal *udata, BoutReal *dudata)
         *Solver::run_rhs(BoutReal t)
         *int PhysicsModel::runRHS(BoutReal time)
         *int CelmaCurMom::rhs(BoutReal t)
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
    // Use the sheath boundary condition with constant Te for jPar at SE
    // Here we have phiRef = Lambda
    ownBC.jParSheath (jPar, uEPar, uIPar, phi, n, Lambda, Lambda);
    ownBC.uEParSheath(uEPar, phi, Lambda, Lambda);
    // Set the BC for momDensPar
    ownBC.parDensMomSheath(momDensPar, uIPar, n);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(comGroup);

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
    ddt(lnN) =
          lnNAdv
        + gradUEPar
        + lnNUeAdv
        + srcN
        ;
    // Filtering highest modes
    ddt(lnN) = ownFilter->ownFilter(ddt(lnN));
    // ************************************************************************


    // Terms in jPar
    // ************************************************************************
    jParAdv        = - invJ*bracket(phi, jPar, bm);
    uIParAdvSum    = - Vpar_Grad_par(uIPar, n*(uIPar + uEPar));
    uEParDoubleAdv = 2.0*Vpar_Grad_par(uEPar, n*uEPar);
    jParRes        = - 0.51*nuEI*jPar;
    muElPressure   = mu*DDY(n);
    elField        = - mu*n*DDY(phi);
    neutralERes    = n*nuEN*uEPar;
    neutralIRes    = - n*nuIN*uIPar;

    ddt(jPar) =
          jParAdv
        + uIParAdvSum
        + uEParDoubleAdv
        + jParRes
        + muElPressure
        + elField
        + neutralERes
        + neutralIRes
        ;
    // Filtering highest modes
    ddt(jPar) = ownFilter->ownFilter(ddt(jPar));
    // ************************************************************************


    // Terms in momDensPar
    // ************************************************************************
    momDensAdv   = - invJ*bracket(phi, momDensPar, bm);
 // uIParAdvSum calculated above
    elPressure   = - DDY(n);
 // neutralIRes calculated above

    ddt(momDensPar) =
        + momDensAdv
        + uIParAdvSum
        + elPressure
        + neutralIRes
        ;
    // Filtering highest modes
    ddt(momDensPar) = ownFilter->ownFilter(ddt(momDensPar));
    // ************************************************************************


    // Preparation
    // ************************************************************************
    DivUIParNGradPerpPhi = ownOp->div_f_GradPerp_g(uIPar*n, phi);
    // Set the ghost points in order to take DDY
    ownBC.extrapolateYGhost(DivUIParNGradPerpPhi);
    // We must communicate as we will take DDY
    mesh->communicate(DivUIParNGradPerpPhi);
    // ************************************************************************


    // Terms in vorticity
    // ************************************************************************
    vortNeutral = - nuIN*n*vort;
    potNeutral  = - nuIN*ownOp->Grad_perp(phi)*ownOp->Grad_perp(n);
    vortDAdv    = - ownOp->vortDAdv(phi, vortD);
    kinEnAdvN   = - ownOp->kinEnAdvN(phi, n);
    parDerDivUIParNGradPerpPhi = - DDY(DivUIParNGradPerpPhi);
    divParCur                  =   DDY(jPar);

    ddt(vortD) =
          vortNeutral
        + potNeutral
        + parDerDivUIParNGradPerpPhi
        + divParCur
        + vortDAdv
        + kinEnAdvN
        ;

    // Filtering highest modes
    ddt(vortD) = ownFilter->ownFilter(ddt(vortD));
    // ************************************************************************
    return 0;
}

/* The diffusive part:
 * Should contain all the fast dynamics.
 *
 * What are the fast dynamics
 * --------------------------
 * Other than knowing that diffusion have infinite fast dynamics (only limited
 * by the grid size), and that higher order terms are usually stiff, it is hard
 * to know apriori what terms will contribute to the fast part of the dynamics.
 * What one can do then, is in principle to find the max(ddt(field)), and see
 * if there are any fields which are  much faster than others. These fields can
 * then be put into the diffusive part.
 *
 * The diffusive part will be treated implicitly. This part is usually called
 * several times during a timestep. as there usually will be several newton
 * iterations pr timestep. Therefore, one should have the inversion algorithm
 * here as well
 */
int CelmaCurMom::diffusive(BoutReal t, bool linear)
{
    TRACE("Halt in CelmaCurMom::diffusive");

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
    // uEPar at SE
    // Here we have phiRef = Lambda
    ownBC.jParSheath (jPar, uEPar, uIPar, phi, n, Lambda, Lambda);
    ownBC.uEParSheath(uEPar, phi, Lambda, Lambda);
    // Set the BC for momDensPar
    ownBC.parDensMomSheath(momDensPar, uIPar, n);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(comGroup);

    // Terms in lnNPar
    // ************************************************************************
    lnNRes         = (0.51*nuEI/mu)*(Laplace_perp(lnN)+gradPerpLnN*gradPerpLnN);
    lnNParArtVisc  = artViscParLnN*D2DY2(lnN);
    lnNPerpArtVisc = artViscPerpLnN*Laplace_perp(lnN);
    ddt(lnN) =
          lnNRes
        + lnNParArtVisc
        + lnNPerpArtVisc
        ;
    // Filtering highest modes
    ddt(lnN) = ownFilter->ownFilter(ddt(lnN));
    // ************************************************************************


    // Terms in jPar
    // ************************************************************************
    jParParArtVisc  = (artViscParJpar)*D2DY2(jPar);
    jParPerpArtVisc = (artViscPerpJPar)*Laplace_perp(jPar);

    ddt(jPar) =
          jParParArtVisc
        + jParPerpArtVisc
        ;
    // Filtering highest modes
    ddt(jPar) = ownFilter->ownFilter(ddt(jPar));
    // ************************************************************************


    // Terms in momDensPar
    // ************************************************************************
    densDiffusion      = 0.51*nuEI*(uIPar/mu)*Laplace_perp(n);
    momDensParArtVisc  = (artViscParMomDens)*D2DY2(momDensPar);
    momDensPerpArtVisc = (artViscPerpMomDens)*Laplace_perp(momDensPar);

    ddt(momDensPar) =
          densDiffusion
        + momDensParArtVisc
        + momDensPerpArtVisc
        ;
    // Filtering highest modes
    ddt(momDensPar) = ownFilter->ownFilter(ddt(momDensPar));
    // ************************************************************************


    // Terms in vorticity
    // ************************************************************************
    vortDParArtVisc  =   artViscParVortD*D2DY2(vortD);
    vortDPerpArtVisc =   artViscPerpVortD*Laplace_perp(vortD);
    vortDhyperVisc   = - artHyperAzVortD*D4DZ4(vortD);

    ddt(vortD) =
          vortDParArtVisc
        + vortDPerpArtVisc
        + vortDhyperVisc
        ;
    // Filtering highest modes
    ddt(vortD) = ownFilter->ownFilter(ddt(vortD));
    // ************************************************************************
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(CelmaCurMom);
