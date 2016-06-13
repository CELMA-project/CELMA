// *************** Simulation of CELMA *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "celma.hxx"

// Initialization of the physics
// ############################################################################
int Celma::init(bool restarting) {
    TRACE("Halt in Celma::init");

    // Create the solver
    // ************************************************************************
    /* NOTE: Calls createOperators without making an object of OwnOperators.
     *       The child is typecasted to the parent
     */
    ownOp = OwnOperators::createOperators();
    ownLapl.create(ownOp, ownBC);
    // ************************************************************************

    // Prepare cauchy
    // ************************************************************************
    ownBC.prepareCauchy("lnN");
    // ************************************************************************

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Get the constants
    // ************************************************************************
    // Get the section of the variables from [cst] specified in BOUT.inp
    // or in the command-line arguments
    Options *cst = options->getSection("cst");
    // Storing the variables with the following syntax
    // sectionName->get("nameInInp", varNameInCxx, defaultVal)
    cst->get("mu",               mu,               0.0);
    cst->get("Lambda",           Lambda,           0.0);
    cst->get("nuEI",             nuEI,             0.0);
    cst->get("nuEN",             nuEN,             0.0);
    cst->get("nuIN",             nuIN,             0.0);
    cst->get("artViscParLnN",    artViscParLnN,    0.0);
    cst->get("artViscParUEPar",  artViscParUEPar,  0.0);
    cst->get("artViscParUIPar",  artViscParUIPar,  0.0);
    cst->get("artViscParVortD",  artViscParVortD,  0.0);
    cst->get("artViscPerpLnN",   artViscPerpLnN,   0.0);
    cst->get("artViscPerpUEPar", artViscPerpUEPar, 0.0);
    cst->get("artViscPerpUIPar", artViscPerpUIPar, 0.0);
    cst->get("artViscPerpVortD", artViscPerpVortD, 0.0);
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

    // Calculate diffusion from grid size
    // ************************************************************************
    // SQ is squaring the expression
    // dx and dy are Field2D (0th index is ghost, but gives no problems)
    artViscParLnN   *= SQ(mesh->dy(0,0));
    artViscParUEPar *= SQ(mesh->dy(0,0));
    artViscParUIPar *= SQ(mesh->dy(0,0));
    artViscParVortD *= SQ(mesh->dy(0,0));

    // The perpednicular diffusion is in our model
    artViscPerpLnN   *= SQ(mesh->dx(0,0) + mesh->dz);
    artViscPerpUEPar *= SQ(mesh->dx(0,0) + mesh->dz);
    artViscPerpUIPar *= SQ(mesh->dx(0,0) + mesh->dz);
    artViscPerpVortD *= SQ(mesh->dx(0,0) + mesh->dz);
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

    // The radial profile (obtained from the input)
    dampingProfile = FieldFactory::get()
      ->create3D("dampProf:profile", Options::getRoot(), mesh, CELL_CENTRE, 0);

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
    comGroup.add(uIPar);
    comGroup.add(uEPar);
    comGroup.add(phi);
    comGroup.add(vortD);
    // ************************************************************************

    // Specify what values should be stored in the .dmp file
    // ************************************************************************
    // Variables to be saved once
    // Constants
    SAVE_ONCE2(mu, Lambda);
    SAVE_ONCE3(nuEI, nuEN, nuIN);
    SAVE_ONCE5(a, bRho, bZ, cRho, cZ);
    SAVE_ONCE2(Lx, Ly);
    SAVE_ONCE4(artViscParLnN, artViscParUEPar, artViscParUIPar, artViscParVortD);
    SAVE_ONCE4(artViscPerpLnN, artViscPerpUEPar, artViscPerpUIPar, artViscPerpVortD);
    // Fields
    SAVE_ONCE2(S, dampingProfile);
    // Variables to be saved repeatedly
    // vort and phi
    SAVE_REPEAT2(vort, phi);
    // lnN terms
    SAVE_REPEAT3(lnNAdv, lnNRes, gradUEPar);
    SAVE_REPEAT4(lnNUeAdv, srcN, lnNParArtVisc, lnNPerpArtVisc);
    // uEPar terms
    SAVE_REPEAT3(uEParAdv, uEParParAdv, gradPhiLnN);
    SAVE_REPEAT4(uEParRes, ueSrc, uEParParArtVisc, uEParPerpArtVisc);
    SAVE_REPEAT (ueNeutral);
    // uIPar terms
    SAVE_REPEAT3(uIParAdv, uIParParAdv, gradPhi);
    SAVE_REPEAT4(uIParRes, uiSrc, uIParParArtVisc, uIParPerpArtVisc);
    SAVE_REPEAT (uiNeutral);
    // Vorticity terms
    SAVE_REPEAT2(vortNeutral, potNeutral);
    SAVE_REPEAT2(divExBAdvGradPerpPhiN, parDerDivUIParNGradPerpPhi)
    SAVE_REPEAT4(nGradUiUe, uiUeGradN, vortDParArtVisc, vortDPerpArtVisc);
    // Variables to be solved for
    SOLVE_FOR4(vortD, lnN, uIPar, uEPar);
    //*************************************************************************
    return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int Celma::rhs(BoutReal t) {
    TRACE("Halt in Celma::rhs");

    // Manually specifying rho inner ghost points
    // ************************************************************************
    ownBC.innerRhoCylinder(lnN);
    ownBC.innerRhoCylinder(uEPar);
    ownBC.innerRhoCylinder(uIPar);
    ownBC.innerRhoCylinder(phi);    // Used to set BC in lapalace inversion
    ownBC.innerRhoCylinder(vortD);  // Later taken derivative of
    // The inner boundaries of phi is set in the inversion procedure
    // ************************************************************************

    // Preparations
    // ************************************************************************
    /* NOTE: Preparations
     * 1. gradPerp(lnN) is input to the NaulinSolver , thus
     *    a) gradPerpLnN needs to be calculated...
     *    b) ... which requires that the boundary values are set
     *    c) ... and that lnN has been communicated
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
    // Use caucy boundary condition for the density
    ownBC.cauchyYDown(lnN);
    mesh->communicate(lnN);
    gradPerpLnN = ownOp->Grad_perp(lnN);
    n = DDY(lnN);
    n = D2DY2(lnN);
    n = exp(lnN);
    // ************************************************************************

    // Laplace inversion
    // ************************************************************************
    /* NOTE: Solve for phi
     * 1. phi and vort is output
     * 2. The solver takes care of perpendicular boundaries
     */
    phi = ownLapl.NaulinSolver(gradPerpLnN, n, vortD, phi, vort);
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
    // Use the sheath boundary condition with constant Te for uEPar at SE
    // Here we have phiRef = Lambda
    ownBC.uEParSheath(uEPar, phi, Lambda, Lambda, dampingProfile);
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
    lnNRes         =   (0.51*nuEI/mu) *
                         (Laplace_perp(lnN)+gradPerpLnN*gradPerpLnN);
    gradUEPar      = - Grad_par(uEPar);
    lnNUeAdv       = - Vpar_Grad_par(uEPar, lnN);
    srcN           =   S/n;
    lnNParArtVisc  =   artViscParLnN*D2DY2(lnN);
    lnNPerpArtVisc =   artViscPerpLnN*Laplace_perp(lnN);
    ddt(lnN) =   lnNAdv
                + lnNRes
                + gradUEPar
                + lnNUeAdv
                + srcN
                + lnNParArtVisc
                + lnNPerpArtVisc
                ;
    // ************************************************************************


    // Terms in uEPar
    // ************************************************************************
    uEParAdv         = - invJ*bracket(phi, uEPar, bm) ;
    uEParParAdv      = - Vpar_Grad_par(uEPar, uEPar) ;
    gradPhiLnN       =   mu*Grad_par(phi - lnN) ;
    uEParRes         = - 0.51*nuEI *( uEPar - uIPar) ;
    ueNeutral        = - nuEN*uEPar;
    ueSrc            = - S*uEPar/n ;
    uEParParArtVisc  =   artViscParUEPar*D2DY2(uEPar)/n;
    uEParPerpArtVisc =   (artViscPerpUEPar/n)*Laplace_perp(uEPar);

    ddt(uEPar) =
          uEParAdv
        + uEParParAdv
        + gradPhiLnN
        + uEParRes
        + ueNeutral
        + ueSrc
        + uEParParArtVisc
        + uEParPerpArtVisc
        ;
    // ************************************************************************


    // Terms in uIPar
    // ************************************************************************
    uIParAdv         = - invJ*bracket(phi, uIPar, bm);
    uIParParAdv      = - Vpar_Grad_par(uIPar, uIPar);
    gradPhi          = - Grad_par(phi);
    uIParRes         = - (0.51*nuEI/mu) *(uIPar - uEPar);
    uiNeutral        = - nuIN*uIPar;
    uiSrc            = - S*uIPar/n;
    uIParParArtVisc  =   artViscParUIPar*D2DY2(uIPar)/n;
    uIParPerpArtVisc =   (artViscPerpUIPar/n)*Laplace_perp(uIPar);

    ddt(uIPar) =
          uIParAdv
        + uIParParAdv
        + gradPhi
        + uIParRes
        + uiNeutral
        + uiSrc
        + uIParParArtVisc
        + uIParPerpArtVisc
        ;
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
    vortNeutral                = - nuIN*n*vort;
    potNeutral                 = - nuIN*ownOp->Grad_perp(phi)*ownOp->Grad_perp(n);
    divExBAdvGradPerpPhiN      = - ownOp->div_uE_dot_grad_n_GradPerp_phi(n, phi);
    parDerDivUIParNGradPerpPhi = - DDY(DivUIParNGradPerpPhi);
    nGradUiUe                  =   Vpar_Grad_par(n, uIPar - uEPar);
    uiUeGradN                  =   Vpar_Grad_par(uIPar - uEPar, n);
    vortDParArtVisc            =   artViscParVortD*D2DY2(vortD);
    vortDPerpArtVisc           =   artViscPerpVortD*Laplace_perp(vortD);

    ddt(vortD) =
              vortNeutral
            + potNeutral
            + divExBAdvGradPerpPhiN
            + parDerDivUIParNGradPerpPhi
            + nGradUiUe
            + uiUeGradN
            + vortDParArtVisc
            + vortDPerpArtVisc
            ;
    // ************************************************************************
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(Celma);
