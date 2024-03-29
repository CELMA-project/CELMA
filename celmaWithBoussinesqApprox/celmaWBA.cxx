// ****************** Celma with Boussinesq approximation *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "celmaWBA.hxx"

// Initialization and solving of the physics
// ############################################################################
int CelmaWBA::init(bool restarting) {
  TRACE("Halt in CelmaWBA::init");

  Coordinates *coord = mesh->getCoordinates();

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
  phi = FieldFactory::get()->create3D("phi:function", Options::getRoot(), mesh,
                                      CELL_CENTRE, 0);

  // The metric coefficient (needed in front of the arakawa bracket)
  invJ = (1.0 / coord->J);
  // ************************************************************************

  // Specifying the brackets to the arakawa scheme
  // (see examples/MMS/hw for a elegant way to choose from the input file)
  // ************************************************************************
  bm = BRACKET_ARAKAWA;
  // ************************************************************************

  // Add a FieldGroup to communicate
  // ************************************************************************
  // NOTE: We only communicate variables we are taking derivatives of
  // NOTE: vort is communicated when taking the laplace inversion
  comGroup.add(lnN);
  comGroup.add(n);
  comGroup.add(momDensPar);
  comGroup.add(jPar);
  comGroup.add(uEPar);
  comGroup.add(uIPar);
  comGroup.add(phi);
  // ************************************************************************

  // Specify BC for n and uIPar (used to set BC for jPar and momDensPar)
  // ************************************************************************
  n.setBoundary("n");
  uIPar.setBoundary("uIPar");
  // ************************************************************************

  // Specify what values should be stored in the .dmp file
  // ************************************************************************
  // Variables to be saved repeatedly
  SAVE_REPEAT(phi);
  if (saveTerms) {
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
    SAVE_REPEAT4(divParCur, divSourcePhi, vortParArtVisc, vortPerpArtVisc);
    SAVE_REPEAT3(DDYGradPerpPhiGradPerpUI, vortAdv, vortParAdv);
    // Helping fields
    SAVE_REPEAT2(uIPar, uEPar);

    if (saveDdt) {
      SAVE_REPEAT4(ddt(vort), ddt(lnN), ddt(momDensPar), ddt(jPar));
    }
  }
  // Monitor variables to be solved for
  if (monitorEnergy) {
    for (std::map<std::string, BoutReal>::iterator it = kinE.begin();
         it != kinE.end(); ++it) {
      dump.add(it->second, it->first.c_str(), 1);
    }
  }
  if (monitorParticleNumber) {
    for (std::map<std::string, BoutReal>::iterator it = particleNumber.begin();
         it != particleNumber.end(); ++it) {
      dump.add(it->second, it->first.c_str(), 1);
    }
  }
  // Variables to be solved for
  SOLVE_FOR4(vort, lnN, momDensPar, jPar);
  //*************************************************************************

  return 0;
}

int CelmaWBA::rhs(BoutReal t) {
  TRACE("Halt in CelmaWBA::rhs");

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
  lnNAdv = -invJ * bracket(phi, lnN, bm);
  gradUEPar = -DDY(uEPar);
  lnNUeAdv = -Vpar_Grad_par(uEPar, lnN);
  srcN = S / n;
  lnNRes = (nuEI / mu) * (Laplace_perp(lnN) + gradPerpLnN * gradPerpLnN);
  lnNParArtVisc = artViscParLnN * D2DY2(lnN);
  lnNPerpArtVisc = artViscPerpLnN * Laplace_perp(lnN);
  ddt(lnN) = lnNAdv + gradUEPar + lnNUeAdv + srcN + lnNRes + lnNParArtVisc +
             lnNPerpArtVisc;
  // Filtering highest modes
  ddt(lnN) = ownFilter->ownFilter(ddt(lnN));
  // ************************************************************************

  // Terms in jPar
  // ************************************************************************
  jParAdv = -invJ * bracket(phi, jPar, bm);
  uEParAdv = n * (Vpar_Grad_par(uEPar, uEPar));
  uIParAdv = -n * (Vpar_Grad_par(uIPar, uIPar));
  jParParAdv = -Vpar_Grad_par((jPar / n), uEPar);
  gradPhiLnN = mu * n * DDY(lnN - phi);
  jParRes = -0.51 * nuEI * jPar;
  neutralERes = n * nuEN * uEPar;
  neutralIRes = -n * nuIN * uIPar;
  jParParArtVisc = (artViscParJpar)*D2DY2(jPar);
  jParPerpArtVisc = (artViscPerpJPar)*Laplace_perp(jPar);

  ddt(jPar) = jParAdv + uEParAdv + uIParAdv + jParParAdv + gradPhiLnN +
              jParRes + neutralERes + neutralIRes + jParParArtVisc +
              jParPerpArtVisc;
  // Filtering highest modes
  ddt(jPar) = ownFilter->ownFilter(ddt(jPar));
  // ************************************************************************

  // Terms in momDensPar
  // ************************************************************************
  momDensAdv = -invJ * bracket(phi, momDensPar, bm);
  // UIParAdv calculated above
  tmp = n * uEPar;
  mesh->communicate(tmp);
  uIFluxAdv = -Vpar_Grad_par(uIPar, tmp);
  elPressure = -DDY(n);
  densDiffusion = nuEI * (uIPar / mu) * Laplace_perp(n);
  // neutralIRes calculated above
  neutralEResMu = neutralERes / mu;
  momDensParArtVisc = (artViscParMomDens)*D2DY2(momDensPar);
  momDensPerpArtVisc = (artViscPerpMomDens)*Laplace_perp(momDensPar);

  ddt(momDensPar) = +momDensAdv + uIFluxAdv + uIParAdv + elPressure +
                    densDiffusion + neutralIRes + neutralEResMu +
                    momDensParArtVisc + momDensPerpArtVisc;
  // Filtering highest modes
  ddt(momDensPar) = ownFilter->ownFilter(ddt(momDensPar));
  // ************************************************************************

  // Preparation
  // ************************************************************************
  // Saving gradPerpPhi for use in monitors
  gradPerpPhi = ownOp->Grad_perp(phi);
  // Set the ghost points in order to take DDY
  ownBC.extrapolateYGhost(gradPerpPhi.x);
  ownBC.extrapolateYGhost(gradPerpPhi.y);
  ownBC.extrapolateYGhost(gradPerpPhi.z);
  // Set the ghost points in order to take gradPerp
  ownBC.innerRhoCylinder(gradPerpPhi.x);
  ownBC.innerRhoCylinder(gradPerpPhi.y);
  ownBC.innerRhoCylinder(gradPerpPhi.z);
  ownBC.extrapolateXOutGhost(gradPerpPhi.x);
  ownBC.extrapolateXOutGhost(gradPerpPhi.y);
  ownBC.extrapolateXOutGhost(gradPerpPhi.z);
  mesh->communicate(gradPerpPhi);
  // ************************************************************************

  // Terms in vorticity
  // ************************************************************************
  vortNeutral = -nuIN * n * vort;
  potNeutral = -nuIN * gradPerpPhi * ownOp->Grad_perp(n);
  // NOTE: Basic brackets => vortDAdv is normal bracket adv
  vortAdv = -ownOp->vortDAdv(phi, vort);
  // NOTE: Upwinding could be used on this, but use with care as
  //       dissipation may introduce spurious vorticity
  vortParAdv = -uIPar * DDY(vort);
  DDYGradPerpPhiGradPerpUI = -ownOp->DDY(gradPerpPhi) * ownOp->Grad_perp(uIPar);
  divSourcePhi = -S * vort - gradPerpPhi * ownOp->Grad_perp(S);
  divParCur = DDY(jPar);
  vortParArtVisc = artViscParVort * D2DY2(vort);
  vortPerpArtVisc = artViscPerpVort * Laplace_perp(vort);
  vortHyperVisc = -artHyperAzVort * D4DZ4(vort);

  ddt(vort) = vortNeutral + potNeutral + vortAdv + vortParAdv + divParCur +
              DDYGradPerpPhiGradPerpUI + divSourcePhi + vortParArtVisc +
              vortPerpArtVisc + vortHyperVisc;

  // Filtering highest modes
  ddt(vort) = ownFilter->ownFilter(ddt(vort));
  // ************************************************************************
  return 0;
}
// ############################################################################

// Constructor
// ############################################################################
CelmaWBA::CelmaWBA()
    : kinE{{"perpKinEE", 0.0},
           {"parKinEE", 0.0},
           {"perpKinEI", 0.0},
           {"parKinEI", 0.0}},
      particleNumber{{"particleNumber", 0.0}} {
  TRACE("Halt in CelmaWBA::CelmaWBA");
}
// ############################################################################

// Monitor
// ############################################################################
int CelmaWBA::outputMonitor(BoutReal simtime, int iter, int NOUT) {
  TRACE("Halt in CelmaWBA::outputMonitor");

  if (monitorEnergy) {
    ownMon.kinEnergy(n, gradPerpPhi, uEPar, uIPar, &kinE);
  }
  if (monitorParticleNumber) {
    ownMon.numberOfParticles(n, &particleNumber);
  }

  return 0;
}
// ############################################################################

// Initialization helpers
// ############################################################################
void CelmaWBA::initializeOwnObjects() {
  TRACE("CelmaWBA::initializeOwnObjects");

  // Create OwnOperators
  // ************************************************************************
  // Get the option (before any sections) in the BOUT.inp file
  Options *options = Options::getRoot();

  /* NOTE: Calls createOperators without making an object of OwnOperators.
   *       The child is typecasted to the parent
   */
  ownOp = OwnOperators::createOperators();

  // Create the laplace object
  // ************************************************************************
  // The laplace object will look in the section phiSolver in the BOUT.inp
  // file
  Options *phiSol_opt = options->getSection("phiSolver");
  phiSolver = Laplacian::create(phiSol_opt);
  // Set the coefficients manually (should also be set to this by default)
  phiSolver->setCoefD(1.0);
  phiSolver->setCoefC(1.0);
  phiSolver->setCoefA(0.0);
  // ************************************************************************

  // As a small hack we find the own operator type to figure out how to
  // treat the div(ue*grad(n grad_perp phi))
  Options *ownOp = options->getSection("ownOperators");
  ownOp->get("type", ownOpType, "BasicBrackets");
  if ((lowercase(ownOpType) == lowercase("simpleStupid")) ||
      (lowercase(ownOpType) == lowercase("onlyBracket"))) {
    throw BoutException(
        "simpleStupid and onlyBracket not available in this model");
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

void CelmaWBA::setAndSaveParameters() {
  TRACE("Halt in CelmaWBA::setAndSaveParameters");

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
  input->get("n0", n0, 0.0);
  input->get("Te0", Te0, 0.0);
  input->get("Ti0", Ti0, 0.0);
  input->get("B0", B0, 0.0);
  input->get("Sn", Sn, 0.0);
  input->get("nn", nn, 0.0);
  input->get("gas", gas, "H");
  // Cast to upper
  gas[0] = toupper(gas[0]);

  input->get("warningForException", warningForException, false);

  Parameters params(radius, length, n0, Te0, Ti0, B0, Sn, nn, gas,
                    warningForException);
  // ************************************************************************

  // Get the variables
  // ************************************************************************
  LxParam = params.getLx();
  LyParam = params.getLy();
  nuEI = params.getNuEINorm();
  nuEN = params.getNuENNorm();
  nuIN = params.getNuINNorm();
  SNorm = params.getSNorm();
  mu = params.getMu();
  Lambda = params.getLambda();
  beta = params.getBeta();
  omCI = params.getOmCI();
  rhoS = params.getRhoS();
  eta0INorm = params.getEta0INorm();
  // ************************************************************************

  // Check that Lx and LxParams is the same up until the fourt decimal point
  // ************************************************************************
  std::ostringstream stream;
  int precision = 4;
  bool throwError = false;
  if ((fabs(round(LxParam * 1e4) / 1e4 - round(Lx * 1e4) / 1e4)) >
      DBL_EPSILON) {
    stream << "Mismatch between 'Lx' calculated from 'radius' "
           << "and input 'Lx'\n"
           << "Calculated = " << std::fixed << std::setprecision(precision)
           << round(LxParam * 1e4) / 1e4 << "\n"
           << "Input      = " << std::fixed << std::setprecision(precision)
           << round(Lx * 1e4) / 1e4 << "\n\n";
    throwError = true;
  }
  if ((fabs(round(LyParam * 1e4) / 1e4 - round(Ly * 1e4) / 1e4)) >
      DBL_EPSILON) {
    stream << "Mismatch between 'Ly' calculated from 'length' "
           << "and input 'Ly'\n"
           << "Calculated = " << std::fixed << std::setprecision(precision)
           << round(LyParam * 1e4) / 1e4 << "\n"
           << "Input      = " << std::fixed << std::setprecision(precision)
           << round(Ly * 1e4) / 1e4 << "\n\n";
    throwError = true;
  }
  if (throwError) {
    std::string str = stream.str();
    // Cast the stream to a const char in order to use it in BoutException
    const char *message = str.c_str();
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
  SAVE_ONCE(beta);
  SAVE_ONCE2(omCI, rhoS);
  SAVE_ONCE(eta0INorm);
  // ************************************************************************
}

void CelmaWBA::printPointsPerRhoS() {
  TRACE("Halt in CelmaWBA::printPointsPerRhoS");

  // Get the option (before any sections) in the BOUT.inp file
  int MXG;
  BoutReal pointsPerRhoSRadially;
  BoutReal pointsPerRhoSParallely;
  BoutReal pointsPerRhoSAzimuthally;
  BoutReal minPointsPerRhoSXZ;
  BoutReal minPointsPerRhoSY;

  Options *root = Options::getRoot();
  root->get("MXG", MXG, 0);

  if (MXG == 0) {
    throw BoutException("No MXG found, please specify.");
  }

  Coordinates *coord = mesh->getCoordinates();
  // dx = Lx/(nx-2*MXG) => nx = (Lx/dx) + 2*MXG
  pointsPerRhoSRadially = ((Lx / coord->dx(0, 0)) + 2 * MXG) / Lx;
  // dy = Ly/ny => ny = Ly/dy => ny/Ly = 1/dy
  pointsPerRhoSParallely = 1.0 / coord->dy(0, 0);
  // O=2*pi*r, so on edge nz/rho_s = nz/(2*pi*Lx)
  pointsPerRhoSAzimuthally = (mesh->LocalNz) / (2.0 * PI * Lx);

  root->getSection("geom")->get("minPointsPerRhoSXZ", minPointsPerRhoSXZ, 3.0);
  root->getSection("geom")->get("minPointsPerRhoSY", minPointsPerRhoSY, 1.0e-1);

  std::ostringstream stream;
  bool throwError = false;
  if (pointsPerRhoSRadially < minPointsPerRhoSXZ) {
    stream << "Minimum points per rhoS not fulfilled in x.\n"
           << "Limit is         " << minPointsPerRhoSXZ << "\n"
           << "Current value is " << pointsPerRhoSRadially << "\n\n";
    throwError = true;
  }
  if (mesh->LocalNz > 2 && pointsPerRhoSAzimuthally < minPointsPerRhoSXZ) {
    stream << "Minimum points per rhoS not fulfilled on outer circumference.\n"
           << "Limit is         " << minPointsPerRhoSXZ << "\n"
           << "Current value is " << pointsPerRhoSAzimuthally << "\n\n";
    throwError = true;
  }
  if (pointsPerRhoSParallely < minPointsPerRhoSY) {
    stream << "Minimum points per rhoS not fulfilled in y.\n"
           << "Limit is         " << minPointsPerRhoSY << "\n"
           << "Current value is " << pointsPerRhoSParallely << "\n\n";
    throwError = true;
  }
  if (throwError) {
    std::string str = stream.str();
    // Cast the stream to a const char in order to use it in BoutException
    const char *message = str.c_str();
    throw BoutException(message);
  }

  output << '\n' << std::string(51, '*') << std::endl;
  output << "Points per rhoS in x                   = " << pointsPerRhoSRadially
         << std::endl;
  output << "Points per rhoS on outer circumference = "
         << pointsPerRhoSAzimuthally << std::endl;
  output << "Points per rhoS in y                   = "
         << pointsPerRhoSParallely << std::endl;
  output << std::string(51, '*') << '\n' << std::endl;
}

void CelmaWBA::setAndSaveSource() {
  TRACE("Halt in CelmaWBA::setAndSaveSource");

  BoutReal radialWidth, radialCentre, radialSteepness;
  BoutReal parallelWidth, parallelCentre, parallelSteepness;

  // Get the source constants, save them and create the source
  // ************************************************************************
  // Get the option (before any sections) in the BOUT.inp file
  Options *options = Options::getRoot();

  Options *thesource = options->getSection("thesource");
  thesource->get("radialWidth", radialWidth, 0.0);
  thesource->get("radialCentre", radialCentre, 0.0);
  thesource->get("radialSteepness", radialSteepness, 0.0);
  thesource->get("parallelWidth", parallelWidth, 0.0);
  thesource->get("parallelCentre", parallelCentre, 0.0);
  thesource->get("parallelSteepness", parallelSteepness, 0.0);

  SAVE_ONCE3(radialWidth, radialCentre, radialSteepness);
  SAVE_ONCE3(parallelWidth, parallelCentre, parallelSteepness);

  // The source (obtained from the input, and multiplied with SNorm)
  S = SNorm *
      FieldFactory::get()->create3D("theSource:S", Options::getRoot(), mesh,
                                    CELL_CENTRE, 0);

  // Save the source field
  SAVE_ONCE(S);
  // ************************************************************************
}

void CelmaWBA::setSwithces(bool &restarting) {
  TRACE("Halt in CelmaWBA::setSwithces");

  // Get the switches
  // ************************************************************************
  // Get the option (before any sections) in the BOUT.inp file
  Options *options = Options::getRoot();

  Options *switches = options->getSection("switch");
  switches->get("useHyperViscAzVort", useHyperViscAzVort, false);
  switches->get("saveDdt", saveDdt, false);
  switches->get("constViscPar", constViscPar, false);
  switches->get("constViscPerp", constViscPerp, false);
  switches->get("constViscHyper", constViscHyper, false);
  switches->get("saveTerms", saveTerms, true);
  switches->get("monitorEnergy", monitorEnergy, true);
  switches->get("monitorParticleNumber", monitorParticleNumber, true);
  switches->get("viscosityGuard", viscosityGuard, true);
  // ************************************************************************
}

void CelmaWBA::setAndSaveViscosities() {
  TRACE("Halt in CelmaWBA::setAndSaveViscosities");

  Coordinates *coord = mesh->getCoordinates();

  // Get and save the viscosities
  // ************************************************************************
  // Get the option (before any sections) in the BOUT.inp file
  Options *options = Options::getRoot();

  Options *visc = options->getSection("visc");
  visc->get("artViscParLnN", artViscParLnN, 0.0);
  visc->get("artViscParJpar", artViscParJpar, 0.0);
  visc->get("artViscParMomDens", artViscParMomDens, 0.0);
  visc->get("artViscParVort", artViscParVort, 0.0);
  visc->get("artViscPerpLnN", artViscPerpLnN, 0.0);
  visc->get("artViscPerpJPar", artViscPerpJPar, 0.0);
  visc->get("artViscPerpMomDens", artViscPerpMomDens, 0.0);
  visc->get("artViscPerpVort", artViscPerpVort, 0.0);
  visc->get("artHyperAzVort", artHyperAzVort, 0.0);

  // Calculate diffusion from grid size
  if (!constViscPar) {
    // SQ is squaring the expression
    // dx and dy are Field2D (0th index is ghost, but gives no problems)
    artViscParLnN *= SQ(coord->dy(0, 0));
    artViscParJpar *= SQ(coord->dy(0, 0));
    artViscParMomDens *= SQ(coord->dy(0, 0));
    artViscParVort *= SQ(coord->dy(0, 0));
  }

  if (!constViscPerp) {
    // The perpednicular diffusion is in our model
    /* NOTE: Chosen independent of dz
     *       This makes artVisc constant when expanding restarts
     */
    artViscPerpLnN *= SQ(coord->dx(0, 0));
    artViscPerpJPar *= SQ(coord->dx(0, 0));
    artViscPerpMomDens *= SQ(coord->dx(0, 0));
    artViscPerpVort *= SQ(coord->dx(0, 0));
  }

  // Set artificial viscosities to 0 if useHyperViscAzVort is false
  if (!useHyperViscAzVort) {
    output << "Setting artHyperAzVort = 0.0 as useHyperViscAzVort = False"
           << std::endl;
    artHyperAzVort = 0.0;
  }
  if (!constViscHyper) {
    // Azimuthal hyperviscosities
    artHyperAzVort *= SQ(SQ(coord->dz));
  }

  // Print and store the variables
  output << "\nPrinting real artificial viscosity coefficients" << std::endl;
  output << "***********************************************" << std::endl;
  output << "Perpendicular";
  if (!constViscPerp) {
    output << " (SQ(coord->dx(0,0)) = "
           << SQ(coord->dx(0, 0)) << "):";
  }
  output << std::endl;
  output << "    For ln(n)    : " << artViscPerpLnN << std::endl;
  output << "    For j_{\\|}   : " << artViscPerpJPar << std::endl;
  output << "    For nu_{i,\\|}: " << artViscPerpMomDens << std::endl;
  output << "    For vort     : " << artViscPerpVort << std::endl;
  output << "Parallel";
  if (!constViscPar) {
    output << " (SQ(coord->dy(0,0)) = "
           << SQ(coord->dy(0, 0)) << "):";
  }
  output << std::endl;
  output << "    For ln(n)    : " << artViscParLnN << std::endl;
  output << "    For j_{\\|}   : " << artViscParJpar << std::endl;
  output << "    For nu_{i,\\|}: " << artViscParMomDens << std::endl;
  output << "    For vort     : " << artViscParVort << std::endl;
  output << "Azimuthal hyperviscosity";
  if (!constViscHyper) {
    output << "Azimuthal hyperviscosity (SQ(SQ(coord->dz)) = "
           << SQ(SQ(coord->dz)) << "):";
  }
  output << std::endl;
  output << "    For vort    : " << artHyperAzVort << std::endl;
  output << "***********************************************\n" << std::endl;

  // Guards
  std::vector<BoutReal> perpVisc(4, 0.0);
  perpVisc.push_back(artViscPerpLnN);
  perpVisc.push_back(artViscPerpJPar);
  perpVisc.push_back(artViscPerpMomDens);
  perpVisc.push_back(artViscPerpVort);

  for (std::vector<BoutReal>::iterator it = perpVisc.begin();
       it != perpVisc.end(); ++it) {
    if (*it > nuEI) {
      throw BoutException("One of the perpendicular viscosities "
                          "is larger than nuEI");
    }
  }

  std::vector<BoutReal> parVisc(4, 0.0);
  perpVisc.push_back(artViscParLnN);
  perpVisc.push_back(artViscParJpar);
  perpVisc.push_back(artViscParMomDens);
  perpVisc.push_back(artViscParVort);

  if (viscosityGuard) {
    for (std::vector<BoutReal>::iterator it = perpVisc.begin();
         it != perpVisc.end(); ++it) {
      if (*it > 100.0 * eta0INorm) {
        throw BoutException("One of the parallel viscosities "
                            "is 100.0 times larger than eta0I");
      }
    }
  }

  SAVE_ONCE2(artViscParLnN, artViscParJpar);
  SAVE_ONCE2(artViscParMomDens, artViscParVort);
  SAVE_ONCE2(artViscPerpLnN, artViscPerpJPar);
  SAVE_ONCE2(artViscPerpMomDens, artViscPerpVort);
  SAVE_ONCE(artHyperAzVort);
  // ************************************************************************
}
// ############################################################################

// Timestep initialization
// ############################################################################
void CelmaWBA::timestepInitialization() {
  TRACE("Halt in CelmaWBA::timestepInitialization");

  // Manually specifying rho inner ghost points
  // ************************************************************************
  ownBC.innerRhoCylinder(lnN);
  ownBC.innerRhoCylinder(jPar);
  ownBC.innerRhoCylinder(momDensPar);
  ownBC.innerRhoCylinder(phi);  // Used to set BC in lapalace inversion
  ownBC.innerRhoCylinder(vort); // Later taken derivative of
  // The inner boundaries of phi is set in the inversion procedure
  // ************************************************************************

  // Preparations
  // ************************************************************************
  /* NOTE: Preparations
   * 1. Although we are not using the Naulin solver, we still calculate
   *    gradPerpLnN here, which requires lnN first to be communicated
   * 2. So does vort, as it is used in the laplace inversion
   * 3. n is need several places, so it is calculated here
   * 4. uIPar and uEPar are also needed, and thus calculated here
   *
   * Note also:
   * 1. No further derivatives is taken for GradPerp(lnN), so this is not
   *    needed to be communicated.
   */
  mesh->communicate(lnN, vort);
  gradPerpLnN = ownOp->Grad_perp(lnN);
  n = exp(lnN);
  uIPar = momDensPar / n;      // momDensPar = n*uIPar
  uEPar = -(jPar) / n + uIPar; // j = n(uIPar - uEPar)
  // ************************************************************************

  // Laplace inversion
  // ************************************************************************
  /* NOTE: Solve for phi
   * 1. phi and vort is output
   * 2. The solver takes care of perpendicular boundaries
   * 3. Filtering of higher modes is done internally with the filter flag in
   *    BOUT.inp
   */
  phi = phiSolver->solve(vort);
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
  n.applyBoundary();
  uIPar.applyBoundary();
  // Use the sheath boundary condition with constant Te for jPar and
  // uEPar at SE. Note that we have parallel uEPar derivatives
  // Here we have phiRef = Lambda
  ownBC.jParSheath(jPar, uEPar, uIPar, phi, n, Lambda, Lambda);
  ownBC.uEParSheath(uEPar, phi, Lambda, Lambda);
  // Set the BC for momDensPar
  ownBC.parDensMomSheath(momDensPar, uIPar, n);
  // ************************************************************************

  // Communicate before taking derivatives
  mesh->communicate(comGroup);
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(CelmaWBA);
