#ifndef __PARAMETERS_CXX__
#define __PARAMETERS_CXX__

#include "../include/parameters.hxx"

// FIXME: This constructor is long. Should be refactored!
/*!
 * \brief Constructor to Parameters
 *
 * Inputs initialized through list initialization
 *
 * \warning Currently assumes singly ionized atoms
 *
 * \param[in] radius Plasma radius [m]
 * \param[in] length The cylinder length [m]
 * \param[in] n0     The density [m^-3]
 * \param[in] Te0    The electron temperature [eV]
 * \param[in] Ti0    The ion temperature [eV]
 * \param[in] B0     The magnetic field [T]
 * \param[in] S      The source [m^-3s^-1]
 * \param[in] warn   Will throw warnings for exceptions if parameters are
 * outside
 *                   the drift approximation.
 */
Parameters::Parameters(BoutReal const &radius, BoutReal const &length,
                       BoutReal const &n0, BoutReal const &Te0,
                       BoutReal const &Ti0, BoutReal const &B0,
                       BoutReal const &S, BoutReal const &nn,
                       std::string const gas, bool const warningForException)
    : radius_(radius), length_(length), n0_(n0), Te0_(Te0), Ti0_(Ti0), B0_(B0),
      S_(S), nn_(nn), warn_(warningForException), gas_(gas), separatorLen(57),
      separator(' '), nameWidth(23), numberWidth(15), unitsWidth(16),
      precision(4) {
  TRACE("Parameters::Parameters");

  // Obtained from scipy.constants (CODATA2014)
  BoutReal const eps0 = 8.854187817620389e-12;
  BoutReal const mu0 = 1.2566370614359173e-06;
  BoutReal const e = 1.6021766208e-19;
  BoutReal const me = 9.10938356e-31;
  BoutReal const a0 = 5.2917721067e-11;
  BoutReal const u = 1.660539040e-27;
  int const N = 3; // Degrees of freedom

  BoutReal ne = n0_;
  BoutReal ni = n0_;
  BoutReal mi;

  if (gas_ == "H") {
    // Mass of proton
    mi = 1.672621898e-27;
  } else if (gas_ == "Ar") {
    mi = 39.948 * u;
    // Guard
    bool allPass = true;
    std::ostringstream stream;
    if (Te0_ < 0.1 || Te0_ > 10) {
      stream << "Cross section for electron-neutral collision only "
             << "valid in range 0.1-10eV.\n";
      allPass = false;
    }
    if (Ti0_ < 0.1 || Ti0_ > 10) {
      stream << "Cross section for ion-neutral collision only "
             << "valid in range 0.1-10eV.\n";
      allPass = false;
    }
    // Throw error or warning
    if (!allPass) {
      std::string str = stream.str();
      // Cast the stream to a const char in order to use it in BoutException
      const char *message = str.c_str();
      std::string newlines = "\n\n\n";
      std::string exclamation(80, '!');
      if (warn_) {
        output << newlines << exclamation << "\nWARNING: " << message
               << exclamation << newlines << std::endl;
      } else {
        throw BoutException(message);
      }
    }
  } else {
    std::ostringstream stream;
    stream << "The current gases are implemented:"
           << "\nH\nAr\n"
           << "You specified " << gas_ << std::endl;
    std::string str = stream.str();
    // Cast the stream to a const char in order to use it in BoutException
    const char *message = str.c_str();
    throw BoutException(message);
  }

  // Recalculate the temperatures to joules
  Te0J = Te0_ * e;
  Ti0J = Ti0_ * e;

  // Thermal quantities
  /* Pecseli:
   * Low frequency waves and turbulence in magnitized plasmas (draft)
   * page 156
   */
  vA = pow(pow(B0_, 2.0) / (ne * mi * mu0), 0.5);
  cS = pow((Te0J + ((N + 2.0) / N) * Ti0J) / mi, 0.5);
  // Goldston page 12 equation (1.24)
  vThE = pow(Te0J / me, 0.5);
  vThI = pow(Ti0J / mi, 0.5);

  // Frequencies
  omCI = e * B0_ / mi;
  omCE = e * B0_ / me;
  omPE = pow(ne * pow(e, 2.0) / (me * eps0), 0.5);

  // Sizes
  debye = pow(eps0 * Te0J / (ne * pow(e, 2.0)), 0.5);
  eLarmour = me * vThE / (e * B0_);
  iLarmour = mi * vThI / (e * B0_);
  rhoS = cS / omCI;
  Lx = radius_ / rhoS;
  Ly = length_ / rhoS;

  // Collisions
  /* Goldston:
   * Below equation 11.17 page 169
   *
   * Alternatively
   * Friedberg:
   * Plasma Physics and Fusion Energy - equation (9.35)
   */
  coloumbLog = log(12.0 * PI * ne * pow(debye, 3.0));
  /* Braginskii, Goldston and Helander:
   * Maxwellian averaged
   * Braginskii page 215 equation 2.5e)
   * Goldston page 172 equation (11.22)
   * Helander page 5, equation (1.4)
   *
   * NOTE: This differs by approximately 30% from what presented in
   *      Friedberg - Plasma Physics and Fusion Energy - equation (9.49)
   *      which agrees with
   *      Bellan page 489
   */
  nuEI = (pow(2.0, 0.5) * ni * pow(e, 4.0) * coloumbLog) /
         (12.0 * pow(PI, 1.5) * pow(eps0, 2.0) * pow(me, 0.5) * pow(Te0J, 1.5));
  /* Braginskii, Goldston and Helander:
   * Maxwellian averaged
   * Braginskii page 215 equation 2.5i)
   * Goldston page 173 equation (11.24)
   * Helander page 5, equation (1.5)
   */
  nuII = (ni * pow(e, 4.0) * coloumbLog) /
         (12.0 * pow(PI, 1.5) * pow(eps0, 2.0) * pow(mi, 0.5) * pow(Ti0J, 1.5));

  if (gas_ == "H") {
    output << "Im H" << std::endl;
    /* Own calculations:
     * See thesis
     */
    nuEN = (8.0 * pow(2.0, 0.5) / 3.0) * pow(PI, 0.5) * nn * pow(a0, 2.0) *
           pow(Te0J / me, 0.5);
    nuIN = (8.0 / 3.0) * pow(PI, 0.5) * nn * pow(a0, 2.0) * pow(Ti0J / mi, 0.5);
  } else if (gas_ == "Ar") {
    /* Own calculations:
     * See thesis and calculations in additionalPyScripts
     */
    // NOTE: Some quantities are given in eV
    nuEN = (nn / 2.5e19) *
           (33640.349990 * pow(Te0_, 0) - 33174.059200 * pow(Te0_, 1) +
            642273.100111 * pow(Te0_, 2) - 188328.743082 * pow(Te0_, 3) +
            25742.288823 * pow(Te0_, 4) - 1784.118597 * pow(Te0_, 5) +
            50.336945 * pow(Te0_, 6));

    if (Ti0_ > 0.1) {
      nuIN = (8.0 * pow(2.0, 0.5) / 3.0) * (nn / pow(PI, 0.5)) *
             pow(Ti0J / mi, 0.5) * 9.6e-19 *
             pow((1.0 + 0.14 * log(1.0 / Ti0_)), 2.0);
    } else {
      nuIN = 0.0;
    }
  }

  /* Friedberg:
   * Plasma Physics and Fusion Energy - equation (9.52)
   * NOTE: These are probably not Maxwellian averaged!
   */
  nuEE = 1.0 / (2.0 * PI) *
         (pow(e, 4.0) * ne / (pow(eps0, 2.0) * pow(me, 2.0))) * coloumbLog *
         (1.0 / pow(vThE, 3.0));
  nuIE = 1.0 / (4.0 * PI) * (pow(e, 4.0) * ne / (pow(eps0, 2.0) * me * mi)) *
         coloumbLog * (1.0 / pow(vThE, 3.0));

  // Additional parameters
  beta = ne * (Te0J + ((N + 2.0) / N) * Ti0J) / (pow(B0_, 2.0) / (2.0 * mu0));
  mu = mi / me;
  Lambda = log(pow(mu / (2.0 * PI), 0.5));

  // Parallel viscosities
  /* Helander, Sigmar:
   * Collisional transport in magnetized plasmas
   * Page 82
   */
  eta0I = 0.96 * n0 * Ti0J / (pow(2.0, 0.5) * nuII);
  eta2I = 4.0 * (3.0 * n0 * Ti0J * (pow(2.0, 0.5) * nuII)) /
          (10.0 * pow(omCI, 2.0));
  eta4I = 2.0 * (n0 * Ti0J) / (2.0 * omCI);
  eta0E = 0.73 * n0 * Te0J / nuEI;
  eta2E = 4.0 * (0.51 * n0 * Te0J * nuEI) / (pow(omCE, 2.0));
  eta4E = 2.0 * (n0 * Te0J) / (2.0 * omCE);

  // Normalized parameters
  nuEINorm = nuEI / omCI;
  SNorm = S / (n0_ * omCI);
  nuENNorm = nuEN / omCI;
  nuINNorm = nuIN / omCI;
  /* Normalization can be found by looking at parallel momentum
   * equation
   */
  eta0INorm = eta0I / (mi * n0 * rhoS * cS);
  eta2INorm = eta2I / (mi * n0 * rhoS * cS);
  eta4INorm = eta4I / (mi * n0 * rhoS * cS);
  eta0ENorm = eta0E / (mi * n0 * rhoS * cS);
  eta2ENorm = eta2E / (mi * n0 * rhoS * cS);
  eta4ENorm = eta4E / (mi * n0 * rhoS * cS);

  // Print the summary
  printTable();

  // Guard
  bool allPass = true;
  std::ostringstream stream;
  if (coloumbLog <= 1.0) {
    // Huba, J.D. - NRL PLASMA FORMULARY 2013
    stream << "Coloumb Logarithm under 1.0. Theory fails.\n";
    allPass = false;
  }

  // Throw error or warning
  if (!allPass) {
    std::string str = stream.str();
    // Cast the stream to a const char in order to use it in BoutException
    const char *message = str.c_str();
    std::string newlines = "\n\n\n";
    std::string exclamation(80, '!');
    if (warn_) {
      output << newlines << exclamation << "\nWARNING: " << message
             << exclamation << newlines << std::endl;
    } else {
      throw BoutException(message);
    }
  }
}

/*!
 * A printer which prints the variables in a nice way using SI units
 */
void Parameters::printTable() const {
  TRACE("Parameters::printTable");

  output << "\n" << std::endl;
  output << "GIVEN INPUT" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("n0", n0_, "m^-3");
  printVar("B0", B0_, "T");
  printVar("Ti0", Ti0_, "eV");
  printVar("Te0", Te0_, "eV");
  printVar("S", S_, "m^-3s^-1");
  printVar("radius", radius_, "m");
  printVar("length", length_, "m");
  printVar("nn", nn_, "m^-3");
  printVar("gas", gas_, "-");
  output << std::string(separatorLen, '-') << std::endl;
  output << "CONVERTED UNITS" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("Ti0", Ti0J, "J");
  printVar("Te0", Te0J, "J");
  printVar("Ly/Lx", Ly / Lx, "-");
  output << std::string(separatorLen, '-') << std::endl;
  output << "CODE INPUT" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("Lx", Lx, "-");
  printVar("Ly", Ly, "-");
  printVar("mu", mu, "-");
  printVar("Lambda", Lambda, "-");
  printVar("beta", beta, "-");
  printVar("nuEI/omCI", nuEINorm, "-");
  printVar("S/n*omCI ", SNorm, "-");
  printVar("nuEN/omCI", nuENNorm, "-");
  printVar("nuIN/omCI", nuINNorm, "-");
  printVar("eta0I/(mi*n0*rhoS*cS)", eta0INorm, "-");
  printVar("eta0E/(mi*n0*rhoS*cS)", eta0ENorm, "-");
  output << std::string(separatorLen, '-') << std::endl;
  output << "PLOT SPECIFIC" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("rhoS", rhoS, "m");
  printVar("omCI", omCI, "s^-1");
  output << std::string(separatorLen, '-') << std::endl;
  output << std::endl;
  output << std::string(separatorLen, '/') << std::endl;
  output << std::endl;
  output << "THERMAL QUANTITIES" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("cS", cS, "m s^-1");
  printVar("vA", vA, "m s^-1");
  printVar("vThE", vThE, "m s^-1");
  printVar("vThI", vThI, "m s^-1");
  output << std::string(separatorLen, '-') << std::endl;
  output << "FREQUENCIES" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("omCI", omCI, "s^-1");
  printVar("omCE", omCE, "s^-1");
  printVar("omPE", omPE, "s^-1");
  output << std::string(separatorLen, '-') << std::endl;
  output << "SIZES" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("debye", debye, "m");
  printVar("iLarmour", iLarmour, "m");
  printVar("eLarmour", eLarmour, "m");
  printVar("rhoS", rhoS, "m");
  output << std::string(separatorLen, '-') << std::endl;
  output << "COLLISIONS AND VISCOSITIES" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("coloumbLog", coloumbLog, "-");
  printVar("nuEI", nuEI, "s^-1");
  printVar("nuII", nuII, "s^-1");
  printVar("nuIEApprox", nuIE, "s^-1");
  printVar("nuEEApprox", nuEE, "s^-1");
  printVar("nuEN", nuEN, "s^-1");
  printVar("nuIN", nuIN, "s^-1");
  printVar("eta0I", eta0I, "kg m^-1s^-1");
  printVar("eta2I", eta2I, "kg m^-1s^-1");
  printVar("eta4I", eta4I, "kg m^-1s^-1");
  printVar("eta0E", eta0E, "kg m^-1s^-1");
  printVar("eta2E", eta2E, "kg m^-1s^-1");
  printVar("eta4E", eta4E, "kg m^-1s^-1");
  output << std::string(separatorLen, '-') << std::endl;
  output << "ADDITIONAL PARAMETERS" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("beta", beta, "-");
  printVar("mu", mu, "-");
  printVar("Lambda", Lambda, "-");
  output << std::string(separatorLen, '-') << std::endl;
  output << "NORMALIZED PARAMETERS" << std::endl;
  output << std::string(separatorLen, '-') << std::endl;
  printVar("nuEI/omCI", nuEINorm, "-");
  printVar("nuEN/omCI", nuENNorm, "-");
  printVar("nuIN/omCI", nuINNorm, "-");
  printVar("S/n0*omCI", SNorm, "-");
  printVar("eta0I/(mi*n0*rhoS*cS)", eta0INorm, "-");
  printVar("eta2I/(mi*n0*rhoS*cS)", eta2INorm, "-");
  printVar("eta4I/(mi*n0*rhoS*cS)", eta4INorm, "-");
  printVar("eta0E/(mi*n0*rhoS*cS)", eta0ENorm, "-");
  printVar("eta2E/(mi*n0*rhoS*cS)", eta2ENorm, "-");
  printVar("eta4E/(mi*n0*rhoS*cS)", eta4ENorm, "-");
  output << std::string(separatorLen, '-') << std::endl;
  output << "\n" << std::endl;
}

/*!
 * Prints each variable
 *
 * \param[in] name      Name  of the variable
 * \param[in] val       Value of the variable
 * \param[in] units     Units of the variable
 */
void Parameters::printVar(std::string const &name, BoutReal const &val,
                          std::string const &units) const {
  TRACE("Parameters::printVars");

  output << "| " << std::left << std::setw(nameWidth) << std::setfill(separator)
         << name << std::setw(numberWidth) << std::setfill(separator)
         << std::scientific << std::setprecision(precision) << val << std::left
         << std::setw(unitsWidth) << std::setfill(separator)
         << "[" + units + "]"
         << "|" << std::endl;
}

/*!
 * Prints each variable
 *
 * \param[in] name      Name  of the variable
 * \param[in] val       Value of the variable
 * \param[in] units     Units of the variable
 */
void Parameters::printVar(std::string const &name, std::string const &val,
                          std::string const &units) const {
  TRACE("Parameters::printVars");

  output << "| " << std::left << std::setw(nameWidth) << std::setfill(separator)
         << name << std::setw(numberWidth) << std::setfill(separator)
         << std::scientific << std::setprecision(precision) << val << std::left
         << std::setw(unitsWidth) << std::setfill(separator)
         << "[" + units + "]"
         << "|" << std::endl;
}

/*!
 * Returns the normalized domain radius
 *
 * \param[out] Lx The normalized version of the radius
 */
BoutReal Parameters::getLx() const {
  TRACE("Parameters::getLx");
  return Lx;
}

/*!
 * Returns the normalized domain length
 *
 * \param[out] Ly The normalized version of the length
 */
BoutReal Parameters::getLy() const {
  TRACE("Parameters::getLy");
  return Ly;
}

/*!
 * Returns the electron ion collision frequency
 *
 * \param[out] nuEI The normalized version of nuEI
 */
BoutReal Parameters::getNuEINorm() const {
  TRACE("Parameters::getNuEI");
  return nuEINorm;
}

/*!
 * Returns the electron neutral collision frequency
 *
 * \param[out] nuEN The normalized version of nuIN
 */
BoutReal Parameters::getNuENNorm() const {
  TRACE("Parameters::getNuEN");
  return nuENNorm;
}

/*!
 * Returns the ion neutral collision frequency
 *
 * \param[out] nuIN The normalized version of nuIN
 */
BoutReal Parameters::getNuINNorm() const {

  TRACE("Parameters::getNuIN");
  return nuINNorm;
}

/*!
 * Returns the normalized particle creation rate
 *
 * \param[out] SNorm The normalized S a.k.a normalized /f$/nu_S/f$
 */
BoutReal Parameters::getSNorm() const {
  TRACE("Parameters::getMu");
  return SNorm;
}

/*!
 * Returns the mass ratio
 *
 * \param[out] mu /f$\frac{m_i}{m_e}/f$
 */
BoutReal Parameters::getMu() const {
  TRACE("Parameters::getMu");
  return mu;
}

/*!
 * Returns the Lambda
 *
 * \param[out] Lambda \f$\ln\left(\sqrt{\frac{\mu}{2\pi}}\right)\f$
 */
BoutReal Parameters::getLambda() const {
  TRACE("Parameters::getLambda");
  return Lambda;
}

/*!
 * Returns the plasma beta
 *
 * \param[out] beta The plasma beta
 */
BoutReal Parameters::getBeta() const {
  TRACE("Parameters::getBeta");
  return beta;
}

/*!
 * Returns the ion gyration frequency
 *
 * \param[out] omCI The ion gyration frequency in [s^-1]
 */
BoutReal Parameters::getOmCI() const {
  TRACE("Parameters::getOmCI");
  return omCI;
}

/*!
 * Returns the hybrid radius
 *
 * \param[out] rhoS The hybrid radius in [m]
 */
BoutReal Parameters::getRhoS() const {
  TRACE("Parameters::getRhoS");
  return rhoS;
}

/*!
 * Returns normalized version of /f$\eta_0^i/f$
 *
 * \param[out] eta0INorm The normalized viscosity coefficient
 */
BoutReal Parameters::getEta0INorm() const {

  TRACE("Parameters::getEta0I");
  return eta0INorm;
}

/*!
 * Returns normalized version of /f$\eta_0^e/f$
 *
 * \param[out] eta0ENorm The normalized viscosity coefficient
 */
BoutReal Parameters::getEta0ENorm() const {

  TRACE("Parameters::getEta0E");
  return eta0ENorm;
}

#endif
