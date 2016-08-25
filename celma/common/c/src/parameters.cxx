#ifndef __PARAMETERS_CXX__
#define __PARAMETERS_CXX__

#include "../include/parameters.hxx"

/*!
 * \brief Constructor to Parameters
 *
 * Inputs initialized through list initialization
 *
 *
 * \param[in] radius Plasma radius [m]
 * \param[in] len    The cylinder length [m]
 * \param[in] n0     The density [m^-3]
 * \param[in] Te0    The electron temperature [eV]
 * \param[in] Ti0    The ion temperature [eV]
 * \param[in] B0     The magnetic field [T]
 * \param[in] S      The source [m^-3s^-1]
 */
Parameters::Parameters(BoutReal const &radius,
                       BoutReal const &len,
                       BoutReal const &n0,
                       BoutReal const &Te0,
                       BoutReal const &Ti0,
                       BoutReal const &B0,
                       BoutReal const &S
                       )
: radius_(radius), len_(len), n0_(n0), Te0_(Te0), Ti0_(Ti0), B0_(B0), S_(S),
  separatorLen(49), separator(' '),
  nameWidth(15), numberWidth(15), unitsWidth(16)
{
    TRACE("Parameters::Parameters");

    // Obtained from scipy.constants
    BoutReal const eps0 = 8.854187817620389e-12;
    BoutReal const mu0  = 1.2566370614359173e-06;
    BoutReal const e    = 1.6021766208e-19;
    BoutReal const me   = 9.10938356e-31;
    BoutReal const mp   = 1.672621898e-27;
    int const N = 3; // Degrees of freedom

    BoutReal ne = n0_;
    BoutReal mi = mp;

    // Recalculate the temperatures to joules
    Te0J = Te0_*e;
    Ti0J = Ti0_*e;

    // Thermal quantities
    /* Pecseli:
     * Low frequency waves and turbulence in magnitized plasmas (draft)
     * page 156
     */
    vA   = pow(pow(B0_, 2.0)/(ne*mi*mu0), 0.5);
    cS   = pow((Te0J+((N+2.0)/N)*Ti0J)/mi, 0.5);
    // https://en.wikipedia.org/wiki/Thermal_velocity
    vThE = pow(2.0*Te0J/me, 0.5);
    vThI = pow(2.0*Ti0J/mp, 0.5);

    // Frequencies
    omCI = e*B0_/mi;
    omCE = e*B0_/me;
    omPE = pow(ne*pow(e, 2.0)/(me*eps0), 0.5);

    // Sizes
    debye    = pow(eps0*Te0J/(ne*pow(e, 2.0)), 0.5);
    eLarmour = me*vThE/(e*B0_);
    iLarmour = mi*vThI/(e*B0_);
    rhoS     = cS/omCI;
    Lx       = radius_/rhoS;
    Ly       = len_/rhoS;

    // Collisions
    /*Friedberg:
     * Plasma Physics and Fusion Energy - equation (9.35)
     */
    coloumbLog = log(12.0*PI*pow(eps0*Te0J, 1.5)/(pow(ne, 0.5)*pow(e, 3.0)));
    /* Friedberg:
     * Plasma Physics and Fusion Energy - equation (9.49)
     * Agrees with Bellan page 489
     */
    nuEI = 1.0/(4.0*PI)*(pow(e, 4.0)*ne/(pow(eps0, 2.0)*pow(me,2.0)))*
           coloumbLog*(1.0/pow(vThE, 3.0));
    /* Friedberg:
     * Plasma Physics and Fusion Energy - equation (9.52)
     */
    nuEE = 1.0/(2.0*PI)*(pow(e, 4.0)*ne/(pow(eps0, 2.0)*pow(me, 2.0)))*
           coloumbLog*(1.0/pow(vThE, 3.0));
    nuII = 1.0/(2.0*PI)*(pow(e, 4.0)*ne/(pow(eps0, 2.0)*pow(mi, 2.0)))*
           coloumbLog*(1.0/pow(vThI, 3.0));
    nuIE = 1.0/(4.0*PI)*(pow(e, 4.0)*ne/(pow(eps0, 2.0)*me*mi))*
           coloumbLog*(1.0/pow(vThE, 3.0));

    // Additional parameters
    nuS  = S/ne;
    beta = ne*(Te0J+((N+2.0)/N)*Ti0J)/(pow(B0_, 2.0)/(2.0*mu0));
    mu   = mp/me;
    nuS      = S/ne;

    // Parallel viscosities
    /* Helander, Sigmar:
     * Collisional transport in magnetized plasmas
     * Page 82
     */
    eta0I = 0.96*n0*Ti0J/(pow(2.0, 0.5)*nuII);
    eta2I = 4.0*(3.0*n0*Ti0J*(pow(2.0, 0.5)*nuII))/(10.0*pow(omCI, 2.0));
    eta4I = 2.0*(n0*Ti0J)/(2.0*omCI);
    eta0E = 0.73*n0*Te0J/nuEI;
    eta2E = 4.0*(0.51*n0*Te0J*nuEI)/(pow(omCE, 2.0));
    eta4E = 2.0*(n0*Te0J)/(2.0*omCE);

    // Normalized parameters
    nuEINorm = nuEI/omCI;
    nuSNorm  = nuS/omCI;
    /* Normalization can be found by looking at parallel momentum
     * equation
     */
    eta0INorm = eta0I/(mi*n0*rhoS*cS);
    eta2INorm = eta2I/(mi*n0*rhoS*cS);
    eta4INorm = eta4I/(mi*n0*rhoS*cS);
    eta0ENorm = eta0E/(mi*n0*rhoS*cS);
    eta2ENorm = eta2E/(mi*n0*rhoS*cS);
    eta4ENorm = eta4E/(mi*n0*rhoS*cS);

    // Print the summary
    printTable();

    // Guard
    if(nuEINorm >= 1.0){
        throw BoutException("Normalized nuEI broke drift approximation");
    }
    if(nuSNorm >= 1.0){
        throw BoutException("Normalized nuSNorm broke drift approximation");
    }
}

/*!
 * A printer which prints the variables in a nice way using SI units
 */
void Parameters::printTable() const
{
    TRACE("Parameters::printTable");

    output << "\n" << std::endl;
    output << "GIVEN INPUT" << std::endl;
    output << std::string(separatorLen, '-') << std::endl;
    printVar("n0"    , n0_     , "m^-3"    );
    printVar("B0"    , B0_     , "T"       );
    printVar("Ti0"   , Ti0_    , "eV"      );
    printVar("Te0"   , Te0_    , "eV"      );
    printVar("S"     , S_      , "m^-3s^-1");
    printVar("radius", radius_ , "m"       );
    printVar("len"   , len_    , "m"       );
    output << std::string(separatorLen, '-') << std::endl;
    output << "CONVERTED UNITS" << std::endl;
    output << std::string(separatorLen, '-') << std::endl;
    printVar("Lx"   , Lx    , "-" );
    printVar("Ly"   , Ly    , "-" );
    printVar("Ly/Lx", Ly/Lx , "-" );
    printVar("Ti0"  , Ti0J  , "J");
    printVar("Te0"  , Te0J  , "J");
    output << std::string(separatorLen, '-') << std::endl;
    output << "CODE INPUT" << std::endl;
    output << std::string(separatorLen, '-') << std::endl;
    printVar("mu"  , mu  , "-"   );
    printVar("nuEI", nuEI, "s^-1");
    printVar("nuS" , nuS , "s^-1");
    printVar("beta", beta, "-"   );
    output << std::string(separatorLen, '-') << std::endl;
    output << "PLOT SPECIFIC" << std::endl;
    output << std::string(separatorLen, '-') << std::endl;
    printVar("rhoS", rhoS, "-"   );
    printVar("omCI", omCI, "s^-1");
    output << std::string(separatorLen, '-') << std::endl;
    output << std::string(separatorLen, '-') << std::endl;
    output << "THERMAL QUANTITIES" << std::endl;
    output << std::string(separatorLen, '-') << std::endl;
    printVar("cS"  , cS  , "m s^-1");
    printVar("vA"  , vA  , "m s^-1");
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
    printVar("debye"   , debye   , "m");
    printVar("iLarmour", iLarmour, "m");
    printVar("eLarmour", eLarmour, "m");
    printVar("rhoS"    , rhoS    , "m");
    output << std::string(separatorLen, '-') << std::endl;
    output << "COLLISIONS AND VISCOSITIES" << std::endl;
    output << std::string(separatorLen, '-') << std::endl;
    printVar("coloumbLog", coloumbLog, "-"          );
    printVar("nuEI"      , nuEI      , "s^-1"       );
    printVar("nuIE"      , nuIE      , "s^-1"       );
    printVar("nuEE"      , nuEE      , "s^-1"       );
    printVar("nuII"      , nuII      , "s^-1"       );
    printVar("eta0I"     , eta0I     , "kg m^-1s^-1");
    printVar("eta2I"     , eta2I     , "kg m^-1s^-1");
    printVar("eta4I"     , eta4I     , "kg m^-1s^-1");
    printVar("eta0E"     , eta0E     , "kg m^-1s^-1");
    printVar("eta2E"     , eta2E     , "kg m^-1s^-1");
    printVar("eta4E"     , eta4E     , "kg m^-1s^-1");
    output << std::string(separatorLen, '-') << std::endl;
    output << "ADDITIONAL PARAMETERS" << std::endl;
    output << std::string(separatorLen, '-') << std::endl;
    printVar("nuS" , nuS , "s^-1");
    printVar("beta", beta, "-"   );
    printVar("mu"  , mu  , "-"   );
    output << std::string(separatorLen, '-') << std::endl;
    output << "NORMALIZED PARAMETERS" << std::endl;
    output << std::string(separatorLen, '-') << std::endl;
    printVar("nuEI/rhoS", nuEINorm , "-");
    printVar("nuS/omCI" , nuSNorm  , "-");
    printVar("eta0INorm", eta0INorm, "-");
    printVar("eta2INorm", eta2INorm, "-");
    printVar("eta4INorm", eta4INorm, "-");
    printVar("eta0ENorm", eta0ENorm, "-");
    printVar("eta2ENorm", eta2ENorm, "-");
    printVar("eta4ENorm", eta4ENorm, "-");
    output << std::string(separatorLen, '-') << std::endl;
    output << "\n" << std::endl;
}

/*!
 * Prints each variable
 *
 * \param[in] name  Name  of the variable
 * \param[in] val   Value of the variable
 * \param[in] units   Units of the variable
 */
void Parameters::printVar(std::string const &name,
                          BoutReal    const &val,
                          std::string const &units)
                          const
{
    TRACE("Parameters::printVars");

    output << "| "
           << std::left
           << std::setw(nameWidth)
           << std::setfill(separator)
           << name
           << std::setw(numberWidth)
           << std::setfill(separator)
           << std::scientific
           << std::setprecision(3)
           << val
           << std::left
           << std::setw(unitsWidth)
           << std::setfill(separator)
           << "[" + units + "]"
           << "|"
           << std::endl;
}

/*!
 * Returns the normalized domain radius
 *
 * \param[out] Lx in [-]
 */
BoutReal Parameters::getLx() const
{
    TRACE("Parameters::getLx");
    return Lx;
}

/*!
 * Returns the normalized domain frequency
 *
 * \param[out] Ly in [-]
 */
BoutReal Parameters::getLy() const
{
    TRACE("Parameters::getLy");
    return Ly;
}

/*!
 * Returns the ion gyration frequency
 *
 * \param[out] nuEI in [s^-1]
 */
BoutReal Parameters::getNuEI() const
{
    TRACE("Parameters::getNuEI");
    return nuEI;
}

/*!
 * Returns the mass ratio
 *
 * \param[out] mu
 */
BoutReal Parameters::getMu() const
{
    TRACE("Parameters::getMu");
    return mu;
}

/*!
 * Returns the particle creation frequency
 *
 * \param[out] nuS in [s^-1]
 */
BoutReal Parameters::getNuS() const
{
    TRACE("Parameters::getNuS");
    return nuS;
}

/*!
 * Returns the plasma beta
 *
 * \param[out] beta
 */
BoutReal Parameters::getBeta() const
{
    TRACE("Parameters::getBeta");
    return beta;
}

/*!
 * Returns the ion gyration frequency
 *
 * \param[out] omCI in [s^-1]
 */
BoutReal Parameters::getOmCI() const
{
    TRACE("Parameters::getOmCI");
    return omCI;
}

/*!
 * Returns the hybrid radius
 *
 * \param[out] rhoS in [m]
 */
BoutReal Parameters::getRhoS() const
{
    TRACE("Parameters::getRhoS");
    return rhoS;
}

#endif
