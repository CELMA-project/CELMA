#ifndef __OWNMONITORS_CXX__
#define __OWNMONITORS_CXX__

#include "../include/ownMonitors.hxx"

/* FIXME: Monitors are not printing traces or BoutException at the moment
 *        (see BOUT-dev #355)
 *        As a quickfix std::cout is used for the errors
 */

/*!
 * \brief Constructor
 *
 * Constructor which sets the private member data
 */
OwnMonitors::OwnMonitors()
 :
    polAvgN_    (0.0),
    polAvgLogN_ (0.0),
    polAvgUEPar_(0.0),
    polAvgUIPar_(0.0)
{
    TRACE("Halt in OwnMonitors::OwnMonitors");

    polAvgGradPerpPhi_ = 0.0;
}

/*!
 * Calculates the poloidal average of n.
 *
 * The result is stored in the member data polAvgN_
 *
 * \param[in] n The density
 */
void OwnMonitors::calcPolAvgN(Field3D const &n)
{
    TRACE("Halt in OwnMonitors::calcPolAvgN");

    polAvgN_ = polAvg_.poloidalAverage(n);
}

/*!
 * Calculates the kinetic energy and the poloidal averaged energy of the system
 *
 *  ## CAUTION!!!
 *
 * \warning calcPolAvgN must be called in advance of this function.
 * \warning This function does not multiply with \f$\frac{m_{\alpha}}{m_i}\f$.
 *
 * ## Derivation
 * \note A tilde in this section denotes a normalized quantity, **NOT** a
 * fluctuation.
 *
 * For a single particle, we have
 *
 * \f{eqnarray}{
 * E_{kin} = \frac{1}{2}m\mathbf{v}^2,
 * \f}
 *
 * which means that for a fluid, we have
 *
 * \f{eqnarray}{
 * E_{kin}=\frac{1}{2}m\int n\mathbf{u}^2 dV
 * \f}
 *
 * Due to gyroviscous cancelation, we have that to first order, only the
 * ExB-drift is carrying particles. This gives
 *
 * \f{eqnarray}{
 * E_{kin,\alpha} =& \frac{1}{2}m_{\alpha}\int
 *                      n\mathbf{u}_E^2
 *                      + n\mathbf{u}_{\alpha,\parallel}^2 dV\\
 *                =& \frac{1}{2}m_{\alpha}\int
 *                      n\left(\frac{-\nabla_\perp\phi
 *                             \times\mathbf{b}}{B}\right)^2
 *                      + n\mathbf{u}_{\alpha,\parallel}^2 dV\\
 *                =& \frac{1}{2}m_{\alpha}\int
 *                      n\left(\left[\frac{-\nabla_\perp\phi}{B}\right]^2
 *                      + [\mathbf{u}_{\alpha,\parallel}]^2\right) dV\\
 *                =& \frac{1}{2}m_{\alpha}\int\int\int
 *                      n\left(\left[\frac{-\nabla_\perp\phi}{B}\right]^2
 *                      + [\mathbf{u}_{\alpha,\parallel}]^2\right)
 *                      J d\rho d\theta dz\\
 *                =& \frac{1}{2}m_i\frac{m_{\alpha}}{m_i}
 *                      n_0c_s^2\rho_s^3
 *                      \int\int\int
 *                      \tilde{n}\tilde{u}_\alpha^2
 *                      \tilde{J} \widetilde{d\rho} d\theta \widetilde{dz}\\
 *                =& m_in_0c_s^2\rho_s^3 \tilde{E}_{kin,\alpha}
 *                =& n_0T_e\rho_s^3 \tilde{E}_{kin,\alpha}
 * \f}
 *
 * where we have used (V.4) in D'Haeseleer, and where \f$\alpha\f$ is the
 * particle species. This function calculates
 * \f$\frac{m_i}{m_{\alpha}}\tilde{E}_{kin,\alpha}\f$
 *
 * ## Poloidal averaging
 * \note A tilde in this section denotes a fluctuation, **NOT** a normalized
 * quantity.
 *
 * We will here use the notation
 *
 * \f{eqnarray}{
 *      \langle f \rangle =& \frac{\int_0^{2\pi} f d\theta}{2\pi}\\
 *      \widetilde{f} =& f - \langle f \rangle
 * \f}
 *
 * this gives
 *
 * \f{eqnarray}{
 *      nu^2 =& (\langle n \rangle + \widetilde{n})
 *              (\langle u \rangle + \widetilde{u})^2\\
 *      nu^2 =& (\langle n \rangle + \widetilde{n})
 *              (\langle u \rangle^2
 *              + 2\widetilde{u}\langle u \rangle
 *              + \widetilde{ u }^2
 *              )\\
 *      nu^2 =& \langle n \rangle \langle u \rangle^2
 *              + 2 \langle n \rangle \widetilde{u}\langle u \rangle
 *              + \langle n \rangle\widetilde{ u }^2
 *              +
 *              \widetilde{n} \langle u \rangle^2
 *              + 2 \widetilde{n} \widetilde{u}\langle u \rangle
 *              + \widetilde{n}\widetilde{ u }^2\\
 *      nu^2 - \langle n \rangle \langle u \rangle^2
 *              =& 2 \langle n \rangle \widetilde{u}\langle u \rangle
 *              + \langle n \rangle\widetilde{ u }^2
 *              +
 *              \widetilde{n} \langle u \rangle^2
 *              + 2 \widetilde{n} \widetilde{u}\langle u \rangle
 *              + \widetilde{n}\widetilde{ u }^2
 * \f}
 *
 * as this quantity is integrated over the volume, it will be integrated
 * poloidally, which means that we are going to take the equvialence of an
 * \f$2\pi\langle f \rangle\f$ operation on the above. As
 * \f$\langle f \rangle\f$ is merely a constant, and as
 * \f$\langle \widetilde{f} \rangle = 0\f$ (as many fluctuations above the mean
 * as below, we get
 *
 * \f{eqnarray}{
 *     \langle nu^2 - \langle n \rangle \langle u \rangle^2 \rangle
 *      =& \langle n \rangle \langle u \rangle^2
 *         + \widetilde{n}\widetilde{ u }^2
 * \f}
 *
 * Notice that this means that
 *
 * \f{eqnarray}{
 *      E_{kin} = \langle E_{kin} \rangle + \widetilde{E}_{kin}
 * \f}
 *
 *
 * ## Input output
 *
 * \param[in] n        The density
 * \param[in] uSquared The parallel velocity of species \f$\alpha\f$
 * \param[in] kinE    Variable where the kinetic energy will be stored\n
 *                    Must contain the following keys:\n
 *                    perpKinEE       - Electron perpendicular kinetic energy\n
 *                    parKinEE        - Electron parallel kinetic energy\n
 *                    totKinEE        - Electron total kinetic energy\n
 *                    perpKinEI       - Ion perpendicular kinetic energy\n
 *                    parKinEI        - Ion parallel kinetic energy\n
 *                    totKinEI        - Ion total kinetic energy\n
 *                    polAvgPerpKinEE - Poloidally averaged electron
 *                                      perpendicular kinetic energy\n
 *                    polAvgParKinEE  - Poloidally averaged electron parallel
 *                                      kinetic energy\n
 *                    polAvgPerpKinEE - Poloidally averaged electron total
 *                                      kinetic energy\n
 *                    polAvgPerpKinEI - Poloidally averaged ion perpendicular
 *                                      kinetic energy\n
 *                    polAvgParKinEI  - Poloidally averaged ion parallel kinetic
 *                                      energy\n
 *                    polAvgPerpKinEI - Poloidally averaged ion total kinetic
 *                                      energy
 *
 * \param[out] kinE   Variable where the kinetic energy is stored
 */
void OwnMonitors::kinEnergy(Field3D  const &n                    ,
                            Vector3D const &gradPerpPhi          ,
                            Field3D  const &uEPar                ,
                            Field3D  const &uIPar                ,
                            std::map<std::string, BoutReal> *kinE)
{
    TRACE("Halt in OwnMonitors::kinEnergy");

    // Electron energy
    if((*kinE).count("perpKinEE")){
        (*kinE)["perpKinEE"] =
            volInt_.volumeIntegral(0.5*n*gradPerpPhi*gradPerpPhi);
    }
    else{
        std::cout << "'perpKinEE' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'perpKinEE' was not a key in the input 'kinE'");
    }
    if((*kinE).count("parKinEE")){
        (*kinE)["parKinEE"] =
            volInt_.volumeIntegral(0.5*n*SQ(uEPar));
    }
    else{
        std::cout << "'parKinEE' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'parKinEE' was not a key in the input 'kinE'");
    }
    if((*kinE).count("totKinEE" )){
        (*kinE)["toKinEE"] = (*kinE)["perpKinEE"] + (*kinE)["parKinEE"];
    }
    else{
        std::cout << "'totKinEE' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'totKinEE' was not a key in the input 'kinE'");
    }

    // Ion energy
    if((*kinE).count("perpKinEI")){
        (*kinE)["perpKinEI"] = (*kinE)["perpKinEE"];
    }
    else{
        std::cout << "'perpKinEI' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'perpKinEI' was not a key in the input 'kinE'");
    }
    if((*kinE).count("parKinEI")){
        (*kinE)["parKinEI"] =
            volInt_.volumeIntegral(0.5*n*SQ(uIPar));
    }
    else{
        std::cout << "'parKinEI' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'parKinEI' was not a key in the input 'kinE'");
    }
    if((*kinE).count("totKinEI" )){
        (*kinE)["toKinEI"] = (*kinE)["perpKinEI"] + (*kinE)["parKinEI"];
    }
    else{
        std::cout << "'perpKinEI' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'perpKinEI' was not a key in the input 'kinE'");
    }

    // Calculate the poloidal averages
    /* NOTE: <f>^2 neq <f^2>
     *       As f = <f> + \tilde{f}
     *       As a consequence, we must take the mean before squaring gradPerpPhi
     */
    polAvgGradPerpPhi_.x = polAvg_.poloidalAverage(gradPerpPhi.x);
    polAvgGradPerpPhi_.x = polAvg_.poloidalAverage(gradPerpPhi.x);
    polAvgGradPerpPhi_.y = polAvg_.poloidalAverage(gradPerpPhi.y);
    polAvgGradPerpPhi_.z = polAvg_.poloidalAverage(gradPerpPhi.z);
    polAvgUEPar_         = polAvg_.poloidalAverage(uEPar);
    polAvgUIPar_         = polAvg_.poloidalAverage(uIPar);

    // Poloidally averaged electron energy
     if((*kinE).count("polAvgPerpKinEE")){
        (*kinE)["polAvgPerpKinEE"] =
            volInt_.volumeIntegral(0.5*polAvgN_*polAvgGradPerpPhi_*polAvgGradPerpPhi_);
    }
    else{
        std::cout << "'polAvgPerpKinEE' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'polAvgPerpKinEE' was not a key in the input 'kinE'");
    }
    if((*kinE).count("polAvgParKinEE")){
        (*kinE)["polAvgParKinEE"] =
            volInt_.volumeIntegral(0.5*polAvgN_*SQ(polAvgUEPar_));
    }
    else{
        std::cout << "'polAvgParKinEE' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'polAvgParKinEE' was not a key in the input 'kinE'");
    }
    if((*kinE).count("polAvgTotKinEE" )){
        (*kinE)["polAvgTotKinEE"] = (*kinE)["polAvgPerpKinEE"] + (*kinE)["polAvgParKinEE"];
    }
    else{
        std::cout << "'polAvgTotKinEE' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'polAvgTotKinEE' was not a key in the input 'kinE'");
    }

    // Poloidally averaged ion energy
    if((*kinE).count("polAvgPerpKinEI")){
        (*kinE)["polAvgPerpKinEI"] = (*kinE)["polAvgPerpKinEE"];
    }
    else{
        std::cout << "'polAvgPerpKinEI' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'polAvgPerpKinEI' was not a key in the input 'kinE'");
    }
    if((*kinE).count("polAvgParKinEI")){
        (*kinE)["polAvgParKinEI"] =
            volInt_.volumeIntegral(0.5*polAvgN_*SQ(polAvgUIPar_));
    }
    else{
        std::cout << "'polAvgParKinEI' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'polAvgParKinEI' was not a key in the input 'kinE'");
    }
    if((*kinE).count("polAvgTotKinEI" )){
        (*kinE)["polTotPerpKinEI"] = (*kinE)["polAvgPerpKinEI"] + (*kinE)["polAvgParKinEI"];
    }
    else{
        std::cout << "'polAvgPerpKinEI' was not a key in the input 'kinE'" << std::endl;
        throw BoutException("'polAvgPerpKinEI' was not a key in the input 'kinE'");
    }
}

/*!
 * Calculates the electron potential energy and the poloidal averaged potential
 * energy of the system
 *
 * ## CAUTION !!
 *
 * \warning calcPolAvgN must be called in advance of this function.
 * \warning This function does not multiply with \f$T_e\f$.
 *
 * ## Derivation
 * A proper derivation can be done where one use the variational
 * principle of the Langrangian of the system to find the energy of the
 * system.
 *
 * One can also derive the potential energy in a less formal way by
 * looking at the transfer terms of the equations, and find which energy
 * terms which are not the kinetic energy. If one neglects the small
 * energy from the fields, one is left with. When this is done, one find
 * that (since the system is isotherm)
 *
 * \f{eqnarray}{
 * E_{pot} =& \int nT_{e, 0}\log(n/n_0) dV
 * \f}
 *
 * \note As \f$ T_{i, 0} = 0 \f$, there is no ion potential energy in the
 * system
 *
 * ## Poloidal averaging
 * For derivation
 *
 * \sa kinEnergy
 *
 * (using \f$ n\log(n) \f$ for \f$ u^2 \f$, and \f$ 1 \f$ for \f$ n \f$).
 *
 * As we are dealing with a logarithm here, care must be taken in what we
 * define as the fluctuation and the mean. We still have that
 *
 * \f{eqnarray}{
 *      E_{pot} = \langle E_{pot} \rangle + \widetilde{E}_{pot}
 * \f}
 *
 * if we define (as \f$ T_{e,0} \f$ is constant in this isothermal approach)
 *
 * \f{eqnarray}{
 *      \langle E_{pot} \rangle
 *      = T_{e, 0}\int \langle n \rangle\langle \log(n/n_0)\rangle dV
 * \f}
 *
 * \note \f$ \langle \log(n/n_0) \rangle \neq \log(\langle n \rangle/n_0) \f$
 *
 * \param[in] n       The density
 *
 * \param[in] potE    Variable where the potential will be stored\n
 *                    Must contain the following keys:\n
 *                    potEE       - Electron potential energy\n
 *                    polAvgPotEE - Poloidally averaged electron potential
 *                                  energy
 *
 * \param[out] potE   Variable where the potential electron energy is stored
 */
void OwnMonitors::potEnergy(Field3D const &n, std::map<std::string, BoutReal> *potE)
{
    TRACE("Halt in OwnMonitors::potEnergy");

    // The full part
    if((*potE).count("potEE")){
        (*potE)["potEE"] =
            volInt_.volumeIntegral(n*log(n));
    }
    else{
        std::cout << "'potEE' was not a key in the input 'potE'" << std::endl;
        throw BoutException("'potEE' was not a key in the input 'potE'");
    }

    polAvgLogN_ = polAvg_.poloidalAverage(log(n));

    // The poloidal average
    if((*potE).count("polAvgPotEE")){
        (*potE)["polAvgPotEE"] =
            volInt_.volumeIntegral(polAvgN_*polAvgLogN_);
    }
    else{
        std::cout << "'polAvgPotEE' was not a key in the input 'potE'" << std::endl;
        throw BoutException("'polAvgPotEE' was not a key in the input 'potE'");
    }
}


/*!
 * Calculates the total particle number and the fluctuation of the total
 * particle number in the system
 *
 * \param[in] n       The density
 * \param[in] totN    Variable where the potential will be stored\n
 *                    Must contain the following keys:\n
 *                    totN       - Total particle number\n
 *                    polAvgTotN - Poloidally averaged of total particle number
 *
 * \param[out] totN  Variable where the potential electron energy is stored
 */
void OwnMonitors::totalN(Field3D const &n, std::map<std::string, BoutReal> *totN)
{
    TRACE("Halt in OwnMonitors::totalN");

    // The full part
    if((*totN).count("totN")){
        (*totN)["totN"] =
            volInt_.volumeIntegral(n*log(n));
    }
    else{
        std::cout << "'totN' was not a key in the input 'totN'" << std::endl;
        throw BoutException("'totN' was not a key in the input 'totN'");
    }

    polAvgLogN_ = polAvg_.poloidalAverage(log(n));

    // The poloidal average
    if((*totN).count("polAvgTotN")){
        (*totN)["polAvgTotN"] =
            volInt_.volumeIntegral(polAvgN_*polAvgLogN_);
    }
    else{
        std::cout << "'polAvgTotN' was not a key in the input 'totN'" << std::endl;
        throw BoutException("'polAvgTotN' was not a key in the input 'totN'");
    }

}
#endif
