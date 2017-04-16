#ifndef __OWNMONITORS_CXX__
#define __OWNMONITORS_CXX__

#include "../include/ownMonitors.hxx"

/* FIXME: Monitors are not printing traces or BoutException at the moment
 *        (see BOUT-dev #355)
 *        As a quickfix output is used for the errors
 */

/*!
 * Calculates the kinetic energy and the poloidal averaged energy of the system
 *
 *  ## CAUTION!!!
 * \warning This function does not multiply with \f$\frac{m_{\alpha}}{m_i}\f$.
 *
 * ## Derivation
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
 *                      \breve{n}\breve{u}_\alpha^2
 *                      \breve{J} \widetilde{d\rho} d\theta \widetilde{dz}\\
 *                =& m_in_0c_s^2\rho_s^3 \breve{E}_{kin,\alpha}
 *                =& n_0T_e\rho_s^3 \breve{E}_{kin,\alpha}
 * \f}
 *
 * where we have used (V.4) in D'Haeseleer, and where \f$\alpha\f$ is the
 * particle species. This function calculates
 * \f$\frac{m_i}{m_{\alpha}}\breve{E}_{kin,\alpha}\f$
 *
 * ## Input output
 *
 * \param[in] n        The density
 * \param[in] uSquared The parallel velocity of species \f$\alpha\f$
 * \param[in] kinE    Variable where the kinetic energy will be stored\n
 *                    Must contain the following keys:\n
 *                    perpKinEE       - Electron perpendicular kinetic energy\n
 *                    parKinEE        - Electron parallel kinetic energy\n
 *                    perpKinEI       - Ion perpendicular kinetic energy\n
 *                    parKinEI        - Ion parallel kinetic energy\n
 *
 * \param[out] kinE   Variable where the kinetic energy is stored
 */
void OwnMonitors::kinEnergy(Field3D const &n, Vector3D const &gradPerpPhi,
                            Field3D const &uEPar, Field3D const &uIPar,
                            std::map<std::string, BoutReal> *kinE) {
  TRACE("Halt in OwnMonitors::kinEnergy");

  // Electron energy
  if ((*kinE).count("perpKinEE")) {
    (*kinE)["perpKinEE"] =
        volInt_.volumeIntegral(0.5 * n * gradPerpPhi * gradPerpPhi);
  } else {
    output << "'perpKinEE' was not a key in the input 'kinE'" << std::endl;
    throw BoutException("'perpKinEE' was not a key in the input 'kinE'");
  }
  if ((*kinE).count("parKinEE")) {
    (*kinE)["parKinEE"] = volInt_.volumeIntegral(0.5 * n * SQ(uEPar));
  } else {
    output << "'parKinEE' was not a key in the input 'kinE'" << std::endl;
    throw BoutException("'parKinEE' was not a key in the input 'kinE'");
  }

  // Ion energy
  if ((*kinE).count("perpKinEI")) {
    (*kinE)["perpKinEI"] = (*kinE)["perpKinEE"];
  } else {
    output << "'perpKinEI' was not a key in the input 'kinE'" << std::endl;
    throw BoutException("'perpKinEI' was not a key in the input 'kinE'");
  }
  if ((*kinE).count("parKinEI")) {
    (*kinE)["parKinEI"] = volInt_.volumeIntegral(0.5 * n * SQ(uIPar));
  } else {
    output << "'parKinEI' was not a key in the input 'kinE'" << std::endl;
    throw BoutException("'parKinEI' was not a key in the input 'kinE'");
  }
}

/*!
 * Calculates the total particle number.
 *
 * This can be used to calculate the potential energy, as the volumetric
 * potential energy is given by \f$nT\f$.
 *
 * \param[in] n    The density
 * \param[in] particleNumber    Variable where the potential will be stored\n
 *                              Must contain the following key:\n
 *                              particleNumber       - Total particle number\n
 *
 * \param[out] particleNumber  Variable where the potential electron energy is
 *                             stored
 */
void OwnMonitors::numberOfParticles(
    Field3D const &n, std::map<std::string, BoutReal> *particleNumber) {
  TRACE("Halt in OwnMonitors::numberOfParticles");

  // The full part
  if ((*particleNumber).count("particleNumber")) {
    (*particleNumber)["particleNumber"] = volInt_.volumeIntegral(n);
  } else {
    output << "'particleNumber' was not a key in the input 'particleNumber'"
           << std::endl;
    throw BoutException(
        "'particleNumber' was not a key in the input 'particleNumber'");
  }
}
#endif
