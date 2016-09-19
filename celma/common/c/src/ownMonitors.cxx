#ifndef __OWNMONITORS_CXX__
#define __OWNMONITORS_CXX__

#include "../include/ownMonitors.hxx"

/*!
 * Calculates the kinetic energy of the system
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
 *                      \tilde{n}\tilde{u}_\alpha^2
 *                      \tilde{J} \widetilde{d\rho} d\theta \widetilde{dz}\\
 *                =& m_in_0c_s^2\rho_s^3 \tilde{E}_{kin,\alpha}
 * \f}
 *
 * where we have used (V.4) in D'Haeseleer, and where \f$\alpha\f$ is the
 * particle species. This function calculates
 * \f$\frac{m_i}{m_{\alpha}}\tilde{E}_{kin,\alpha}\f$
 *
 * \warning This function does not multiply with \f$\frac{m_{\alpha}}{m_i}\f$.
 *
 * \param[in] n       The density
 * \param[in] phi     The potential
 * \param[in] uPar    The parallel velocity of species \f$\alpha\f$
 * \param[in] kinE    Variable where the kinetic energy will be stored\n
 *                    kinE[0] - The perpendicular kinetic energy\n
 *                    kinE[1] - The parallel kinetic energy\n
 *                    kinE[2] - The total kinetic energy\n
 *
 * \param[out] kinE   Variable where the kinetic energy is stored
 *                    kinE[0] - The perpendicular kinetic energy\n
 *                    kinE[1] - The parallel kinetic energy\n
 *                    kinE[2] - The total kinetic energy\n
 */
void OwnMonitors::kinEnergy(Field3D  const &n          ,
                            Vector3D const &gradPerpPhi,
                            Field3D  const &uPar       ,
                            std::vector<BoutReal> *kinE)
{
    TRACE("Halt in OwnMonitors::kinEnergy");

    // Guard
    if ((*kinE).size() != 3){
        throw BoutException("'kinE' must have length 3");
    }

    // Reset result
    for (std::vector<BoutReal>::iterator it = kinE->begin();
         it != kinE->end();
         ++it)
    {
        *it = 0.0;
    }

    // Calculate the perpendicular kinetic energy
    volInt.volumeIntegral(0.5*n*gradPerpPhi*gradPerpPhi, (*kinE)[0]);
    // Calculate the parallel kinetic energy
    volInt.volumeIntegral(0.5*n*SQ(uPar), (*kinE)[1]);
    // Calculate the total kinetic energy
    (*kinE)[2] = (*kinE)[0] + (*kinE)[1];
}

/*!
 * Calculates the total particle number of the system
 *
 * \param[in] n       The density
 * \param[in] N       Variable where the total particle number is stored
 *
 * \param[out] N      The total particle number
 */
void OwnMonitors::totalN(Field3D  const &n, BoutReal *N)
{
    TRACE("Halt in OwnMonitors::kinEnergy");

    volInt.volumeIntegral(n, *N);
}
#endif
