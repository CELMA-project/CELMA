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
 *                      + n\mathbf{u}_{\alpha}^2 dV\\
 *                =& \frac{1}{2}m_{\alpha}\int
 *                      n\left(\frac{-\nabla_\perp\times\mathbf{b}}{B}\right)^2
 *                      + n\mathbf{u}_{\alpha}^2 dV\\
 *                =& \frac{1}{2}m_{\alpha}\int
 *                      n\left(\left[\frac{-\nabla_\perp}{B}\right]^2
 *                      + [\mathbf{u}_{\alpha}]^2\right) dV
 * \f}
 *
 * where we have used (V.4) in D'Haeseleer, and where \f$\alpha\f$ is the
 * particle species. Normalizing using
 * \f$\tilde{E}_{kin,\alpha} = \frac{E_{kin,\alpha}}{m_ic_s}\f$ and Bohm
 * normalization yields (dropping tilde)
 *
 * \f{eqnarray}{
 * E_{kin,\alpha} =& \frac{1}{2}\frac{m_{\alpha}}{m_i}\int
 *                  n\left(\left[\nabla_\perp\right]^2
 *                  + [\mathbf{u}_{\alpha}]^2\right) dV
 * \f}
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
                            Vector3D const &GradPerpPhi,
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
    volInt.volumeIntegral(0.5*n*Grad_perp(phi)*Grad_perp(phi), (*kinE)[0]);
    // Calculate the parallel kinetic energy
    volInt.volumeIntegral(0.5*n*SQ(uPar), (*kinE)[1]);
    // Calculate the total kinetic energy
    (*kinE)[2] = (*kinE)[0] + (*kinE)[1];
}

/*!
 * Monitors the outflow rate of a field f
 *
 * Status: Under developement
 *
 * \param[in] f         The field to measure the flow of
 * \param[in] outflowR  Variable to store the outflow
 *
 * \param[out] outflowR  The outflow rate of f
 */
void OwnMonitors::outflowRate(Field3D  const &f, BoutReal *outflowR)
{
    TRACE("OwnMonitors::outflowRate");

    // Guard
    throw BoutException("outflowRate not implemented");
}

#endif
