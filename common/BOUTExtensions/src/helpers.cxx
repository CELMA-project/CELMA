#ifndef __HELPERS_CXX__
#define __HELPERS_CXX__

#include "../include/helpers.hxx"

/*!
 * Returns the poloidal average of a field
 *
 * Returns the numerical equivalence of
 *
 * \f{eqnarray}{
 *      \frac{\int_0^{2\pi} f \rho d\theta}{\int_0^{2\pi} \rho d\theta}
 *      = \frac{\int_0^{2\pi} f d\theta}{\int_0^{2\pi} d\theta}
 *      = \frac{\int_0^{2\pi} f d\theta}{2\pi}
 * \f}
 *
 * by using the rectangular rule for integration
 *
 * \param[in] f     The field to take the average of
 *
 * \returns result  The poloidal average of the field
 */
Field3D const PolAvg::poloidalAverage(const Field3D &f) {
  TRACE("Halt in PolAvg::polAvg");

  Field3D result = 0.0;
  BoutReal avg;

  for (xInd = mesh->xstart; xInd <= mesh->xend; xInd++) {
    for (yInd = mesh->ystart; yInd <= mesh->yend; yInd++) {
      // Find the average
      avg = 0.0;
      for (zInd = 0; zInd < mesh->LocalNz; zInd++) {
        avg += f(xInd, yInd, zInd);
      }
      avg /= mesh->LocalNz;

      // Fill the poloidal points with the average
      for (zInd = 0; zInd < mesh->LocalNz; zInd++) {
        result(xInd, yInd, zInd) = avg;
      }
    }
  }

  return result;
}

/*!
 * Returns the volume integrated of a field using the rectangle rule in 3D, and
 * assuming that the point is in the center of a volume cell.
 *
 * We have that the volume integral is given by (see D'Haeseleer chapter 2.5.4
 * c) )
 *
 * \f{eqnarray}{
 * \int\int\int f d^3R
 * = \int\int\int f J du^1 du^2 du^3
 * \simeq \sum_{u^3}\sum_{u^2}\sum_{u^1} f J du^1 du^2 du^3
 * \f}
 *
 * Where \f$d u^i\f$ denotes the coordinate curves
 *
 * \param[in] f      The field to integrate over
 * \param[in] result Variable to store the result
 *
 * \return result The volume integral of field f
 */
BoutReal VolumeIntegral::volumeIntegral(Field3D const &f) {
  TRACE("Halt in VolumeIntegral::volumeIntegral");

  // Make a local variable (which will be safeCollected by MPI_Allreduce)
  BoutReal result = 0.0;
  BoutReal localResult = 0.0;

  Coordinates *coord = mesh->getCoordinates();

  // We loop over the processor domain
  for (xInd = mesh->xstart; xInd <= mesh->xend; xInd++) {
    for (yInd = mesh->ystart; yInd <= mesh->yend; yInd++) {
      for (zInd = 0; zInd < mesh->LocalNz; zInd++) {
        localResult +=
            f(xInd, yInd, zInd) * coord->J(xInd, yInd) *
            coord->dx(xInd, yInd) *
            coord->dy(xInd, yInd) * coord->dz;
      }
    }
  }

  // Sum the data from the processors, and store it in result
  MPI_Allreduce(&localResult, &result, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get());

  return result;
}

#endif
