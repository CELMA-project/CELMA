#ifndef __HELPERS_CXX__
#define __HELPERS_CXX__

#include "../include/helpers.hxx"

/*!
 * Returns the poloidal average of a field
 *
 * Returns the numerical equivalence of
 *
 * \f{eqnarray}{
 *      \frac{\int_0^{2\pi} f d\theta}{\int_0^{2\pi} d\theta}
 *      = \frac{\int_0^{2\pi} f d\theta}{2\pi}
 * \f}
 *
 * by using the rectangular rule for integration
 *
 * \param[in] f     The field to take the average of
 *
 * \returns result  The poloidal average of the field
 */
Field3D const PolAvg::poloidalAverage(const Field3D &f)
{
    TRACE("Halt in PolAvg::polAvg");

    Field3D result = 0.0;
    BoutReal avg;

    for(xInd = mesh->xstart; xInd <= mesh->xend; xInd++){
        for(yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
            // Find the average
            avg = 0.0;
            for(zInd = 0; zInd < mesh->ngz -1; zInd ++){
                avg += f(xInd, yInd, zInd);
            }
            avg /= (mesh->ngz - 1);

            // Fill the poloidal points with the average
            for(zInd = 0; zInd < mesh->ngz -1; zInd ++){
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
BoutReal VolumeIntegral::volumeIntegral(Field3D const &f)
{
    TRACE("Halt in VolumeIntegral::volumeIntegral");

    // Make a local variable (which will be safeCollected by MPI_Allreduce)
    BoutReal result = 0.0;
    BoutReal localResult = 0.0;

    /* NOTE: Addressing "off by one" looping over the local range
     *       We want to loop starting from (and including) mesh->@start to (and
     *       including) the mesh->@end.
     *       mesh->@end is included by using <= insted of <
     */
    /* NOTE: Addressing "off by one" looping over the z indices
     *       The z plane includes one extra plane which is unused. To account
     *       for this, we subtract by one in the loop. mesh->ngz is the number
     *       of z planes (counting from 1). That is, it is the same as MZ in
     *       the input file. To account for the counting from 0, we simply use
     *       < instead of <= in the loop.
     */
    // We loop over the processor domain
    for (xInd = mesh->xstart; xInd <= mesh->xend; xInd ++){
        for (yInd = mesh->ystart; yInd <= mesh->yend; yInd ++){
            for (zInd = 0; zInd < mesh->ngz - 1; zInd ++){
                localResult +=
                    f       (xInd, yInd, zInd) *
                    mesh->J (xInd, yInd)       *
                    mesh->dx(xInd, yInd)       *
                    mesh->dy(xInd, yInd)       *
                    mesh->dz
                    ;
            }
        }
    }

    // Sum the data from the processors, and store it in result
    MPI_Allreduce(&localResult, &result,
                  1, MPI_DOUBLE, MPI_SUM, BoutComm::get());

    return result;
}

#endif
