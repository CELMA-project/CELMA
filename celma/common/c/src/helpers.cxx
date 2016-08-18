#ifndef __HELPERS_CXX__
#define __HELPERS_CXX__

#include "../include/helpers.hxx"

/*!
 * \brief Constructor
 *
 * Constructor which sets the private member data used in the
 * surfaceEdgeIntegral. Note that localResults will be initialized to a vector
 * with 4 entries with value 0.0
 */
Helpers::Helpers() : vDotDSXin   (0.0),
                     vDotDSXout  (0.0),
                     vDotDSYdown (0.0),
                     vDotDSYup   (0.0),
                     localResults(4, 0.0)
{
    TRACE("Halt in Helpers::Helpers");

    // Referring to the components
    dSXin  .covariant = true;
    dSXout .covariant = true;
    dSYdown.covariant = true;
    dSYup  .covariant = true;

    // Vector is pointing in the direction of negative x
    dSXin  .x = - mesh->J*mesh->dy*mesh->dz;
    dSXin  .y = 0.0;
    dSXin  .z = 0.0;

    // Vector is pointing in the direction of positive x
    dSXout .x = mesh->J*mesh->dy*mesh->dz;
    dSXout .y = 0.0;
    dSXout .z = 0.0;

    // Vector is pointing in the direction of negative y
    dSYdown.x = 0.0;
    dSYdown.y = - mesh->J*mesh->dx*mesh->dz;
    dSYdown.z = 0.0;

    // Vector is pointing in the direction of positive y
    dSYup  .x = 0.0;
    dSYup  .y = mesh->J*mesh->dx*mesh->dz;
    dSYup  .z = 0.0;
}

/*!
 * Returns the poloidal average of a field
 *
 * \param[in] f     The field to take the average of
 *
 * \returns result The poloidal average of the field
 */
Field3D const Helpers::polAvg(const Field3D &f)
{
    TRACE("Halt in Helpers::polAvg");

    Field3D result = 0.0;
    BoutReal avg;

    for(xInd = mesh->xstart+1; xInd <= mesh->xend-1; xInd++){
        for(yInd = mesh->ystart; yInd <= mesh->yend; yInd++){
            // Find the average
            avg = 0.0;
            for(zInd = 0; zInd < mesh->ngz -1; zInd ++){
                avg += f(xInd, yInd, zInd);
            }
            avg /= (mesh->ngz - 1);

            // Subtract the average from the field
            for(zInd = 0; zInd < mesh->ngz -1; zInd ++){
                result(xInd, yInd, zInd) = f(xInd, yInd, zInd) - avg ;
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
 * \param[out] result The volume integral of field f
 */
void Helpers::volumeIntegral(Field3D const &f, BoutReal &result)
{
    TRACE("Halt in Helpers::volumeIntegral");

    // Make a local variable (which will be collected by MPI_Allreduce)
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

}

/*!
 * Returns the surface integral of a vector field at the edge of the domain
 * using the rectangle rule in 3D, and that the surface is half between
 * grid-points.
 *
 * We have that the surface integral is given by (see D'Haeseleer chapter 2.5.4
 * b) )
 *
 * \f{eqnarray}{
 *   \int\int \mathbf{v} \cdot d\mathbf{S}(\mathrm{in} u^i = \mathrm{constant})
 * = \int\int \mathbf{v} \cdot d\mathbf{S}(i)
 * = \int\int \mathbf{v} \cdot (\pm J du^j du^k \mathbf{e}^i)
 * \simeq \sum_{u^k}\sum_{u^j} \mathbf{v} \cdot (\pm J du^j du^k \mathbf{e}^i)
 * \f}
 *
 * where \f$d u^i\f$ denotes the coordinate curves, \f$d \mathbf{e}^i \f$
 * denotes a contravariant basis vector, and that the sign must be adopted for
 * the normal of the surface (from equation 2.5.49 b in D'Haeseleer)
 *
 * \param[in] v         The vector field integrated over
 * \param[in] results   std::vector to store the results in
 *                      result[0] - integral over v on the xin surface
 *                      result[1] - integral over v on the xout surface
 *                      result[2] - integral over v on the ydown surface
 *                      result[3] - integral over v on the yup surface
 *
 *
 * \param[out] results  The result of the integration
 *                      result[0] - integral over v on the xin surface
 *                      result[1] - integral over v on the xout surface
 *                      result[2] - integral over v on the ydown surface
 *                      result[3] - integral over v on the yup surface
 *
 * \note There are no integral on the z-surfaces as these are periodic.
 *
 * \warning Only to be used when the surface is located half between the grid
 *          points.
 */
void Helpers::surfaceEdgeIntegral(Vector3D const &v,
                                 std::vector<BoutReal> &results)
{
    TRACE("Halt in Helpers::surfaceEdgeIntegral");

    // Guard
    if (results.size() != 4 ){
        throw BoutException("'results' must have length 4");
    }

    // Reset results
    for (std::vector<BoutReal>::iterator it = results.begin();
         it != results.end();
         ++it)
    {
        *it = 0.0;
        // it - results.begin() is the current index
        localResults[it - results.begin()] = 0.0;
    }

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
    if(mesh->firstX()){
        vDotDSXin = v*dSXin;
        // Loop through the inner boundary points
        xInd = mesh->xstart;
        for (yInd = mesh->ystart; yInd <= mesh->yend; yInd ++){
            for (zInd = 0; zInd < mesh->ngz - 1; zInd ++){
                // Extrapolate to the edge using a 4th order stencil
                localResults[0] += vDotDSXin(xInd, yInd, zInd)
//                localResults[0] +=   35.0*vDotDSXin(xInd  , yInd, zInd)
//                                   - 35.0*vDotDSXin(xInd+1, yInd, zInd)
//                                   + 21.0*vDotDSXin(xInd+2, yInd, zInd)
//                                   -  5.0*vDotDSXin(xInd+3, yInd, zInd)
                                   ;
            }
        }
        // The last part of the extrapolation
//        localResults[0] /= 16.0;
    }
    if(mesh->lastX()){
        vDotDSXout = v*dSXout;
        // Loop through the outer boundary points
        xInd = mesh->xend;
        for (yInd = mesh->ystart; yInd <= mesh->yend; yInd ++){
            for (zInd = 0; zInd < mesh->ngz - 1; zInd ++){
                // Extrapolate to the edge using a 4th order stencil
                localResults[1] += vDotDSXout(xInd, yInd, zInd)
//                localResults[1] +=   35.0*vDotDSXout(xInd  , yInd, zInd)
//                                   - 35.0*vDotDSXout(xInd-1, yInd, zInd)
//                                   + 21.0*vDotDSXout(xInd-2, yInd, zInd)
//                                   -  5.0*vDotDSXout(xInd-3, yInd, zInd)
                                   ;
            }
        }
        // The last part of the extrapolation
 //       localResults[1] /= 16.0;
    }
    if(mesh->firstY()){
        vDotDSYdown = v*dSYdown;
        // Loop through the lower boundary points
        yInd = mesh->ystart;
        for (xInd = mesh->xstart; xInd <= mesh->xend; xInd ++){
            for (zInd = 0; zInd < mesh->ngz - 1; zInd ++){
                localResults[2] += vDotDSYdown(xInd, yInd, zInd)
//                localResults[2] +=   35.0*vDotDSYdown(xInd, yInd  , zInd)
//                                   - 35.0*vDotDSYdown(xInd, yInd+1, zInd)
//                                   + 21.0*vDotDSYdown(xInd, yInd+2, zInd)
//                                   -  5.0*vDotDSYdown(xInd, yInd+3, zInd)
                                   ;
            }
        }
        // The last part of the extrapolation
//        localResults[2] /= 16.0;
    }
    if(mesh->lastY()){
        vDotDSYup = v*dSYup;
        // Loop through the upper boundary points
        yInd = mesh->yend;
        for (xInd = mesh->xstart; xInd <= mesh->xend; xInd ++){
            for (zInd = 0; zInd < mesh->ngz - 1; zInd ++){
                localResults[3] += vDotDSYup(xInd, yInd, zInd);
//                localResults[3] +=   35.0*vDotDSYup(xInd, yInd  , zInd)
//                                   - 35.0*vDotDSYup(xInd, yInd-1, zInd)
//                                   + 21.0*vDotDSYup(xInd, yInd-2, zInd)
//                                   -  5.0*vDotDSYup(xInd, yInd-3, zInd)
                                   ;
            }
        }
        // The last part of the extrapolation
//        localResults[3] /= 16.0;
    }

    output << "localResult.size()="<< localResults.size() << std::endl;
    // Sum the data from the processors, and store it in results
     for (std::vector<BoutReal>::iterator it = localResults.begin();
         it != localResults.end();
         ++it)
    {
        output << "localResult["<< it - localResults.begin()<<"]=" << *it << std::endl;
    }

    MPI_Allreduce(localResults.data(), results.data(),
                  4, MPI_DOUBLE, MPI_SUM, BoutComm::get());
}

#endif
