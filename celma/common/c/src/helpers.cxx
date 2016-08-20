#ifndef __HELPERS_CXX__
#define __HELPERS_CXX__

#include "../include/helpers.hxx"

/*!
 * Returns the poloidal average of a field
 *
 * \param[in] f     The field to take the average of
 *
 * \returns result The poloidal average of the field
 */
Field3D const PolAvg::polAvg(const Field3D &f)
{
    TRACE("Halt in PolAvg::polAvg");

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
void VolumeIntegral::volumeIntegral(Field3D const &f, BoutReal &result)
{
    TRACE("Halt in VolumeIntegral::volumeIntegral");

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
 * \brief Constructor
 *
 * Constructor which sets the private member data used in the
 * surfaceEdgeIntegral. Note that localResults will be initialized to a vector
 * with 4 entries with value 0.0
 */
SurfaceIntegral::SurfaceIntegral() : vDotDSXin   (0.0),
                                     vDotDSXout  (0.0),
                                     vDotDSYdown (0.0),
                                     vDotDSYup   (0.0),
                                     localResults(4, 0.0),
                                     useInner_   (true),
                                     useOuter_   (true),
                                     useLower_   (true),
                                     useUpper_   (true)
{
    TRACE("SurfaceIntegral::SurfaceIntegral");

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

    MXG = mesh->getMXG();
    MYG = mesh->getMYG();
}

void SurfaceIntegral::setSurfaces(bool useInner,
                                  bool useOuter,
                                  bool useLower,
                                  bool useUpper)
{
    useInner_ = useInner;
    useOuter_ = useOuter;
    useLower_ = useLower;
    useUpper_ = useUpper;
}

/*!
 * Returns the surface integral of a vector field at the edge of the domain
 * using the rectangle rule in 3D. Note that the volume element \f$dS\f$
 * is not spanned by four points, but rather the surface encapsulating
 * the point under consideration (i.e. by the four lines running
 * parallelly with the coordinate lines, but is half a grid size away
 * from the point in the four directions)
 *
 * We have that the surface integral is given by (see D'Haeseleer chapter 2.5.49
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
 * \param[in] xIndInner Global index for the x-index to use for the -x surface
 * \param[in] xIndOuter Global index for the x-index to use for the +x surface
 * \param[in] yIndLower Global index for the y-index to use for the -y surface
 * \param[in] yIndUpper Global index for the y-index to use for the +y surface
 * \param[in] results   std::vector to store the results in\n
 *                      result[0] - integral over v on the xin surface\n
 *                      result[1] - integral over v on the xout surface\n
 *                      result[2] - integral over v on the ydown surface\n
 *                      result[3] - integral over v on the yup surface
 *
 *
 * \param[out] results  The result of the integration\n
 *                      result[0] - integral over v on the xin surface\n
 *                      result[1] - integral over v on the xout surface\n
 *                      result[2] - integral over v on the ydown surface\n
 *                      result[3] - integral over v on the yup surface
 *
 * \note There are no integral on the z-surfaces as these are periodic.
 */
void SurfaceIntegral::surfaceEdgeIntegral(Vector3D const &v,
                                          int const &xIndInner,
                                          int const &xIndOuter,
                                          int const &yIndLower,
                                          int const &yIndUpper,
                                          std::vector<BoutReal> &results)
{
    TRACE("Halt in SurfaceIntegral::surfaceEdgeIntegral");

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
    /* NOTE: YLOCAL/YGLOBAL different than XLOCAL/XGLOBAL
     *       They differ by an extra +/- MYG in YLOCAL/YGLOBAL
     *       As the input in this function are pure global indices
     *       (which counts from the first ghost through the inner
     *       points to the last ghost points), the YLOCAL/YGLOBAL will
     *       be subtracted/added with MYG to add up for the discrepancy
     */
    if(mesh->IS_MYPROC(xIndInner, mesh->YGLOBAL(mesh->ystart + MYG)) && useInner_){
        // Cast to local indices
        xIndInnerLoc = mesh->XLOCAL(xIndInner);
        // Loop through the inner surface
        vDotDSXin = v*dSXin;
        for (yInd = mesh->ystart; yInd <= mesh->yend; yInd ++){
            for (zInd = 0; zInd < mesh->ngz - 1; zInd ++){
                localResults[0] += vDotDSXin(xIndInnerLoc, yInd, zInd);
            }
        }
    }
    if(mesh->IS_MYPROC(xIndOuter, mesh->YGLOBAL(mesh->ystart + MYG)) && useOuter_){
        // Cast to local indices
        xIndOuterLoc = mesh->XLOCAL(xIndOuter);
        // Loop through the outer surface
        vDotDSXout = v*dSXout;
        for (yInd = mesh->ystart; yInd <= mesh->yend; yInd ++){
            for (zInd = 0; zInd < mesh->ngz - 1; zInd ++){
                localResults[1] += vDotDSXout(xIndOuterLoc, yInd, zInd);
            }
        }
    }
    if(mesh->IS_MYPROC(mesh->XGLOBAL(mesh->xstart), yIndLower) && useLower_){
        // Cast to local indices
        yIndLowerLoc = mesh->YLOCAL(yIndLower - MYG);
        // Loop through the lower surface
        vDotDSYdown = v*dSYdown;
        for (xInd = mesh->xstart; xInd <= mesh->xend; xInd ++){
            for (zInd = 0; zInd < mesh->ngz - 1; zInd ++){
                localResults[2] += vDotDSYdown(xInd, yIndLowerLoc, zInd);
            }
        }
    }
    if(mesh->IS_MYPROC(mesh->XGLOBAL(mesh->xstart), yIndUpper) && useUpper_){
        // Cast to local indices
        yIndUpperLoc = mesh->YLOCAL(yIndUpper -MYG);
        // Loop through the upper surface
        vDotDSYup = v*dSYup;
        for (xInd = mesh->xstart; xInd <= mesh->xend; xInd ++){
            for (zInd = 0; zInd < mesh->ngz - 1; zInd ++){
                localResults[3] += vDotDSYup(xInd, yIndUpperLoc, zInd);
            }
        }
    }

    MPI_Allreduce(localResults.data(), results.data(),
                  4, MPI_DOUBLE, MPI_SUM, BoutComm::get());
}

#endif
