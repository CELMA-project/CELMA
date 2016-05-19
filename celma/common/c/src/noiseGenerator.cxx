#ifndef __NOISEGENERATOR_CXX__
#define __NOISEGENERATOR_CXX__

#include "../include/noiseGenerator.hxx"

/*!
 * \brief Constructor
 *
 * Constructor which sets the private member data
 *
 * \param[in] section   What section to get bndryFuncGen from.
 *                      Default is "phi"
 * \param[in] seed      Number to seed the random number generator with
 * \param[in] maxMode   Max mode number of perturbations
 *
 */
NoiseGenerator::NoiseGenerator(const string &section, int const &seed, size_t const &maxMode) :
    seed(seed), maxMode(maxMode)
{
    TRACE("Halt in NoiseGenerator::NoiseGenerator");

    // Get LRho
    Options *geomOptions = Options::getRoot()->getSection(section);
    geomOptions->get("Lx", LRho, 0.0);

    // NOTE: Take care when comparing floats
    // http://stackoverflow.com/questions/17333/most-effective-way-for-float-and-double-comparison
    if ( fabs(LRho - 0.0) < DBL_EPSILON ){
        output << "Searched in section '" << section << "' with Lx = " << Lx << std::endl;
        output << "fabs(LRho - 0.0) = " << fabs(LRho - 0.0) <<std::endl;
        output << "DBL_EPSILON = " << DBL_EPSILON << std::endl;
        throw BoutException("No Lx found, or Lx set to 0.0!\n");
    }

    // LRho read from the input file is the radius. We must thus multiply this
    // with 2 in order to get the Lx length.
    // Note also that we let Lx = Ly in the cartesian grid.
    Lx = LRho*2.0;

    // Get the number of points in x and y direction for the matrix
    /* NOTE: Number of points for the matrix
     *       Chosen large in order to have a lot of
     *       matrix points around the meshpoint for Field3D f.
     */
    numXYPoints = 5*mesh->GlobalNx;

    // Make a linspace of x
    linspace(numXYPoints, Lx, x);

    // Initialize the matrix
    initializeMatrix();

    // Seed the random number generator
    // http://stackoverflow.com/questions/686353/c-random-float-number-generation
    srand (static_cast <unsigned> (seed));

    // Set the noiseField to 0
    noiseField = 0.0;
}

/*!
 * Wrapper fundtion which adds the random noise
 * Explanation of the procedure:
 * 1. Dynamically create a matrix M
 * 2. Fill it with modes from 1 to maxMode (the modes will have random phases)
 * 3. Interpolate the values from the matrix onto Field3D f
 *
 * \param[in] f         The original field
 * \param[in] scale     Multiplicative factor to f
 *
 * \param[out] f        The field after filling it with random phases
 *
 * \note No parallelization of the big matrix
 *       The matrix will be in cartesian coordinates, whereas our mesh is
 *       in cylindrical coordinates. Thus, if we wanted to parallelize the
 *       matrix, we would need to do that according to the coordinates of
 *       the mesh.
 *       In the case of a cylindrical mesh (which would be discretized in
 *       rho) we would need to discretize the matrix as an "annulus" and an
 *       overlapping squares.
 */
void NoiseGenerator::generateRandomPhases(Field3D &f, BoutReal const &scale)
{
    TRACE("Halt in: NoiseGenerator::generateRandomPhases");

    // Add noise to a 2D cartesian grid
    add52Noise();

    // Interpolate the matrix
    interpolCartMatrixToCylField(noiseField);

    // Multiply cos(rho*0.5*pi/rhoLen)
    // The Jacobian in cylinder geometry equals rho
    noiseField *= cos(mesh->J*0.5*PI/LRho);

    // Normalize (the "true" in the end means that the max will be taken for
    // all processors)
    noiseField /= abs(noiseField.max(true));

    // Scale with the proper scale
    noiseField *= scale;

    // Add the old value
    f += noiseField;
}

/*!
 * A linspace function which starts from 0.0, and ends at totalLength
 *
 * \param[in] nPoints       Number of points
 * \param[in] totalLength   Endpoint
 * \param[in] continuousVar Variable to set the linspace to
 *
 * \param[out] continuousVar The linspace given as an array
 *
 * \note Instead calling by reference and manipulating the reference, it would
 *       be possible to return by reference (see for example
 *       http://stackoverflow.com/questions/5033627/static-variable-inside-of-a-function-in-c
 *       this is however not chosen as it seems like the variable is copied to
 *       to the LHS anyways.
 */
void NoiseGenerator::linspace(size_t const &nPoints,
                              BoutReal const &totalLength,
                              std::vector<BoutReal> &continuousVar
                              )
{

    /* NOTE: Returning by reference
     * http://www.tutorialspoint.com/cplusplus/returning_values_by_reference.htm
     */
    TRACE("Halt in NoiseGenerator::linspace");

    BoutReal gridSpacing;

    // Calculate gridSpacing
    gridSpacing = totalLength/(nPoints - 1);

    // Reserve memory for the vector (does not call the constructor)
    continuousVar.reserve(nPoints);

    // As continuousVar currently has no size yet, we cannot use the iterator to loop over
    // the vector, falling back to standard loop over numbers
    for(size_t xInd = 0; xInd < nPoints ; ++xInd){
        continuousVar.push_back(xInd*gridSpacing);
    }

}

/*!
 * Initialized the matrix M with 0's
 */
void NoiseGenerator::initializeMatrix()
{

    TRACE("Halt in NoiseGenerator::initializeMatrix");

    M.reserve(numXYPoints);

    // As matrix currently has no size yet, we cannot use the iterator to loop over
    // the rows, we will therefore fall back to standard loop over numbers
    // Create a column of size numXYPoints, and fill it with 0
    std::vector<BoutReal> cols(numXYPoints, 0.0);
    for(size_t xInd = 0; xInd < numXYPoints ; ++xInd){
        M.push_back(cols);
    }
}

// FIXME: Write why we set it to -5/2
/*!
 * Adds noise to a cartesian grid.
 *
 * The noise here consists of a superposition of maxMode number of waves, each
 * with a random phase.
 *
 * Each mode is scaled with \f$k^{-5/2}\f$, where \f$k\f$ is the mode number
 */
void NoiseGenerator::add52Noise()
{

    TRACE("Halt in NoiseGenerator::add52Noise");

    for(int modeNr = 1; modeNr <= maxMode; modeNr++){
        // Generating the random phases
        phiX = static_cast <BoutReal>
             (rand()) / (static_cast <BoutReal> (RAND_MAX/(2.0*M_PI)));
        phiY = static_cast <BoutReal>
             (rand()) / (static_cast <BoutReal> (RAND_MAX/(2.0*M_PI)));
        for(size_t xInd = 0; xInd < M.size(); ++xInd){
            for(size_t yInd = 0; yInd < M[xInd].size(); ++yInd){
                // We scale the modes with k^(-5/2)
                kScale = pow(modeNr,-5.0/2.0);
                // Making the M
                M[xInd][yInd] +=
                    kScale*sin(modeNr*(x[xInd]*2.0*M_PI/Lx) + phiX)
                  + kScale*sin(modeNr*(x[yInd]*2.0*M_PI/Lx) + phiY);
            }
        }
    }
}

/*!
 * Interpolates the values in the cartesian grid to the cylindrical field.
 * A bilinear interpolation is chosen.
 *
 * \param[in] noiseField     The field to be interpolated
 *
 * \param[out] noiseField    The interpolated field
 */
void NoiseGenerator::interpolCartMatrixToCylField(Field3D &noiseField)
{

    TRACE("Halt in NoiseGenerator::interpolCartMatrixToCylField");

    // Declare variables
    BoutReal circleX; // Corresponding x coordinate from [rho, theta]
    BoutReal circleY; // Corresponding y coordinate from [rho, theta]
    size_t x1; // Index in x closest to (before or equal) current circleX value
    size_t x2; // x1 + 1
    size_t y1; // Index in x closest to (before or equal) current circleY value
    size_t y2; // y1 + 1

    // Loop through all the inner points
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
    for (size_t xInd = mesh->xstart; xInd <= mesh->xend; xInd ++){
        for (size_t yInd = mesh->ystart; yInd <= mesh->yend; yInd ++){
            for (size_t zInd = 0; zInd < mesh->ngz - 1; zInd ++){

                // Calculate the x, y variables from rho, theta variables
                /* NOTE: mesh->GlobalX
                 *       mesh->GloblaX returns the continous variable from dx/2
                 *       to 1 - dx/2
                 */
                /* NOTE: Calculating rho and theta
                 *      rho is given by GlobalX*LRho
                 *      theta is given by zInd*mesh->dz
                */
                // Shift the center by adding LRho
                circleX = mesh->GlobalX(xInd)*LRho
                           *cos(zInd*mesh->dz)
                           + LRho;
                circleY = mesh->GlobalX(xInd)*LRho
                           *sin(zInd*mesh->dz)
                           + LRho;

                // Find index in x closest to (before or equal) cur circleX
                // value
                x1 = getClosestFieldIndex(circleX);
                // If we are on the end of the vector (minus one since indices
                // counts from 0)
                if (x1 == x.size() - 1){
                    // Go one index back
                    x1 -= 1;
                }
                // Find index in x closest to (before or equal) cur circleX
                // value
                x2 = x1 + 1;
                // We find the closest y index relative to the current index (but
                // before the current index)
                y1 = getClosestFieldIndex(circleY);
                // If we are on the end of the vector (minus one since indices
                // counts from 0)
                if (y1 == x.size() - 1){
                    // Go one index back
                    y1 -= 1;
                }
                // We find the closest y index relative to the current index (but
                // after the current index)
                y2 = y1 + 1;

                // Making the interpolation
                noiseField(xInd, yInd, zInd) =
                                  (1.0/(
                                         (x[x2] - x[x1])*
                                         (x[y2] - x[y1])
                                        )
                                   )*
                                   (
                                     M[x1][y1]*(
                                         (x[x2] - circleX)*
                                         (x[y2] - circleY)
                                                  )
                                   + M[x2][y1]*(
                                         (circleX - x[x1])*
                                         (x[y2]  - circleY)
                                                  )
                                   + M[x1][y2]*(
                                         (x[x2]  - circleX)*
                                         (circleY - x[y1])
                                                  )
                                   + M[x2][y2]*(
                                         (circleX - x[x1])*
                                         (circleY - x[y1])
                                                  )
                                   );
            }
        }
    }
}

/*!
 * Gets the closest index in the field corresponding to the input
 *
 * \param[in] matchMe   Either x or y coordinates from the \f$\rho\f$,
 *                      \f$\theta\f$ coordinates
 *
 * \returns closestInd   The closest index
 */
size_t NoiseGenerator::getClosestFieldIndex(BoutReal const &matchMe)
{

    TRACE("Halt in NoiseGenerator::getClosestFieldIndex");
    // Initialize closestInd and currentInd
    size_t closestInd = 0;
    size_t currentInd = 0;

    // Initialize the distances
    // Start with the furthest possible distance
    BoutReal bestDistance = Lx;
    BoutReal newDistance  = 0.0;

    // Loop over the vector
    for (std::vector<BoutReal>::const_iterator it = x.begin();
         it != x.end();
         ++ it){
        // Obtain the current index (does arithmetics on the iterator [like a
        // pointer] and cast it to and integer)
        // http://stackoverflow.com/questions/2152986/best-way-to-get-the-index-of-an-iterator
        currentInd = it - x.begin();
        // The current element is a reference to the adress pointed to by 'it'
        // Changing 'it', will also change the address at that adress. Therefore,
        // we add the const qualifier
        BoutReal const &curElem = *it;
        newDistance = fabs(matchMe - curElem);

        // If new distance is less than the bestDistance
        if (newDistance <= bestDistance){
            // Update the closest index
            closestInd = currentInd;
            // Update the best distance
            bestDistance = newDistance;
        }
        else{
            // The new distance is larger than the bestDistance.
            // We are moving away from target
            break;
        }
    }

    // Return the closest index
    return closestInd;
}

#endif
