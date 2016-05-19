#ifndef __NOISEGENERATOR_H__
#define __NOISEGENERATOR_H__

#include <bout.hxx>           // Includes all necessary classes and types
#include <bout/constants.hxx> // Gives PI and TWOPI
#include <float.h>            // Includes DBL_EPSILON
#include <math.h>             // Includes fabs

/*!
 * \class NoiseGenerator
 *
 * \brief Add noise in cylinder geometry
 *
 * This class provides noise generators for adding noise in cylinder geometry
 * by adding a reasonable perturbation with random phases.
 * The noise is made on a cartesian grid, and then interpolated to a cylinder
 * grid
 *
 * \warning The following is only suitable if using a cylindrical Clebsch
 *          coordinate system
 *
 * \author Michael Løiten
 * \date 2016.18.05 (based on previous files)
 */
class NoiseGenerator
{
    private:
        // Data members
        Field3D noiseField; //!< Field which holds the noise

        /* Note about destructing vectors
         * It seems like vectors have their own memory management.
         * Recall to only use new and delete in pairs
         *
         * http://stackoverflow.com/questions/21568072/should-i-always-call-vector-clear-at-the-end-of-the-function
         * http://www.cplusplus.com/forum/general/7452/
         */
        //! Vector with the continuous variable x
        std::vector<BoutReal> x;
        //! Declare matrix as 2D vector: Should rather call M2D
        std::vector< std::vector <BoutReal> > M;

        int const seed;      //!< Seed to the random noise generator
        size_t const maxMode; //!< Max mode of perturbation
        BoutReal LRho;       //!< Radius is cylindrical geometry
        BoutReal Lx;         //!< x length in cartesian geometry
        BoutReal phiX;       //!< Phase in the (Cartesian) x direction
        BoutReal phiY;       //!< Phase in the (Cartesian) y direction
        BoutReal kScale;     //!< Scaling of the modes

        size_t numXYPoints;  //!< Number of points in x and y direction for the matrix

        // Member functions
        //! Creates a vector from 0 to totalLength with equidistant spacing
        void linspace(size_t const &nPoints,
                      BoutReal const &totalLength,
                      std::vector<BoutReal> &continuousVar
                      );

        //! Initializes the matrix M (called by the constructor)
        void initializeMatrix();

        //! Add noise where the modes are scaled with \f$k^{-5/2}\f$
        void add52Noise();

        //! Interpolation between the grids
        void interpolCartMatrixToCylField(Field3D &noiseField);

        //! Match the closest index in the field
        size_t getClosestFieldIndex(BoutReal const &matchMe);

    public:
        // Constructors
        NoiseGenerator (const string &section = "geom",
                        int const &seed = 1,
                        size_t const &maxMode = 14
                        );

        // Destructors
        // NOTE: No alloc or new called, so no reason to deallocate

        // Member functions
        //! Random phases generator
        void generateRandomPhases(Field3D &f, BoutReal const &scale);

};

// Function bodies of the non-inlined functions are located in the .cxx file
#include "../src/noiseGenerator.cxx"

#endif
