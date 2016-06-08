#ifndef __OWNFILTERS_H__
#define __OWNFILTERS_H__

#include <bout.hxx>
#include <fft.hxx>

/*!
 * \class OwnFilters
 *
 * \brief Lowpass filters in the \f$\theta\f$ direction
 *
 * This class provides filter, which goal is to prevent aliasing and unresolved
 * modes
 *
 * \warning The following is only suitable if using a cylindrical Clebsch
 *          coordinate system
 *
 * \author Michael Løiten
 * \date 2016.08.06
 */
class OwnFilters
{
    private:
        // Data members
        //! Array containing the fourier transform
        dcomplex *fourierArray;
        BoutReal circumference; //!< Current circumference
        BoutReal lambdaMin;     //!< Shortest resolved wave length at outer circumference
        int kMaxCurrent;        //!< Max mode at current \f$\rho\f$
        int ncz;                //!< Number of z-planes

    public:
        // Constructors
        OwnFilters();

        // Destructors
        // NOTE: New is called, so should destruct
        ~OwnFilters();

        // Member functions
        //! Bandpassfilter where the filter depends on \f$\rho\f$
        const Field3D radialLowPass(const Field3D &var);
};

#include "../src/ownFilters.cxx"

#endif
