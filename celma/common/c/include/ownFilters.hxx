#ifndef __OWNFILTERS_H__
#define __OWNFILTERS_H__

#include <bout.hxx>
#include <fft.hxx>

class OwnFilters;
class AllPass;
class RadialLowPass;

// OwnFilters

/*!
 * \class OwnFilters
 *
 * \brief Parent class of filters
 *
 * Contains generic data and function for child filters
 *
 * \author Michael Løiten
 * \date 2016.08.06
 * \date 2016.08.09
 */
class OwnFilters
{
    public:
        // Constructor
        OwnFilters(Options *options);

        // Destructors
        // NOTE: New is called, so should destruct
        // NOTE: Needs to be virtual in order for the child classes to call destruct
        virtual ~OwnFilters();

        // Member functions
        //! Virtual function which is overridden by child classes
        /* NOTE: The = 0 is needed
         *       Tells the compiler that no function body will be given
         *       (pure virtual)
         */
        virtual const Field3D ownFilter(const Field3D &var) = 0;

        static OwnFilters* createFilter(Options *options = NULL);
    protected:
        // Data members
        //! Array containing the fourier transform
        dcomplex *fourierArray;
        BoutReal circumference; //!< Current circumference
        BoutReal lambdaMin;     //!< Shortest resolved wave length at outer circumference
        int ncz;                //!< Number of z-planes
        int kMax;               //!< Global maximum mode number
};

// AllPass

/*!
 * \class AllPass
 *
 * \brief No filtering performed.
 *
 * Simply returns the input field.
 *
 * Inherit from OwnFilters through public inheritance.
 *
 * \author Michael Løiten
 * \date 2016.08.09
 */
class AllPass : public OwnFilters
{
    public:
        //! Constructor does nothing
        AllPass(Options *options) : OwnFilters(options){};

        //! Destructor
        /* NOTE: {} in the end is needed
         *       If else the compiler gives
         *       "udefined reference to `vtable for ...'"
         */
        virtual ~AllPass(){};

        //! AllPass implementation of ownFilter
        const Field3D ownFilter(const Field3D &var);
};

// RadialLowPass

/*!
 * \class RadialLowPass
 *
 * \brief Lowpass filters in the \f$\theta\f$ direction
 *
 * This class provides filter, which goal is to prevent aliasing and unresolved
 * modes.
 *
 * Inherit from OwnFilters through public inheritance.
 *
 * \warning The following is only suitable if using a cylindrical Clebsch
 *          coordinate system
 *
 * \author Michael Løiten
 * \date 2016.08.09
 */
class RadialLowPass : public OwnFilters
{
    public:
        // Constructor
        RadialLowPass(Options *options);

        //! Destructor
        virtual ~RadialLowPass(){};

        //! RadialLowPass implementation of ownFilter
        const Field3D ownFilter(const Field3D &var);

    private:
        int kMaxCurrent;        //!< Max mode at current \f$\rho\f$
};

#include "../src/ownFilters.cxx"

#endif
