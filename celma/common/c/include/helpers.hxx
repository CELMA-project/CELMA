#ifndef __HELPERS_H__
#define __HELPERS_H__

#include <bout.hxx>

class PolAvg;
class VolumeIntegral;

/*!
 * \class Helpers
 *
 * \brief Base class containing helper functions used in monitoring and
 *        post-processing
 *
 * \author Michael Løiten
 * \date 2016.08.18
 */
class Helpers
{
    protected:
        int xInd; //!< x-looping index
        int yInd; //!< y-looping index
        int zInd; //!< z-looping index
};

/*!
 * \class PolAvg
 *
 * \brief Class which performs poloidal average
 *
 * \author Michael Løiten
 * \date 2016.08.19
 */
class PolAvg : private Helpers
{
    public:
        //! Function which returns the poloidal average of a field
        Field3D const polAvg(Field3D const &f);
};

/*!
 * \class volumeIntegral
 *
 * \brief Class containing volume integration
 *
 * \author Michael Løiten
 * \date 2016.08.19
 */
class VolumeIntegral : private Helpers
{
    public:
        //! Volume integral
        void volumeIntegral(Field3D const &f, BoutReal &result);
};

#include "../src/helpers.cxx"

#endif
