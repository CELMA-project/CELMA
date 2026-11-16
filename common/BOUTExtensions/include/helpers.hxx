#ifndef COMMON_BOUTEXTENSIONS_INCLUDE_HELPERS_HXX_
#define COMMON_BOUTEXTENSIONS_INCLUDE_HELPERS_HXX_

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
class Helpers {
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
class PolAvg : private Helpers {
public:
  //! Function which returns the poloidal average of a field
  Field3D const poloidalAverage(Field3D const &f);
};

/*!
 * \class VolumeIntegral
 *
 * \brief Class containing volume integration
 *
 * \author Michael Løiten
 * \date 2016.08.19
 */
class VolumeIntegral : private Helpers {
public:
  //! Volume integral
  BoutReal volumeIntegral(Field3D const &f);
};

// Each model compiles one source file, so the header pulls in its
// implementation; cpplint assumes separately compiled .cxx files
#include "../src/helpers.cxx" // NOLINT(build/include)

#endif // COMMON_BOUTEXTENSIONS_INCLUDE_HELPERS_HXX_
