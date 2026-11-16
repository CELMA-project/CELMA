#ifndef COMMON_BOUTEXTENSIONS_INCLUDE_OWNMONITORS_HXX_
#define COMMON_BOUTEXTENSIONS_INCLUDE_OWNMONITORS_HXX_

#include "helpers.hxx"
#include <bout.hxx>
#include <map>
#include <string>

/*!
 * \class OwnMonitors
 *
 * \brief Class containing monitor functions
 *
 * \author Michael Løiten
 * \date 2016.08.20
 */
class OwnMonitors {
private:
  VolumeIntegral volInt_; //<! Volume integration object
public:
  // Functions
  //! Calculates the kinetic energy
  void kinEnergy(Field3D const &n, Vector3D const &GradPerpPhi,
                 Field3D const &uEPar, Field3D const &uIPar,
                 std::map<std::string, BoutReal> *kinE);

  //! Calculates the kinetic energy
  void numberOfParticles(Field3D const &n, std::map<std::string, BoutReal> *N);
};

// Each model compiles one source file, so the header pulls in its
// implementation; cpplint assumes separately compiled .cxx files
#include "../src/ownMonitors.cxx" // NOLINT(build/include)

#endif // COMMON_BOUTEXTENSIONS_INCLUDE_OWNMONITORS_HXX_
