#ifndef __OWNMONITORS_H__
#define __OWNMONITORS_H__

#include <bout.hxx>
#include <map>
#include "helpers.hxx"

/*!
 * \class OwnMonitors
 *
 * \brief Class containing monitor functions
 *
 * \author Michael Løiten
 * \date 2016.08.20
 */
class OwnMonitors
{
private:
    VolumeIntegral volInt_;      //<! Volume integration object
public:
    // Functions
    //! Calculates the kinetic energy
    void kinEnergy(Field3D  const &n                    ,
                   Vector3D const &GradPerpPhi          ,
                   Field3D  const &uEPar                ,
                   Field3D  const &uIPar                ,
                   std::map<std::string, BoutReal> *kinE);

    //! Calculates the kinetic energy
    void numberOfParticles(Field3D const &n, std::map<std::string, BoutReal> *N);
};


#include "../src/ownMonitors.cxx"

#endif
