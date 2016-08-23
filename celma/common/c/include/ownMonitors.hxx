#ifndef __OWNMONITORS_H__
#define __OWNMONITORS_H__

#include <bout.hxx>
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
    VolumeIntegral volInt;
    SurfaceIntegral surfInt;
public:
    // Functions
    //! Calculates the kinetic energy
    void kinEnergy(Field3D  const &n          ,
                   Vector3D const &GradPerpPhi,
                   Field3D  const &uPar       ,
                   std::vector<BoutReal> *kinE);

    //! Calculates the outflow rate
    void outflowRate(Field3D const &f, BoutReal *outflowR);
};


#include "../src/ownMonitors.cxx"

#endif
