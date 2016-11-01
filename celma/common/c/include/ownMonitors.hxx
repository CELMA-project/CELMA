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
    PolAvg polAvg_;              //<! Poloidal average object
    Field3D polAvgN_;            //<! Poloidally averaged density
    Field3D polAvgLogN_;         //<! Poloidally averaged of the log of n
    Vector3D polAvgGradPerpPhi_; //<! Poloidally averaged $\f\nabla_\perp\phi$\f
    Field3D polAvgUEPar_;        //<! Poloidally averaged parallel el velocity
    Field3D polAvgUIPar_;        //<! Poloidally averaged parallel ion velocity
public:
    // Constructors
    OwnMonitors();
    // Functions
    //! Calculates the poloidal density
    void calcPolAvgN(Field3D  const &n);

    //! Calculates the kinetic energy
    void kinEnergy(Field3D  const &n                    ,
                   Vector3D const &GradPerpPhi          ,
                   Field3D  const &uEPar                ,
                   Field3D  const &uIPar                ,
                   std::map<std::string, BoutReal> *kinE);

    //! Calculates the kinetic energy
    void potEnergy(Field3D const &n, std::map<std::string, BoutReal> *potE);

    //! Calculates the kinetic energy
    void totalN(Field3D const &n, std::map<std::string, BoutReal> *totN);
};


#include "../src/ownMonitors.cxx"

#endif
