#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include <iostream>           // Gives the strings
#include <iomanip>            // Gives setw and setfill
#include <bout.hxx>
#include <bout/constants.hxx> // Gives PI and TWOPI

/*!
 * \class Parameters
 *
 * \brief Class which sets the normalized quantities
 *
 * \warning Assumes Hydrogen plasma
 *
 * \author Michael Løiten
 * \date 2016.08.24
 */
class Parameters
{
    private:
        BoutReal const len_;    //<! Cylinder length [m]
        BoutReal const radius_; //<! Plasma radius [m]
        BoutReal const n0_;     //<! Density normalization [m^-3]
        BoutReal const Te0_;    //<! Electron temperature [eV]
        BoutReal const Ti0_;    //<! Ion temperature [eV]
        BoutReal const B0_;     //<! Magnetic field [T]
        BoutReal const S_ ;     //<! Particle source [m^-3s^-1]
        BoutReal const nuEN_;   //<! Electron neutral frequency
        BoutReal const nuIN_;   //<! Ion neutral frequency

        // Converted
        BoutReal Ti0J;
        BoutReal Te0J;

        // Thermal quantities
        BoutReal cS;   //<! Ion sound speed
        BoutReal vA;   //<! Alfven velovity
        BoutReal vThE; //<! Thermal electron velocity
        BoutReal vThI; //<! Thermal electron velocity

        // Frequencies
        BoutReal omCI; //<! Ion cyclotron velocity
        BoutReal omCE; //<! Electron cyclotron velocity
        BoutReal omPE; //<! Plasma frequency

        // Sizes
        BoutReal debye;    //<! Debye length
        BoutReal eLarmour; //<! Electron Larmour radius
        BoutReal iLarmour; //<! Ion Larmour radius
        BoutReal rhoS;     //<! Hybrid radius
        BoutReal Lx;       //<! Normalized domain radius
        BoutReal Ly;       //<! Normalized domain length

        // Collisions
        BoutReal coloumbLog; //<! Coloumb logarithm
        BoutReal nuEI;       //<! Electron ion frequency
        BoutReal nuIE;       //<! Ion electron frequency
        BoutReal nuEE;       //<! Electron electron frequency
        BoutReal nuII;       //<! Ion ion frequency

        // Additional parameters
        BoutReal beta;   //<! Plasma beta
        BoutReal mu;     //<! Mass ratio
        BoutReal Lambda; //<! \f$\ln\left(\sqrt{\frac{\mu}{2\pi}}\rigth)\f$

        // Parallel viscosities
        BoutReal eta0I; //<! Viscosity parameter ion eta 0
        BoutReal eta2I; //<! Viscosity parameter ion eta 2
        BoutReal eta4I; //<! Viscosity parameter ion eta 4
        BoutReal eta0E; //<! Viscosity parameter electron eta 0
        BoutReal eta2E; //<! Viscosity parameter electron eta 2
        BoutReal eta4E; //<! Viscosity parameter electron eta 4

        // Normalized parameters
        BoutReal nuEINorm;  //<! Normalized nuEI
        BoutReal nuENNorm;  //<! Normalized nuEI
        BoutReal nuINNorm;  //<! Normalized nuEI
        BoutReal SNorm;     //<! Normalized S
        BoutReal eta0INorm; //<! Normalized eta0I
        BoutReal eta2INorm; //<! Normalized eta2I
        BoutReal eta4INorm; //<! Normalized eta4I
        BoutReal eta0ENorm; //<! Normalized eta0E
        BoutReal eta2ENorm; //<! Normalized eta2E
        BoutReal eta4ENorm; //<! Normalized eta4E

        // Text
        int const separatorLen; //<! Length of text separator
        char const separator;   //<! Separator character
        int const nameWidth;    //<! Text width for names
        int const numberWidth;  //<! Text width for numbers
        int const unitsWidth;   //<! Text width for units
        int const precision;    //<! The precision of the prints

    public:
        // Constructor
        Parameters(BoutReal const &radius,
                   BoutReal const &len,
                   BoutReal const &n0,
                   BoutReal const &Te0,
                   BoutReal const &Ti0,
                   BoutReal const &B0,
                   BoutReal const &S,
                   BoutReal const &nuEN,
                   BoutReal const &nuIN
                   );

        //! Prints the table
        void printTable() const;
        //! Variable printer
        void printVar(std::string const &name,
                      BoutReal const &val,
                      std::string const &units)
                      const;

        // Getters
        //! Obtain the domain radius
        BoutReal getLx() const;
        //! Obtain the domain length
        BoutReal getLy() const;
        //! Obtain the electron ion collision frequency
        BoutReal getNuEI() const;
        //! Obtain the electron neutral collision frequency
        BoutReal getNuEN() const;
        //! Obtain the ion neutral collision frequency
        BoutReal getNuIN() const;
        //! Obtain the normalized particle creation rate
        BoutReal getSNorm() const;
        //! Obtain the plasma beta
        BoutReal getBeta() const;
        //! Obtain the mass ratio
        BoutReal getMu() const;
        //! Obtain Lambda
        BoutReal getLambda() const;
        //! Obtain the ion gyration frequency
        BoutReal getOmCI() const;
        //! Obtain the hybrid radius
        BoutReal getRhoS() const;
};

#include "../src/parameters.cxx"

#endif
