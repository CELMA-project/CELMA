#ifndef __OWNBCS_H__
#define __OWNBCS_H__

#include <bout.hxx>           // Includes all necessary classes and types
#include <bout/constants.hxx> // Gives PI and TWOPI

/*!
 * \class OwnBCs
 *
 * \brief Class which has member functions to set ghost points
 *
 * This class is used to set the boundary of Fields by specifying the ghost
 * points.
 *
 * \author Michael Løiten
 * \date 2016.12.04
 */
class OwnBCs
{
    private:
        // Data members
        int piIndex;          //!< Index corresponding to \f$\pi\f$
        int firstOuterXGhost; //!< Given that the processor has a boundary
        int firstUpperYGhost; //!< Given that the processor has a boundary
        int firstLowerYGhost; //!< Given that the processor has a boundary

        /*! If a warning is given rather than throwing an error if
         *  insufficient number of points is found
         */
        bool warnPoints;

        // Member functions
        //! A function which loops over the inner boundary
        void innerRhoCylinderLoop(Field3D &the_field,
                                  const int &y_start,
                                  const int &y_end);

        /**@{*/
        /*!
         * \brief Fieldgenerator for boundary values
         *
         * Analytic function used for boundary getting value at boundary
         * The member function "generate" of field data takes (x,y,z,t) as an
         * input, and returns the function evaluated in the given point in
         * time-space
         */
        FieldGenerator *aBndryFuncGen;
        FieldGenerator *bBndryFuncGen;
        /**@}*/
        BoutReal a; //!< The current a value
        BoutReal b; //!< The current b value
        //! The value of y at the boundary
        BoutReal yValAtYDownBndry;
        BoutReal x; //!< The value of x at the current index
        BoutReal z; //!< The value of z at the current index
    public:
        // Constructors
        OwnBCs();

        // Destructors
        // Not needed as we have no dynamic memory allocation yet

        // Member functions
        //! Specify ghost at inner \f$\rho\f$
        void innerRhoCylinder(Field3D &the_field);
        //! Specify ghost at outer \f$\rho\f$
        void extrapolateXOutGhost(Field3D &f);
        //! Extrapolate to first outer ghost
        void extrapolateYGhost(Field3D &f);
        //! Extrapolate to first upper ghosts
        void extrapolateYUp(Field3D &f);
        //! Extrapolate to second upper ghosts
        void extrapolateYUpScondGhost(Field3D &f);

        //! Extrapolate to first lower ghosts
        void extrapolateYDown(Field3D &f);
        //! Use sheath condition to set ghost point for ue
        void uEParSheath(Field3D &uEPar,
                      const Field3D &phi,
                      const BoutReal &Lambda,
                      const BoutReal &phiRef = 0.0,
                      const Field3D &profile = 1.0);
        //! Use sheath condition to set ghost point for the current
        void jParSheath(Field3D &jPar,
                        const Field3D &uEPar,
                        const Field3D &uIPar,
                        const Field3D &phi,
                        const Field3D &n,
                        const BoutReal &Lambda,
                        const BoutReal &phiRef = 0.0);
        //! Use a profiled sheath condition to set ghost point for the current
        void jParSheathProfiled(Field3D &jPar,
                                const Field3D &uEPar,
                                const Field3D &uIPar,
                                const Field3D &phi,
                                const Field3D &n,
                                const BoutReal &Lambda,
                                const BoutReal &phiRef = 0.0,
                                const Field3D &profile = 1.0);
        //! Use sheath condition to set ghost point for the par dens mom
        void parDensMomSheath(Field3D &momDensPar,
                              const Field3D &uIPar,
                              const Field3D &n);
        //! Profiled sheath condition to set ghost point for the par dens mom
        void parDensMomSheathProfiled(Field3D &momDensPar,
                                      const Field3D &uIPar,
                                      const Field3D &n,
                                      const Field3D &profile = 1.0);
        //! Specify the lower ghost from the Cauchy boundary condition
        void cauchyYDown(Field3D &f,
                         const BoutReal &t = 0.0,
                         bool const &yUpExtrapolate = true);
        /**@{*/
        //! Wrapper which initialize Cauchy
        void prepareCauchy(const string &section);
        //! Specify the lower ghost from the Cauchy boundary condition
        void getAFunction(const string &section);
        void getBFunction(const string &section);
        /**@}*/
        //! Get the value of y at the boundary
        void getYValAtYDownBndry();
};

// Function bodies of the non-inlined functions are located in the .cxx file
#include "../src/ownBCs.cxx"

#endif
