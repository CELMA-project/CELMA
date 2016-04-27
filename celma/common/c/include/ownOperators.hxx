#ifndef __OWNOPERATORS_H__
#define __OWNOPERATORS_H__

#include <bout.hxx>           // Includes all necessary classes and types
#include <bout/constants.hxx> // Gives PI and TWOPI

/*!
 * \class OwnOperators
 *
 * \brief Class with additional discretizing operators
 *
 * This class provides additional discretizing operators using a cylindrical
 * Clebsch coordinate system
 *
 * \warning The following is only suitable if using a cylindrical Clebsch
 *          coordinate system
 *
 * \author Michael Løiten
 * \date 2016.12.04
 */
class OwnOperators
{
    private:
        // Data members
        //! Should xboundary be calculated in D3DZ3
        bool incXBndry;
        //! Value of last X, used in D3DZ3
        BoutReal globalLastXVal;
        //! Function generator for D3DZ3
        /*!
         * Analytic function used for boundary getting value at boundary
         * The member function "generate" of field data takes (x,y,z,t) as an
         * input, and returns the function evaluated in the given point in
         * time-space
         */
        FieldGenerator *bndryFuncGen;
        Field2D J ; //!< The Jacobian
        Field2D J2; //!< The Jacobian raised to power 2
        Field2D J3; //!< The Jacobian raised to power 3
        Field2D J4; //!< The Jacobian raised to power 4
        // Derivatives used in div_uE_dot_grad_n_GradPerp_phi
        // These are declared here due to optimization
        /**@{*/
        //! \f$\rho\f$ derivatives of \f$\phi\f$
        Field3D phi_x,  phi_xx,  phi_xxx;
        /**@}*/
        //! \f$\theta\f$ derivatives of \f$\phi\f$
        /**@{*/
        Field3D phi_z,  phi_zz,  phi_zzz;
        /**@}*/
        /**@{*/
        //! Mixed derivatives of \f$\phi\f$
        Field3D phi_xz, phi_xxz, phi_xzz;
        /**@}*/
        /**@{*/
        Field3D n_x,    n_xx; //!< \f$\rho\f$ derivatives of \f$n\f$
        /**@}*/
        /**@{*/
        Field3D n_z,    n_zz; //!< \f$\theta\f$ derivatives of \f$n\f$
        /**@}*/
        Field3D n_xz; //!< Mixed derivatives of \f$n\f$

    public:
        // Constructors
        OwnOperators (const string &section = "phi");

        // Destructors
        // Unsure if this is needed here

        // Member functions
        //! Discretizes \f$\partial_\rho^3\f$
        Field3D D3DX3(const Field3D &f, const BoutReal &t = 0.0);
        //! Discretizes \f$\partial_\theta^3\f$
        Field3D D3DZ3(const Field3D &f);
        //! Operator for \f$\nabla\cdot_(f \nabla_\perp g)\f$
        Field3D div_f_GradPerp_g(const Field3D &f, const Field3D &g);
        //! Operator for \f$\nabla\cdot_(\mathbf{u}_e \nabla_\perp \phi)\f$
        Field3D div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                               const Field3D &phi);
        //! Operator for \f$\nabla_\perp f\f$ in cylinder geometry
        Vector3D Grad_perp(const Field3D &f);

        // Auxiliary functions
        //! Getter for IncXbndry
        bool getIncXbndry();
        //! Setter for IncXbndry
        void setIncXbndry(bool option);
};

// Function bodies of the non-inlined functions are located in the .cxx file
#include "../src/ownOperators.cxx"

#endif
