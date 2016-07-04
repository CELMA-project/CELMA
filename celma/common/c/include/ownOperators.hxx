#ifndef __OWNOPERATORS_H__
#define __OWNOPERATORS_H__

#include <bout.hxx>           // Includes all necessary classes and types
#include <bout/constants.hxx> // Gives PI and TWOPI
#include "ownBCs.hxx"         // Gives inner rho boundaries

class OwnOperators;
class OwnOpSimpleStupid;
class OwnOpOnlyBracket;
class OwnOp2Brackets;
class OwnOp3Brackets;
class OwnOpBasicBrackets;

// OwnOperators

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
        /*! If a warning is given rather than throwing an error if
         *  insufficient number of points is found
         */
        bool warnPoints;
    protected:
        // Data members
        FieldGenerator *bndryFuncGen;
        Field2D J ; //!< The Jacobian
        Field2D J2; //!< The Jacobian raised to power 2
        Field2D J3; //!< The Jacobian raised to power 3

        Field2D invJ ; //!< The inverse of the Jacobian
        Field2D invJ2; //!< The inverse of the Jacobian raised to power 2
        Field2D invJ3; //!< The inverse of the Jacobian raised to power 3
        Field2D invJ4; //!< The inverse of the Jacobian raised to power 4

        int xInd; //!< x-index
        int yInd; //!< y-index
        int zInd; //!< z-index
    public:
        // Constructors
        OwnOperators (Options *options);

        // Destructors
        // NOTE: New is called, so should destruct
        // NOTE: Needs to be virtual in order for the child classes to call destruct
        virtual ~OwnOperators(){};

        // Member functions
        //! Discretizes \f$\partial_\theta^3\f$
        Field3D D3DZ3(const Field3D &f);
        //! Discretizes \f$\partial_\rho^4\f$
        Field3D D4DX4(const Field3D &f);
        //! Discretizes \f$\partial_\|^4\f$
        Field3D D4DY4(const Field3D &f);
        //! Operator for \f$\nabla\cdot_(f \nabla_\perp g)\f$
        Field3D div_f_GradPerp_g(const Field3D &f, const Field3D &g);
        //! Operator for \f$\nabla_\perp f\f$ in cylinder geometry
        Vector3D Grad_perp(const Field3D &f);

        /* NOTE: Child classes can have new memberfunctions...
         *       which are not declared in the parent class.
         *       Pure virtual functions (has an = 0.0 in the end of the
         *       declaration) are used for functions which MUST be
         *       implemented
         *       If no default implementation is given, the declaration can end
         *       with an empty function body {}
         */
        /*! Operator for
         * \f$\nabla\cdot(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$
         * (only used in the simpleStupid implementation)
         */
        virtual Field3D div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                                       const Field3D &phi) = 0;
        //! Operator for \f$\{\phi, \Omega^D\}\f$
        virtual Field3D vortDAdv (const Field3D &phi, const Field3D &vortD) = 0;
        /*! Operator for
         * \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$
         * (only used in 2Brackets)
         */
        virtual Field3D kinEnAdvN(const Field3D &phi, const Field3D &n) = 0;
        /*! Operator for D3DX3 (must be declared virtual if it is to be
         * publicly accessable in just one of the child classes)
         */
        virtual Field3D D3DX3(const Field3D &f, const BoutReal &t = 0.0);
        /*! Operator for \f$\{(\partial_\rho\phi)^2, n\}\f$ (must be declared virtual if it is to be
         * publicly accessable in just one of the child classes)
         */
        virtual Field3D ArakawaOfDDXPhi2N(Field3D const &phi, Field3D const &n);

        // Auxiliary functions
        //! Getter for IncXbndry
        bool getIncXbndry();
        //! Setter for IncXbndry
        void setIncXbndry(bool option);

        //! Factory which chooses child class
        static OwnOperators* createOperators(Options *options = NULL);
};

// OwnOpSimpleStupid

/*!
 * \class OwnOpSimpleStupid
 *
 * \brief Simple stupid implementation of
 *        \f$\nabla\cdot_(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$
 *
 * The implementation is just using normal finite differences. No Arakawa
 * brackets is used. There are indications that this scheme is not energy
 * conserving.
 *
 * Inherit from OwnOperators through public inheritance.
 *
 * \warning This implementation has been found to create high \f$k\f$ structures
 *
 * \warning The following is only suitable if using a cylindrical Clebsch
 *          coordinate system
 *
 * \author Michael Løiten
 * \date 2016.07.13
 */
class OwnOpSimpleStupid : public OwnOperators
{
    private:
        //! Value of last X, used in D3DZ3
        BoutReal globalLastXVal;
        //! Function generator for D3DZ3
        /*!
         * Analytic function used for boundary getting value at boundary
         * The member function "generate" of field data takes (x,y,z,t) as an
         * input, and returns the function evaluated in the given point in
         * time-space
         */
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
        Field3D n_xz;   //!< Mixed derivatives of \f$n\f$
    public:
        // Constructors
        OwnOpSimpleStupid(Options *options, string &phiBndrySec);

        // Member functions
        /*! Operator for
         * \f$\nabla\cdot_(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$
         */
        Field3D div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                               const Field3D &phi);
        //! Will throw error
        Field3D vortDAdv (const Field3D &phi, const Field3D &vortD);
        //! Will throw error
        Field3D kinEnAdvN(const Field3D &phi, const Field3D &n);

        // Member functions
        //! Discretizes \f$\partial_\rho^3\f$
        Field3D D3DX3(const Field3D &f, const BoutReal &t = 0.0);

        //! Destructor
        /* NOTE: {} in the end is needed
         *       If else the compiler gives
         *       "udefined reference to `vtable for ...'"
         */
        virtual ~OwnOpSimpleStupid(){};
};

// OwnOpOnlyBracket

/*!
 * \class OwnOpOnlyBracket
 *
 * \brief Inconsistent implementation of
 *        \f$\nabla\cdot_(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$
 *
 * In this implementation, only \f$\{\phi, Omega^D\}\f$ is used, which is not
 * consistent.
 *
 * Inherit from OwnOperators through public inheritance.
 *
 * \warning This implementation is not consistent.
 *
 * \warning The following is only suitable if using a cylindrical Clebsch
 *          coordinate system
 *
 * \author Michael Løiten
 * \date 2016.07.13
 */
class OwnOpOnlyBracket : public OwnOperators
{
    private:
        BRACKET_METHOD bm;   //!< The bracket method
    public:
        // Constructors
        OwnOpOnlyBracket(Options *options);

        //! Will throw error
        Field3D div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                               const Field3D &phi);


        //! Operator for \f$\{\phi, Omega^D\}\f$
        Field3D vortDAdv (const Field3D &phi, const Field3D &vortD);
        //! Will throw error
        Field3D kinEnAdvN(const Field3D &phi, const Field3D &n);

        //! Destructor
        /* NOTE: {} in the end is needed
         *       If else the compiler gives
         *       "udefined reference to `vtable for ...'"
         */
        virtual ~OwnOpOnlyBracket(){};
};

// OwnOp2Brackets

/*!
 * \class OwnOp2Brackets
 *
 * \brief Implementation of
 *        \f$\nabla\cdot_(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$
 *        using 2 brackets, and one non-bracket term
 *
 * This implementation uses brackets for advection of \f$\Omega^D\f$, and
 * splits the \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$ into an
 * Arakawa part and a finite differences part. Finite differences are used as
 * two consecutive radial derivatives are not converging (but diverging) when
 * the grid space is decreasing.
 *
 * Inherit from OwnOperators through public inheritance.
 *
 * \note Looks like this implementation gives small scales (what is called
 *       "fake noodling" in Arakawas original paper)
 *
 * \warning The following is only suitable if using a cylindrical Clebsch
 *          coordinate system
 *
 * \author Michael Løiten
 * \date 2016.07.13
 */
class OwnOp2Brackets : public OwnOperators
{
    private:
        BRACKET_METHOD bm;   //!< The bracket method
    public:
        // Constructors
        OwnOp2Brackets(Options *options);

        //! Will throw error
        Field3D div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                               const Field3D &phi);
        //! Operator for \f$\{\phi, Omega^D\}\f$
        Field3D vortDAdv (const Field3D &phi, const Field3D &vortD);
        /*! Operator for
         * \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$
         * (only used in 2Brackets)
         */
        Field3D kinEnAdvN(const Field3D &phi, const Field3D &n);

        //! Destructor
        /* NOTE: {} in the end is needed
         *       If else the compiler gives
         *       "udefined reference to `vtable for ...'"
         */
        virtual ~OwnOp2Brackets(){};
};

// OwnOp3Brackets

/*!
 * \class OwnOp3Brackets
 *
 * \brief Implementation of
 *        \f$\nabla\cdot_(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$
 *        using 3 brackets
 *
 * This implementation uses brackets for advection of \f$\Omega^D\f$.
 * Further
 *
 * \f{eqnarray}{
 * \frac{B}{2}
 * \{\mathbf{u}_E\cdot\mathbf{u}_E, n\}
 * =
 * \frac{1}{J2}
 * \{(\partial_\rho\phi)^2+\frac{1}{\rho^2}(\partial_\theta \phi)^2, n\}
 * \f}
 *
 * This implemtation uses a modified Arakawa bracket for the first term. This
 * approach makes the \f$\partial_\rho^2 f\f$ term in the bracket resemble the
 * D2DX2(f) operator rather than DDX(DDX(f)). A standard Arakawa bracket for
 * the second term.
 *
 * Inherit from OwnOperators through public inheritance.
 *
 * \warning The following is only suitable if using a cylindrical Clebsch
 *          coordinate system
 *
 * \author Michael Løiten
 * \date 2016.07.22
 */
class OwnOp3Brackets : public OwnOperators
{
    private:
        BRACKET_METHOD bm;   //!< The bracket method

        /**@{*/
        /*! Index used in Arakawa loop
         *  xP1 denotes x + 1 from the point under consideration, whereas
         *  zM1 denotes z - 1 from the point under consideration
         */
        int xInd;
        int xIndP1;
        int xIndM1;
        int ncz;
        int zIndP1;
        int zIndM1;
        int xstart;
        int xend;
        /**@}*/
        /**@{*/
        /*! Stencil points of the Arakawa bracket
         *  xP1 denotes x + 1 from the point under consideration, whereas
         *  zM1 denotes z - 1 from the point under consideration
         */
        BoutReal phi_xyzP1;
        BoutReal phi_xyzM1;
        BoutReal phi_xP1yz;
        BoutReal phi_xM1yz;
        BoutReal phi_xP1yzP1;
        BoutReal phi_xP1yzM1;
        BoutReal phi_xM1yzP1;
        BoutReal phi_xM1yzM1;
        /**@}*/
        /**@{*/
        //! Ways of combining the stencils (explained in Arakawa's paper)
        BoutReal Jpp;
        BoutReal Jpx;
        BoutReal Jxp;
        /**@}*/

    public:
        // Constructors
        OwnOp3Brackets(Options *options);

        //! Will throw error
        Field3D div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                               const Field3D &phi);
        //! Operator for \f$\{\phi, \Omega^D\}\f$
        Field3D vortDAdv (const Field3D &phi, const Field3D &vortD);
        /*! Operator for
         * \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$
         * (only used in 2Brackets)
         */
        Field3D kinEnAdvN(const Field3D &phi, const Field3D &n);

        //! Arakawa bracket for \f$\{(\partial_\rho\phi)^2, n\}\f$
        Field3D ArakawaOfDDXPhi2N(Field3D const &f, Field3D const &g);

        //! Destructor
        /* NOTE: {} in the end is needed
         *       If else the compiler gives
         *       "udefined reference to `vtable for ...'"
         */
        virtual ~OwnOp3Brackets(){};
};

// OwnOpBasicBrackets

/*!
 * \class OwnOpBasicBrackets
 *
 * \brief Implementation of
 *        \f$\nabla\cdot_(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$
 *        using 3 basic brackets
 *
 * This implementation uses brackets for advection of \f$\Omega^D\f$.
 * Further
 *
 * \f{eqnarray}{
 * \frac{B}{2}
 * \{\mathbf{u}_E\cdot\mathbf{u}_E, n\}
 * =
 * \frac{1}{J2}
 * \{(\partial_\rho\phi)^2+\frac{1}{\rho^2}(\partial_\theta \phi)^2, n\}
 * \f}
 *
 * This implemtation uses an Arakawa bracket for these two terms.
 * A field \f$f=\partial_\rho\phi\f$ is used as an input for the first term,
 * where the boundary condition \f$f\f$ has been re-applied.
 *
 * Inherit from OwnOperators through public inheritance.
 *
 * \warning The following is only suitable if using a cylindrical Clebsch
 *          coordinate system
 *
 * \author Michael Løiten
 * \date 2016.07.23
 */
class OwnOpBasicBrackets : public OwnOperators
{
    private:
        BRACKET_METHOD bm;  //!< The bracket method
        Field3D DDXPhi;     //!< \f$\partial_\rho \phi\f$
        OwnBCs ownBC;       //!< Needed for setting the inner rho
        int ghostIndX;      //!< Index for first outer ghostpoint in x
    public:
        // Constructors
        OwnOpBasicBrackets(Options *options);

        //! Will throw error
        Field3D div_uE_dot_grad_n_GradPerp_phi(const Field3D &n,
                                               const Field3D &phi);
        //! Operator for \f$\{\phi, \Omega^D\}\f$
        Field3D vortDAdv (const Field3D &phi, const Field3D &vortD);
        /*! Operator for
         * \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$
         * (only used in 2Brackets)
         */
        Field3D kinEnAdvN(const Field3D &phi, const Field3D &n);

        //! Destructor
        /* NOTE: {} in the end is needed
         *       If else the compiler gives
         *       "udefined reference to `vtable for ...'"
         */
        virtual ~OwnOpBasicBrackets(){};
};

// Function bodies of the non-inlined functions are located in the .cxx file
#include "../src/ownOperators.cxx"

#endif
