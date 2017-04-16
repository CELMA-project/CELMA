#ifndef __OWNOPERATORS_H__
#define __OWNOPERATORS_H__

#include "ownBCs.hxx"         // Gives inner rho boundaries
#include <bout.hxx>           // Includes all necessary classes and types
#include <bout/constants.hxx> // Gives PI and TWOPI

class OwnOperators;
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
class OwnOperators {
private:
  // Data members
  /*! If a warning is given rather than throwing an error if
   *  insufficient number of points is found
   */
  bool warnPoints;

protected:
  // Data members
  FieldGenerator *bndryFuncGen;
  Field2D J;    //!< The Jacobian
  Field2D J2;   //!< The Jacobian raised to power 2
  Field2D invJ; //!< The inverse of the Jacobian

  int xInd; //!< x-index
  int yInd; //!< y-index
  int zInd; //!< z-index
public:
  // Constructors
  OwnOperators();

  // Destructors
  // NOTE: New is called, so should destruct
  // NOTE: Needs to be virtual in order for the child classes to call destruct
  virtual ~OwnOperators(){};

  // Member functions
  //! Operator for \f$\nabla\cdot_(f \nabla_\perp g)\f$
  Field3D div_f_GradPerp_g(const Field3D &f, const Field3D &g);
  //! Operator for \f$\nabla_\perp f\f$ in cylinder geometry
  Vector3D Grad_perp(const Field3D &f);
  //! Operator for \f$\partial_z \mathbf{v}\f$ in cylinder geometry
  Vector3D DDY(const Vector3D &f);

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
  //! Operator for \f$\{\phi, \Omega^D\}\f$
  virtual Field3D vortDAdv(const Field3D &phi, const Field3D &vortD) = 0;
  /*! Operator for
   * \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$
   * (only used in 2Brackets)
   */
  virtual Field3D kinEnAdvN(const Field3D &phi, const Field3D &n) = 0;

  //! Factory which chooses child class
  static OwnOperators *createOperators(Options *options = NULL);
};

// OwnOpBasicBrackets

/*!
 * \class OwnOpBasicBrackets
 *
 * \brief Implementation of
 *        \f$\nabla\cdot_(\mathbf{u}_e \cdot \nabla[n\nabla_\perp \phi])\f$
 *        using 2 basic brackets
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
 * \note An alternative derivation is given in appendix B of
 *       P. Popovich, M. Umansky, T. A. Carter, and B. Friedman - Phys. Plasmas
 * 17, 102107 2010
 *
 * \author Michael Løiten
 * \date 2016.07.23
 */
class OwnOpBasicBrackets : public OwnOperators {
private:
  BRACKET_METHOD bm; //!< The bracket method
  Field3D DDXPhi;    //!< \f$\partial_\rho \phi\f$
  OwnBCs ownBC;      //!< Needed for setting the inner rho
  int ghostIndX;     //!< Index for first outer ghostpoint in x
public:
  // Constructors
  OwnOpBasicBrackets();

  //! Operator for \f$\{\phi, \Omega^D\}\f$
  Field3D vortDAdv(const Field3D &phi, const Field3D &vortD);
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
