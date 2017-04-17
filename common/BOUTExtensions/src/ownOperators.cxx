#ifndef __OWNOPERATORS_CXX__
#define __OWNOPERATORS_CXX__

#include "../include/ownOperators.hxx"
#include <fft.hxx>           //Includes the FFT
#include <interpolation.hxx> //Includes the interpolation

// OwnOperators

/*!
 * \brief Constructor
 *
 * Constructor which sets the private member data
 *
 * \param[in] section What section to get bndryFuncGen from.
 *                    Default is "phi"
 */
OwnOperators::OwnOperators() {
  TRACE("Halt in OwnOperators::OwnOperators");

  // Calculate the powers of the Jacobian
  // ************************************************************************
  J = mesh->coordinates()->J;
  J2 = mesh->coordinates()->J ^ (2.0);
  invJ = 1.0 / (mesh->coordinates()->J);
  // ************************************************************************
}

/*!
 * This function now works as a constructor of the child-classes of OwnFilters
 */
OwnOperators *OwnOperators::createOperators(Options *options) {
  TRACE("Halt in OwnOperators::createOperators");

  // The filter option is by defualt found in the ownFilter section
  if (options == NULL) {
    options = Options::getRoot()->getSection("ownOperators");
  }

  string type;
  options->get("type", type, "BasicBrackets");

  if (lowercase(type) == lowercase("BasicBrackets")) {
    output << "OwnOperators type set to 'BasicBrackets'" << std::endl;
    return new OwnOpBasicBrackets();
  } else {
    // Create a stream which we cast to a string
    std::ostringstream stream;
    stream << "OwnOperators '" << type << "' not implemented\n"
           << "Available operators:\n"
           << "BasicBrackets - Consistent implementation using brackets "
              "with basic Arakawa.\n";
    std::string str = stream.str();
    // Cast the stream to a const char in order to use it in BoutException
    const char *message = str.c_str();

    throw BoutException(message);
  }
}

/*!
 * Operator for \f$\nabla\cdot_(f \nabla_\perp g)\f$ in cylindrical geometry.
 * We have that
 *
 * \f{eqnarray}{
 *   \nabla\cdot(f\nabla_\perp g) = f\nabla_\perp^2g + \nabla
 *   f\cdot \nabla_\perp g = f\nabla_\perp^2g + \nabla_\perp f\cdot
 *   \nabla_\perp g
 * \f}
 *
 * The expression for the perpendicular Laplacian can be found in the
 * coordinates manual. Note that in cylinder coordinates
 *
 * \f{eqnarray}{
 *   G^x &=& \frac{1}{J}\\
 *   G^y &=& 0\\
 *   G^z &=& 0\\
 *   g^{zz} &=& \frac{1}{\rho^2}\\
 *   g^{yy}\partial_y^2
 *   - \frac{1}{J}\partial_y\left(\frac{J}{g^{yy}}\partial_y\right)
 *   &=& 0
 * \f}
 *
 * \param[in] f The f field
 * \param[in] g The g field
 *
 * \return result The result of the operation
 */
Field3D OwnOperators::div_f_GradPerp_g(const Field3D &f, const Field3D &g) {
  TRACE("Halt in OwnOperators::div_f_GradPerp_g");

  Field3D result;

  result = f * D2DX2(g) + (f / J) * DDX(g) + (f / J2) * D2DZ2(g) +
           DDX(f) * DDX(g) + (1 / J2) * DDZ(f) * DDZ(g);

  return result;
}

/*!
 * \f$\nabla_\perp f\f$, equivalent to Grad_perp in vecops.cxx, but in
 * cylindrical geometry.  This means that there the y-component of the vector
 * is set to 0 (as there are no off-diagonal elements). The advantage of
 * introducing this operator is:
 *
 *      1. No need for setting the y-boundaries on the operand
 *      2. Reduced calculation time
 *
 * \param[in] f The original field
 *
 * \return result The result of the operation written in a covariant basis
 */
Vector3D OwnOperators::Grad_perp(const Field3D &f) {
  TRACE("Halt in OwnOperators::own_Grad_perp");
  Vector3D result;

  result.x = DDX(f);
  result.y = 0.0;
  result.z = DDZ(f);

  result.covariant = true;

  return result;
}

/*!
 * \f$\partial_z \mathbf{v}\f$. We note that in cylindrical geometry, all the
 * metrics are independant of \f$z\f$. Hence, this operator only calls the DDY
 * for each vector component.
 *
 * \param[in] f The original vector field
 *
 * \warning Only siutable for cylindrical geometry
 *
 * \return result The result of the operation.
 */
Vector3D OwnOperators::DDY(const Vector3D &v) {
  TRACE("Halt in OwnOperators::DDY");
  Vector3D result;

  result.covariant = v.covariant;

  /* NOTE: DDY is the name of the member function, but we would like to call
   *       the function from the global scope. Using the scope operation
   *       without a namespace fixes this
   */
  result.x = ::DDY(v.x);
  result.y = ::DDY(v.y);
  result.z = ::DDY(v.z);

  return result;
}

// OwnOpBasicBrackets

/*!
 * \brief Constructor
 *
 * Constructor which calls parent constructor and sets the bracket method
 */
OwnOpBasicBrackets::OwnOpBasicBrackets() : OwnOperators() {
  TRACE("Halt in OwnOpBasicBrackets::OwnOpBasicBrackets");

  bm = BRACKET_ARAKAWA;
}

/*!
 * Calculates \f$\{\phi, \Omega^D\}\f$
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns result The result of the operation
 */
Field3D OwnOpBasicBrackets::vortDAdv(const Field3D &phi, const Field3D &vortD) {
  TRACE("Halt in OwnOpBasicBrackets::vortDAdv");

  return invJ * bracket(phi, vortD, bm);
}

/*!
 * Calculates \f$\frac{1}{J2}\{\mathbf{u}_E\cdot\mathbf{u}_E, n\} \f$
 *
 * \param[in] phi The potential
 * \param[in] n The density (not used here)
 *
 * \returns result The result of the operation
 */
Field3D OwnOpBasicBrackets::kinEnAdvN(const Field3D &phi, const Field3D &n) {
  TRACE("Halt in OwnOpBasicBrackets::kinEnAdvN");

  Field3D result;

  // Calculate the derivative of phi
  DDXPhi = DDX(phi);
  // Reset inner boundary
  ownBC.innerRhoCylinder(DDXPhi);
  // Reset outer boundary
  if (mesh->lastX()) {
    /* NOTE: xend
     *       xend = index value of last inner point on this processor
     *       xend+1 = first guard point
     */
    ghostIndX = mesh->xend + 1;
    // Newton polynomial of fourth order (including boundary) evaluated at ghost
    for (yInd = mesh->ystart; yInd <= mesh->yend; yInd++) {
      for (zInd = 0; zInd < mesh->ngz - 1; zInd++) {
        DDXPhi(ghostIndX, yInd, zInd) =
            -DDXPhi(ghostIndX - 4, yInd, zInd) +
            4.0 * DDXPhi(ghostIndX - 3, yInd, zInd) -
            6.0 * DDXPhi(ghostIndX - 2, yInd, zInd) +
            4.0 * DDXPhi(ghostIndX - 1, yInd, zInd);
      }
    }
  }

  // Communicate before taking new derivative
  mesh->communicate(DDXPhi);

  result =
      bracket(((DDXPhi) ^ (2.0)) + ((invJ * DDZ(phi, true)) ^ (2.0)), n, bm);

  // Multiply with B/2
  return 0.5 * invJ * result;
}

#endif
