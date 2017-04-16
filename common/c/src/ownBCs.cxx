#ifndef __OWNBCS_CXX__
#define __OWNBCS_CXX__

#include "../include/ownBCs.hxx"

/*!
 * \brief Constructor
 *
 * Constructor which sets the private member data
 *
 * \warning firstUpperYGhost and firstLowerGhost
 *          should only be used for the processors using the respective
 *          boundaries
 */
OwnBCs::OwnBCs() {
  TRACE("Halt in OwnBCs::OwnBCs");

  // Set the piIndex
  /* NOTE: The index corresponding to pi
   *       Since z in [0, 2 pi[, the z index corresponding to pi is
   *       (mesh->ngz -1) / 2, where mesh->ngz - 1 is the last actual z point
   *       (in addition there is one extra z point never used)
   */
  piIndex = (mesh->ngz - 1) / 2;

  /* NOTE: xend
   *       xend = index value of last inner point on this processor
   *       xend+1 = first guard point
   */
  firstOuterXGhost = mesh->xend + 1;
  /* NOTE: yend
   *       yend = index value of last inner point on this processor
   *       yend+1 = first guard point
   */
  firstUpperYGhost = mesh->yend + 1;
  /* NOTE: ystart
   *       ystart = index value of first inner point on this processor
   *       ystart-1 = first guard point
   */
  firstLowerYGhost = mesh->ystart - 1;

  /* Check that there are enough points
   * ngy is the size of local mesh including guard cells
   */
  Options *switchOptions = Options::getRoot()->getSection("switch");
  switchOptions->get("warnPoints", warnPoints, false);
  if (mesh->ngy - 2 * mesh->ystart < 4) {

    // Create a stream which we cast to a string
    std::ostringstream stream;
    stream << "Not enough inner points i the y-direction\n"
           << "The cauchy and uEPar BC needs 3 inner points in y\n"
           << "extrapolateYUp and extrapolateYDown needs 4 inner points "
           << "in y\n"
           << "Currently the number of inner points is "
           << mesh->ngy - 2 * mesh->ystart;

    if (warnPoints) {
      output << "\n\n!!! WARNING\n" << stream.str() << "\n\n" << std::endl;
    } else {
      std::string str = stream.str();
      // Cast the stream to a const char in order to use it in BoutException
      const char *message = str.c_str();
      throw BoutException(message);
    }
  }
}

// Member functions
/*!
 * Sets the inner x ghost points of a Field3D in cylindrical coordinates.
 *
 * \param[in] f The original field
 * \param[out] f The field after setting the inner boundary
 *
 * ## Idea:
 * We will set the boundary at the inner \f$\rho\f$, so that the value of the
 * n'th ghost point = value of the diametrically n'th inner point to the n'th
 * point closest to the boundary
 *
 * ## Explanation of the procedure:
 * The innermost ghost points in x (closest to the boundary) with a z value
 * lower than (and excluding) \f$\pi\f$ will be set to the value of the
 * innermost internal point (meaning all points excluding the ghost points)
 * which lies \f$\theta + \pi\f$ away. The next ghost point will be set to the
 * value of the second innermost internal point which lies \f$\theta + \pi\f$
 * away, and so on.
 *
 * This means that the index of the point diametrically opposite to the
 * current ghost point will be on the form
 *
 * ~~~{.cpp}
 * corresponding_xIndex = (2*mesh->xstart) - (current_ghost_xIndex + 1)
 * corresponding_zIndex = current_ghost_zIndex + (mesh->ngz -1)/2
 * ~~~
 *
 * where `x` maps to \f$\rho\f$, and z maps to \f$\theta\f$.
 *
 * When we are at a z index corresponding to a theta angle equal to \f$\pi\f$ or
 * above, we must connect the ghost point at the current z index with the
 * internal point laying z index corresponding to a theta angle \f$\pi\f$ behind
 * (since we have no index corresponding to an angle above \f$2\pi[\f$)
 *
 * This means that the index of the point diametrically opposite to the
 * current ghost point will be on the form
 *
 * ~~~{.cpp}
 * corresponding_xIndex = (2*mesh->xstart) - (current_ghost_xIndex + 1)
 * corresponding_zIndex = current_ghost_yIndex - (mesh->ngz -1)/2
 * ~~~
 *
 * \sa innerRhoCylinderLoop
 *
 * \note We are also setting the inner x-y corner points
 */
void OwnBCs::innerRhoCylinder(Field3D &f) {
  TRACE("Halt in OwnBCs::innerRhoCylinder");

  if (mesh->firstX()) {
    // Set the boundary for the inner y points
    innerRhoCylinderLoop(f, mesh->ystart, mesh->yend);
    // Do the same for the ghost points in y
    if (mesh->firstY()) {
      innerRhoCylinderLoop(f, 0, firstLowerYGhost);
    }
    if (mesh->lastY()) {
      // Note that ngy starts counting from 1
      innerRhoCylinderLoop(f, firstUpperYGhost, mesh->ngy - 1);
    }
  }
}

/*!
 * Extrapolates to the first upper ghost point using a 4th order Newton
 * polynomial
 *
 * \param[in] f The original field
 * \param[out] f The field after extrapolating to the first upper ghost point
 *
 * \sa extrapolateYGhost
 */
void OwnBCs::extrapolateXOutGhost(Field3D &f) {
  TRACE("Halt in OwnBCs::extrapolateXOutGhost");

  if (mesh->lastX()) {
    for (int yInd = mesh->ystart; yInd <= mesh->xend; yInd++) {
      for (int zInd = 0; zInd < mesh->ngz - 1; zInd++) {
        f(firstOuterXGhost, yInd, zInd) =
            4.0 * f(firstOuterXGhost - 1, yInd, zInd) -
            6.0 * f(firstOuterXGhost - 2, yInd, zInd) +
            4.0 * f(firstOuterXGhost - 3, yInd, zInd) -
            f(firstOuterXGhost - 4, yInd, zInd);
      }
    }
  }
}

/*!
 * Extrapolates to the first lower ghost point and the first upper ghost
 * point using a 4th order Newton polynomial
 *
 * \param[in] f The original field
 * \param[out] f The field after extrapolating to the first upper ghost point
 * and the first
 *
 * \sa extrapolateYUp
 * \sa extrapolateDown
 */
void OwnBCs::extrapolateYGhost(Field3D &f) {
  TRACE("Halt in OwnBCs::extrapolateYGhost");

  extrapolateYUp(f);
  extrapolateYDown(f);
}

/*!
 * Extrapolates to the first upper ghost point using a 4th order Newton
 * polynomial
 *
 * \param[in] f The original field
 * \param[out] f The field after extrapolating to the first upper ghost point
 *
 * \sa extrapolateYGhost
 */
void OwnBCs::extrapolateYUp(Field3D &f) {
  TRACE("Halt in OwnBCs::extrapolateYUp");

  if (mesh->lastY()) {
    for (int xInd = mesh->xstart; xInd <= mesh->xend; xInd++) {
      for (int zInd = 0; zInd < mesh->ngz - 1; zInd++) {
        f(xInd, firstUpperYGhost, zInd) =
            4.0 * f(xInd, firstUpperYGhost - 1, zInd) -
            6.0 * f(xInd, firstUpperYGhost - 2, zInd) +
            4.0 * f(xInd, firstUpperYGhost - 3, zInd) -
            f(xInd, firstUpperYGhost - 4, zInd);
      }
    }
  }
}

/*!
 * Extrapolates to the first lower ghost point using a 4th order Newton
 * polynomial
 *
 * \param[in] f The original field
 * \param[out] f The field after extrapolating to the first lower ghost point
 *
 * \sa extrapolateYGhost
 */
void OwnBCs::extrapolateYDown(Field3D &f) {
  TRACE("Halt in OwnBCs::extrapolateYDown");

  if (mesh->firstY()) {
    for (int xInd = mesh->xstart; xInd <= mesh->xend; xInd++) {
      for (int zInd = 0; zInd < mesh->ngz - 1; zInd++) {
        f(xInd, firstLowerYGhost, zInd) =
            4.0 * f(xInd, firstLowerYGhost + 1, zInd) -
            6.0 * f(xInd, firstLowerYGhost + 2, zInd) +
            4.0 * f(xInd, firstLowerYGhost + 3, zInd) -
            f(xInd, firstLowerYGhost + 4, zInd);
      }
    }
  }
}

/*!
 * This function will set the ghost point of uEPar according to
 * the sheath boundary condition
 *
 * \param[in] uEPar   The field to set the ghost point for
 * \param[in] phi     The current potential (must contain a valid yup ghost
 *                    point)
 * \param[in] Lambda  \$f\ln\left(\frac{\mu}{2\pi}\right)\$f
 * \param[in] phiRef  The reference potential compared to the ground
 *                    (0 by default)
 *
 * \param[out] uEPar The field after the ghost point has been set
 *
 * ## Explanation of the procedure:
 *
 * \f{eqnarray}{
 * u_{e, \|, B} = (c_s\exp(\Lambda-(\phi_{Ref} + \phi_B)/T_e))
 * \f}
 *
 * where
 *      * \f$c_s = \sqrt{\frac{T_e}{m_i}} = 1\f$ (\f$c_s\f$ is normalized)
 *      * \f$_B\f$ denotes a value at the boundary
 * This gives
 *
 * \f{eqnarray}{
 * u_{e, \|, B} = (\exp(\Lambda-(\phi_{Ref} + \phi_B)))
 * \f}
 *
 * We will use a 4th order boundary polynomial to interpolate to the value
 * at the boundary, as this resides half between grid points. For a field
 * \f$f\f$, this gives
 *
 * \f{eqnarray}{
 * f_B = (5f_{n} + 15f_{n-1} - 5f_{n-2} + f_{n-3})/16
 * \f}
 *
 * where \f$f_{n}\f$ denotes the first upper ghost point in y
 *
 * This gives
 *
 * \f{eqnarray}{
 * u_{e, \|,B} =
 * (5u_{e, \|, n} + 15 u_{e, \|, n-1} - 5 u_{e, \|, n-2} + u_{e, \|, n-3})/16
 * =
 * \exp(\Lambda
 *      - ((\phi_{Ref} + 5 \phi_{n} + 15 \phi_{n-1} - 5\phi_{n-2} + \phi_{n-3})
 *         /16)
 * )
 * \f}
 *
 * Which rearranged gives
 *
 * \f{eqnarray}{
 * u_{e, \|, n}
 * =
 * (16/5)\exp(\Lambda
 *      - ((\phi_{Ref} + 5 \phi_{n} + 15 \phi_{n-1} - 5\phi_{n-2} + \phi_{n-3})
 *         /16)
 * )
 * - 3 u_{e, \|, n-1} + u_{e, \|, n-2} - (1/5)u_{e, \|, n-3}
 * \f}
 *
 * \par Sources:
 * 1. Eq (26) in Loizu et al Pop 19-2012, using
 *      * \f$\theta_N = \theta_{\phi} = \theta_{Te} = 0 \f$ due to
 *        \f$\alpha = \pi/2\f$
 *      * \f$\Lambda = \ln(\mu/2\pi)\f$
 *      * \f$\eta_m = (\phi_{MPE} - \phi_W)/Te\f$
 *      * \f$\phi_{MPE} = \phi_{CSE}\f$ since we do not have a magnetic
 * presheath
 * 2. Eq (26) in Naulin et al PoP 15-2008
 * 3. Equation F.6 in Tiago's PhD 2007
 */
void OwnBCs::uEParSheath(Field3D &uEPar, const Field3D &phi,
                         const BoutReal &Lambda, const BoutReal &phiRef) {
  TRACE("Halt in OwnBCs::uEParSheath");

  if (mesh->lastY()) {
    for (int xInd = mesh->xstart; xInd <= mesh->xend; xInd++) {
      for (int zInd = 0; zInd < mesh->ngz - 1; zInd++) {
        uEPar(xInd, firstUpperYGhost, zInd) =
            exp(Lambda - (phiRef +
                          (+5.0 * phi(xInd, firstUpperYGhost, zInd) +
                           15.0 * phi(xInd, firstUpperYGhost - 1, zInd) -
                           5.0 * phi(xInd, firstUpperYGhost - 2, zInd) +
                           phi(xInd, firstUpperYGhost - 3, zInd)) /
                              16.0)) *
                (16.0 / 5.0) -
            3.0 * uEPar(xInd, firstUpperYGhost - 1, zInd) +
            uEPar(xInd, firstUpperYGhost - 2, zInd) -
            (1.0 / 5.0) * uEPar(xInd, firstUpperYGhost - 3, zInd);
      }
    }
  }
}

/*!
 * This function will set the ghost point of jPar according to
 * the sheath boundary condition.
 *
 * \param[in] jPar    The field to set the ghost point for
 * \param[in] uEPar   The parallel electron velocity
 * \param[in] uIPar   The parallel ion velocity (must contain a valid yup
 *                    ghost point)
 * \param[in] phi     The current potential (must contain a valid yup ghost)
 * \param[in] n       The density (must contain a valid yup ghost)
 * \param[in] Lambda  \$f\ln\left(\frac{\mu}{2\pi}\right)\$f
 * \param[in] phiRef  The reference potential compared to the ground
 *                    (0 by default)
 *
 * \param[out] jPar The field after the ghost point has been set
 *
 * \note Although we know \f$u_{i, \|, B}\f$ exact, we are here setting the
 *       ghost point, so we need the value at the ghost point as an input
 * \note The value of \f$u_{e, \|, B}\f$ is calculated in the code
 * \note The ghost point of \f$n\f$, \f$\phi\f$ and \f$u_{i, \|, B}\f$ must be
 *       set before calling this function, whereas the value of \f$u_{e, \|}\f$
 *       is calculated within the function.
 *
 * ## Explanation of the procedure:
 *
 * We have that
 *
 * \f{eqnarray}{
 * u_{e, \|, B} &= (c_s\exp(\Lambda-(\phi_{Ref} + \phi_B)/T_e))\\
 * u_{i, \|, B} &= c_s
 * \f}
 *
 * where
 *      * \f$c_s = \sqrt{\frac{T_e}{m_i}} = 1\f$ (\f$c_s\f$ is normalized)
 *      * \f$_B\f$ denotes a value at the boundary
 * This gives
 *
 * \f{eqnarray}{
 * j_B =
 * n_B(u_{i, \|, B} - u_{e, \|, B}) =
 * n_B(1-\exp(\Lambda-(\phi_{Ref} + \phi_B)))
 * \f}
 *
 * We will use a 4th order boundary polynomial to interpolate \f$u_{e,\|}\f$
 * and \f$\phi\f$ to the value at the boundary, as this resides half between
 * grid points. For a field \f$f\f$, this gives
 *
 * \f{eqnarray}{
 * f_B = (5f_{k} + 15f_{k-1} - 5f_{k-2} + f_{k-3})/16
 * \f}
 *
 * where \f$f_{k}\f$ denotes the value at the first upper ghost point in y
 *
 * This gives
 *
 * \f{eqnarray}{
 * u_{e, \|, B} = (5u_{e, \|, k} + 15 u_{e, \|, k-1} - 5 u_{e, \|, k-2} + u_{e,
 * \|, k-3})/16
 * =
 * \exp(\Lambda
 *      - ((\phi_{Ref} + 5 \phi_{k} + 15 \phi_{k-1} - 5\phi_{k-2} + \phi_{k-3})
 *         /16)
 * )p
 * \f}
 *
 * Which rearranged gives
 *
 * \f{eqnarray}{
 * u_{e, \|, k}
 * =
 * (16/5)\exp(\Lambda
 *      - ((\phi_{Ref} + 5 \phi_{k} + 15 \phi_{k-1} - 5\phi_{k-2} + \phi_{k-3})
 *         /16)
 * )p
 * - 3 u_{e, \|, k-1} + u_{e, \|, k-2} - (1/5)u_{e, \|, k-3}
 * \f}
 *
 * and inserted in the original equation yields
 *
 * \f{eqnarray}{
 * j_{\|, k}
 * =
 * n_k(
 * u_{i, \|, k}
 * -
 * (
 * (16/5)\exp(\Lambda
 *      - ((\phi_{Ref} + 5 \phi_{k} + 15 \phi_{k-1} - 5\phi_{k-2} + \phi_{k-3})
 *         /16)
 * )p
 * - 3 u_{e, \|, k-1} + u_{e, \|, k-2} - (1/5)u_{e, \|, k-3}
 * )
 * \f}
 *
 *
 * \par Sources:
 * 1. Eq (26) in Loizu et al Pop 19-2012, using
 *      * \f$\theta_N = \theta_{\phi} = \theta_{Te} = 0 \f$ due to
 *        \f$\alpha = \pi/2\f$
 *      * \f$\Lambda = \ln(\mu/2\pi)\f$
 *      * \f$\eta_m = (\phi_{MPE} - \phi_W)/Te\f$
 *      * \f$\phi_{MPE} = \phi_{CSE}\f$ since we do not have a magnetic
 * presheath
 * 2. Eq (26) in Naulin et al PoP 15-2008
 * 3. Equation F.6 in Tiago's PhD 2007
 */
void OwnBCs::jParSheath(Field3D &jPar, const Field3D &uEPar,
                        const Field3D &uIPar, const Field3D &phi,
                        const Field3D &n, const BoutReal &Lambda,
                        const BoutReal &phiRef) {
  TRACE("Halt in OwnBCs::jParSheath");

  if (mesh->lastY()) {
    for (int xInd = mesh->xstart; xInd <= mesh->xend; xInd++) {
      for (int zInd = 0; zInd < mesh->ngz - 1; zInd++) {
        jPar(xInd, firstUpperYGhost, zInd) =
            n(xInd, firstUpperYGhost, zInd) *
            (uIPar(xInd, firstUpperYGhost, zInd) -
             (exp(Lambda - (phiRef +
                            (+5.0 * phi(xInd, firstUpperYGhost, zInd) +
                             15.0 * phi(xInd, firstUpperYGhost - 1, zInd) -
                             5.0 * phi(xInd, firstUpperYGhost - 2, zInd) +
                             phi(xInd, firstUpperYGhost - 3, zInd)) /
                                16.0)) *
                  (16.0 / 5.0) -
              3.0 * uEPar(xInd, firstUpperYGhost - 1, zInd) +
              uEPar(xInd, firstUpperYGhost - 2, zInd) -
              (1.0 / 5.0) * uEPar(xInd, firstUpperYGhost - 3, zInd)));
      }
    }
  }
}

/*!
 * This function will set the ghost point of momDensPar according to the sheath
 * boundary condition.
 *
 * \param[in] momDensPar The field to set the ghost point for
 * \param[in] uIPar      The parallel ion velocity
 * \param[in] n          The density
 *
 * \param[out] momDensPar The field after the ghost point has been set
 *
 * \note Although we know \f$u_{i, \|, B}\f$ exact, we are here setting the
 *       ghost point, so we need the value at the ghost point as an input
 *
 * ## Explanation of the procedure:
 * As we have sat BC's on both \f$u_{i, \|, B}\f$ and \f$n\f$, we simply need
 * to multiply them together at the ghost point
 */
void OwnBCs::parDensMomSheath(Field3D &momDensPar, const Field3D &uIPar,
                              const Field3D &n) {
  TRACE("Halt in OwnBCs::parDensMomSheath");

  if (mesh->lastY()) {
    for (int xInd = mesh->xstart; xInd <= mesh->xend; xInd++) {
      for (int zInd = 0; zInd < mesh->ngz - 1; zInd++) {
        momDensPar(xInd, firstUpperYGhost, zInd) =
            n(xInd, firstUpperYGhost, zInd) *
            uIPar(xInd, firstUpperYGhost, zInd);
      }
    }
  }
}

/*!
 * This function will set the ghost point of f according to
 * the Cauchy boundary condition.
 *
 * \param[in] f The field to set the ghost point for
 * \param[in] a The value at the boundary \$f f(0)=a \$f
 * \param[in] b The derivative at the boundary \$f \partial_z f\big|_0=b \$f
 * \param[in] t The current simulation time
 * \param[in] yUpExtrapolate True if we want to extrapolate to yUp
 *
 * \param[out] f The field after the ghost point has been set
 *
 * ## Explanation of the procedure:
 * By 4th order extrapolation \$f f(0)=a \$f, solving it for the first inner
 * point, substituting it into the 2nd order approximation of the derivative
 * and solve it for the ghost point, we find that
 *
 * \f{eqnarray}{
 * f_{g} = \frac{4 a}{5} - \frac{3 b}{4} h + \frac{f_{1}}{4} - \frac{f_{2}}{20}
 * \f}
 *
 * \sa extrapolateYUp
 * \warning This stencil is only 1st order accurate
 */
void OwnBCs::cauchyYDown(Field3D &f, const BoutReal &t,
                         bool const &yUpExtrapolate) {
  TRACE("Halt in OwnBCs::cauchyYDown");

  if (mesh->firstY()) {
    for (int xInd = mesh->xstart; xInd <= mesh->xend; xInd++) {
      // Calculate the x for the current xIndex
      /* NOTE:
       * For a local index, globalX returns the global value between 0
       * and 1 corresponding to that index.
       * When evaluating the function given in the function generator,
       * the x variable need to be in range 0-1 (even if it contains
       * expressions like geom:xl) in order to be consistent with the
       * rest of the code.
       */
      x = mesh->GlobalX(xInd);
      for (int zInd = 0; zInd < mesh->ngz - 1; zInd++) {
        // Calculating the z value
        z = TWOPI * zInd / (mesh->ngz - 1);
        // Calculate a and b
        a = aBndryFuncGen->generate(x, yValAtYDownBndry, z, t);
        b = bBndryFuncGen->generate(x, yValAtYDownBndry, z, t);
        // Set the ghost point
        f(xInd, firstLowerYGhost, zInd) =
            (4.0 / 5.0) * a -
            (3.0 / 4.0) * b * mesh->dy(xInd, firstLowerYGhost) +
            (1.0 / 4.0) * f(xInd, firstLowerYGhost + 2, zInd) -
            (1.0 / 20.0) * f(xInd, firstLowerYGhost + 3, zInd);
      }
    }
  }

  // Extrapolate if true
  if (yUpExtrapolate) {
    extrapolateYUp(f);
  }
}

/*!
 * Prepares the Cauchy boundary operator by calling
 * getAFunction, getBFunction and getYValAtYDownBndry
 *
 * \param[in] section The section to find 'a' and 'b' in
 *
 * \param[out] aBndryFuncGen The function generator for 'a'
 * \param[out] bBndryFuncGen The function generator for 'b'
 * \param[out] yValAtYDownBndry The value of the y coordinate at the boundary
 *
 * \sa getAFunction
 * \sa getBFunction
 * \sa getYValAtYDownBndry
 */
void OwnBCs::prepareCauchy(const string &section) {
  TRACE("Halt in OwnBCs::prepareCauchy");

  getAFunction(section);
  getBFunction(section);
  getYValAtYDownBndry();
}

/*!
 * Reads the function 'a' from the input file, and make a function generator
 * out of it
 *
 * \param[in] section The section to find 'a' in
 *
 * \param[out] aBndryFuncGen The function generator for 'a'
 *
 * \sa getBFunction
 */
void OwnBCs::getAFunction(const string &section) {
  TRACE("Halt in OwnBCs::getAFunction");
  // Get the function
  Options *varOptions = Options::getRoot()->getSection(section);
  string bndryFuncString;
  // Last argument in get is the default
  varOptions->get("a", bndryFuncString, "");
  if (bndryFuncString == "") {
    // Create a stream which we cast to a string
    std::ostringstream stream;
    stream << "'a' not found in section '" << section << "' "
           << "but is needed when setting the Cauchy BC\n";
    std::string str = stream.str();
    // Cast the stream to a const char in order to use it in BoutException
    const char *message = str.c_str();

    throw BoutException(message);
  }
  aBndryFuncGen = FieldFactory::get()->parse(bndryFuncString);
}

/*!
 * Reads the function 'b' from the input file, and make a function generator
 * out of it
 *
 * \param[in] section The section to find 'b' in
 *
 * \param[out] aBndryFuncGen The function generator for 'b'
 *
 * \sa getAFunction
 */
void OwnBCs::getBFunction(const string &section) {
  TRACE("Halt in OwnBCs::getBFunction");
  // Get the function
  Options *varOptions = Options::getRoot()->getSection(section);
  string bndryFuncString;
  // Last argument in get is the default
  varOptions->get("b", bndryFuncString, "");
  if (bndryFuncString == "") {
    // Create a stream which we cast to a string
    std::ostringstream stream;
    stream << "'b' not found in section '" << section << "' "
           << "but is needed when setting the Cauchy BC";
    std::string str = stream.str();
    // Cast the stream to a const char in order to use it in BoutException
    const char *message = str.c_str();

    throw BoutException(message);
  }
  bBndryFuncGen = FieldFactory::get()->parse(bndryFuncString);
}

/*!
 * Calculates the value of the y coordinate at the boundary
 *
 * \param[out] yValAtYDownBndry The value of the y coordinate at the boundary
 */
void OwnBCs::getYValAtYDownBndry() {
  TRACE("Halt in OwnBCs::getYValAtYDownBndry");
  // Get the first y value
  if (mesh->firstY()) {
    /* NOTE:
     * For a local index, globalY returns the global value between 0 and 1
     * corresponding to that index.
     * When evaluating the function given in a function generator, the y
     * variable need to be in range 0-2*pi (even if it contains expressions
     * like geom:yl) in order to be consistent with the rest of the code.
     */
    yValAtYDownBndry =
        PI * // 0.5*TWOPI = PI
        (mesh->GlobalY(mesh->ystart - 1) + mesh->GlobalY(mesh->ystart));
  }
}

// Auxiliary
/*!
 * The actual loops which sets the ghost points for inner \f$\rho\f$
 *
 * \param[in] f The field to loop over
 * \param[in] yStart The start (inclusive) of the y index we are looping over
 * \param[in] yEnd The end (inclusive) of the y index we are looping over
 *
 * \param[out] f The field which has the inner boundaries set
 *
 * \sa innerRhoCylinder
 */
void OwnBCs::innerRhoCylinderLoop(Field3D &f, const int &yStart,
                                  const int &yEnd) {
  /* NOTE: Addressing "off by one" for the inner ghost points in x
   *       The first point on the current processor is 0. Since we are doing
   *       this procedure for all the ghost points, we must start from zero.
   *       We will loop all the way up to (but excluding) the first inner
   *       point. This can be achieved by using < instead of <=
   */
  /* NOTE: Addressing "off by one" for the y index
   *       We want to loop over all the inner points in y. Thus, we start
   *       on mesh->ystart, and loop until (and including) mesh->yend. We
   *       do this by using <= instead of < in the loop
   */
  /* NOTE: No need for old_field and new_field
   *       Only ghost points are dealt with on the LHS
   *       Only inner points are dealt with on the RHS
   */
  /* NOTE: Addressing "off by one" for the z points
   *       We loop up to (but not including) the pi index
   */
  TRACE("Halt in OwnBCs::innerRhoCylinderLoop Field3D");

  // For all z indices corresponding to a theta angle below pi
  for (int xInd = 0; xInd < mesh->xstart; xInd++) {
    for (int yInd = yStart; yInd <= yEnd; yInd++) {
      for (int zInd = 0; zInd < piIndex; zInd++) {
        // Set the value on the ghost point
        f(xInd, yInd, zInd) =
            f(2 * mesh->xstart - (xInd + 1), yInd, zInd + piIndex);
      }
    }
  }
  /* NOTE: Addressing "off by one" for the z points
   *        We loop up over the rest of the z points. Note however that
   *        ngz is a number that starts counting on 1. Thus we need to
   *        subtract by one since we count arrays starting from 0.
   */
  // For all z indices corresponding to a theta value including and above
  // pi
  for (int xInd = 0; xInd < mesh->xstart; xInd++) {
    for (int yInd = yStart; yInd <= yEnd; yInd++) {
      for (int zInd = piIndex; zInd < mesh->ngz - 1; zInd++) {
        // Set the value on the ghost point
        f(xInd, yInd, zInd) =
            f(2 * mesh->xstart - (xInd + 1), yInd, zInd - piIndex);
      }
    }
  }
}

#endif
