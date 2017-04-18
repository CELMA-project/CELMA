#ifndef __OWNFILTERS_CXX__
#define __OWNFILTERS_CXX__

#include "../include/ownFilters.hxx"

// OwnFilters

/*!
 * \brief Constructor
 *
 * Constructor which calculates kMax.
 * kMax is found from the Nyquist sampling theorem which basically states that
 * a given mode number k is properly resolved with 2*k points. In addition, in
 * order not to get aliasing from non-linear wave coupling, Orzsag 2/3 rule is
 * used.
 *
 * \param[in] *options Pointer to options parser
 */
OwnFilters::OwnFilters(Options *options) {
  TRACE("Halt in OwnFilters::OwnFilters");

  ncz = mesh->LocalNz - 1;

  // Calculate the kMax from the Nyquist sampling theorem
  kMax = int(ncz / 2.0);

  // Use the Orzag's 2/3 rule
  /* NOTE: Floor selects the closest lower number.
   *       This ensures that we do not get one spurious mode.
   */
  kMax = int(floor((2.0 / 3.0) * kMax));

  // Allocate fourierArray
  fourierArray = new dcomplex[ncz / 2 + 1];
}

/*!
 * This function now works as a constructor of the child-classes of OwnFilters
 */
OwnFilters *OwnFilters::createFilter(Options *options) {
  TRACE("Halt in OwnFiltersFactory::newOwnFilter");

  // The filter option is by defualt found in the ownFilter section
  if (options == NULL) {
    options = Options::getRoot()->getSection("ownFilters");
  }

  string type;
  options->get("type", type, "none");

  if (lowercase(type) == lowercase("none")) {
    output << "Filter type set to 'none'" << std::endl;
    return new OwnFiltAllPass(options);
  } else if (lowercase(type) == lowercase("radialLowPass")) {
    output << "Filter type set to 'radialLowPass'" << std::endl;
    return new OwnFiltRadialLowPass(options);
  } else {
    // Create a stream which we cast to a string
    std::ostringstream stream;
    stream << "Filtertype '" << type << "' not implemented\n"
           << "Available filters:\n"
           << "none          - No filtering will be performed\n"
           << "radialLowPass - Filtering dependant on rho\n";
    std::string str = stream.str();
    // Cast the stream to a const char in order to use it in BoutException
    const char *message = str.c_str();

    throw BoutException(message);
  }
}

/*!
 * \brief Destructor
 *
 * Deallocates the fourierArray
 */
OwnFilters::~OwnFilters() {
  TRACE("Halt in OwnFilters::~OwnFilters");

  // Deallocate fourierArray
  delete[] fourierArray;
}

// OwnFiltAllPass

/*!
 * Filter modified from the standard LowPass filter in BOUT++.
 * The shortest wavelength poloidally at any \f$\rho\f$ is being set from the
 * Nyquist sampling theorem and Orszag's 2/3 rule at the outermost radius.
 *
 * The corresponding mode number to the shortest wavelength is calculated for a
 * particular \f$\rho\f$, and is used as the maximum k at that \f$\rho\f$
 *
 * \param[in] var  Variable to be filtered.
 *
 * \returns result The filtered variable
 */
const Field3D OwnFiltAllPass::ownFilter(const Field3D &var) {
  TRACE("Halt in OwnFiltAllPass::ownFilter");

  return var;
}

// OwnFiltRadialLowPass

/*!
 * \brief Constructor
 *
 * Constructor which calculates lmabdaMin from kMax.
 * The list initializer calls the parent constructor with an argument.
 *
 * \note Do not confuse the mode number with the inverse wavelength. Instead,
 *       the relation \f$\lambda = \frac{C}{k}\f$ holds, where \f$C\f$ is the
 *       circumference
 * \warning Only works on grid equidistant in \f$\rho\f$
 */
OwnFiltRadialLowPass::OwnFiltRadialLowPass(Options *options)
    : OwnFilters(options) {
  TRACE("Halt in OwnFiltRadialLowPass::OwnFiltRadialLowPass");

  // Get MXG
  int MXG;
  options = Options::getRoot();
  options->get("MXG", MXG, 2);

  /* NOTE: outerRho
   * GlobalNx includes the ghost points and counts from one
   * The length made up by the inner points is therefore
   * (mesh->GlobalNx - 2*MXG - 1)*dx
   * However, we also need the line segment from the boundary (in our case
   * the origin of the cylinder) to the first inner point, located dx/2 away.
   * Thus:
   */
  BoutReal outerRho = (mesh->GlobalNx - 2 * MXG - 0.5) * mesh->coordinates()->dx(0, 0);
  BoutReal outerCircumference = TWOPI * outerRho;

  // Calculate the corresponding minimum resolvable wavelength
  lambdaMin = outerCircumference / kMax;

  // Whether or not to throw exception if the only mode present in the center is
  // the offset mode
  bool throwWarning;
  options->get("throwWarning", throwWarning, true);

  circumference =
      TWOPI * mesh->coordinates()->J(mesh->xstart, 0); // Lowest circumference (inner point)
  kMaxCurrent = int(floor(circumference / lambdaMin));
  if (kMaxCurrent <= 0) {
    if (throwWarning) {
      std::cout << "WARNING: At inner rho, kMax = circumference/lambdaMin = "
                << circumference / lambdaMin
                << ", so floor(circumference/lambdaMin) <= 0" << std::endl;
    } else {
      std::cout << "At inner rho circumference/lambdaMin = "
                << circumference / lambdaMin << std::endl;
      throw BoutException("kMax at inner rho is equal or below 0."
                          "\nFix it by increasing nz, decreasing nx or "
                          "decreasing Lx");
    }
  }
}

/*!
 * Filter modified from the standard LowPass filter in BOUT++.
 * The shortest wavelength poloidally at any \f$\rho\f$ is being set from the
 * Nyquist sampling theorem and Orszag's 2/3 rule at the outermost radius.
 *
 * The corresponding mode number to the shortest wavelength is calculated for a
 * particular \f$\rho\f$, and is used as the maximum k at that \f$\rho\f$
 *
 * \param[in] var  Variable to be filtered.
 *
 * \returns result The filtered variable
 */
const Field3D OwnFiltRadialLowPass::ownFilter(const Field3D &var) {
  TRACE("Halt in OwnFiltRadialLowPass::ownFilter");

  Field3D result;

  if (!var.isAllocated()) {
    return var;
  }

  result.allocate();

  for (int xInd = 0; xInd < mesh->LocalNx; xInd++) {
    // Set the current kMax (J = rho in cylinder coordinates)
    circumference = TWOPI * mesh->coordinates()->J(xInd, 0);
    // Abs since the inner ghost point of the Jacobian can be negative
    kMaxCurrent = int(abs(floor(circumference / lambdaMin)));
    for (int yInd = 0; yInd < mesh->LocalNy; yInd++) {
      // Take the FFT for a given radius at a given parallel plane
      rfft(&(var(xInd, yInd, 0)), ncz, fourierArray);

      // Filter in z
      for (int zInd = kMaxCurrent + 1; zInd <= ncz / 2; zInd++) {
        // NOTE: This also sets the imaginary part to 0
        fourierArray[zInd] = 0.0;
      }

      // Reverse FFT
      irfft(fourierArray, ncz, &(result(xInd, yInd, 0)));
      result(xInd, yInd, ncz) = result(xInd, yInd, 0);
    }
  }

  result.setLocation(var.getLocation());

  return result;
}

#endif
