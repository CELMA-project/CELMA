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
OwnFilters::OwnFilters(Options *options)
{
    TRACE("Halt in OwnFilters::OwnFilters");

    ncz = mesh->ngz-1;

    /* FIXME: Write the Nyquist sampling theorem
     *        Found that can be expressed as
     *        dz = (1/(2k))*2*pi => nz = 2k
     *        but need to formalize it
     */
    // Calculate the kMax from the Nyquist sampling theorem
    kMax = int(ncz/2.0);

    // Use the Orzag's 2/3 rule
    kMax = int (floor((2.0/3.0)*kMax));

    // Allocate fourierArray
    fourierArray = new dcomplex[ncz/2 + 1];
}

/*!
 * This function now works as a constructor of the child-classes of OwnFilters
 */
OwnFilters* OwnFilters::createFilter(Options *options)
{
    TRACE("Halt in OwnFiltersFactory::newOwnFilter");

    // The filter option is by defualt found in the ownFilter section
    if(options == NULL){
        options = Options::getRoot()->getSection("ownFilters");
    }

    string type;
    options->get("type", type, "none");

    if(lowercase(type) == lowercase("none")){
        output << "Filter type set to 'none'" << std::endl;
        return new OwnFiltAllPass(options);
    }
    else if(lowercase(type) == lowercase("radialLowPass")){
        output << "Filter type set to 'radialLowPass'" << std::endl;
        return new OwnFiltRadialLowPass(options);
    }
    else {
        // Create a stream which we cast to a string
        std::ostringstream stream;
        stream << "Filtertype '"<< type << "' not implemented\n"
               << "Available filters:\n"
               << "none          - No filtering will be performed\n"
               << "radialLowPass - Filtering dependant on rho\n"
               ;
        std::string str =  stream.str();
        // Cast the stream to a const char in order to use it in BoutException
        const char* message = str.c_str();

        throw BoutException(message);
    }
}

/*!
 * \brief Destructor
 *
 * Deallocates the fourierArray
 */
OwnFilters::~OwnFilters()
{
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
const Field3D OwnFiltAllPass::ownFilter(const Field3D &var)
{
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
 */
OwnFiltRadialLowPass::OwnFiltRadialLowPass(Options *options) : OwnFilters(options)
{
    TRACE("Halt in OwnFiltRadialLowPass::OwnFiltRadialLowPass");

    // Get MXG
    int MXG;
    options = Options::getRoot();
    options->get("MXG", MXG, 2);

    // GlobalNx includes the ghost points
    BoutReal outerRho = (mesh->GlobalNx - 2*MXG - 0.5)*mesh->dx(0,0);
    BoutReal outerCircumference = TWOPI*outerRho;

    // Calculate the corresponding minimum wavelength
    lambdaMin = outerCircumference/kMax;

    // Throw exception if the only mode present in the center is the offset mode
    circumference = TWOPI*mesh->J(1, 0);  // Lowest circumference
    kMaxCurrent   = int(floor(circumference/lambdaMin));
    if(kMaxCurrent <= 0){
        throw BoutException("kMaxCurrent at inner rho is equal or below 0.\n"
                            "Fix by increasing nz");
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
const Field3D OwnFiltRadialLowPass::ownFilter(const Field3D &var)
{
    TRACE("Halt in OwnFiltRadialLowPass::ownFilter");

    Field3D result;

    if(!var.isAllocated()){
        return var;
    }

    result.allocate();

    for(int xInd=0; xInd<mesh->ngx; xInd++) {
        // Set the current kMax (J = rho in cylinder coordinates)
        circumference = TWOPI*mesh->J(xInd, 0);
        kMaxCurrent   = int(abs(floor(circumference/lambdaMin)));
        for(int yInd=0; yInd<mesh->ngy; yInd++) {
            // Take the FFT for a given radius at a given parallel plane
            rfft(&(var(xInd, yInd, 0)), ncz, fourierArray);

            // Filter in z
            for(int zInd=kMaxCurrent+1; zInd<=ncz/2; zInd++){
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
