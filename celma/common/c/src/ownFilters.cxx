#ifndef __OWNFILTERS_CXX__
#define __OWNFILTERS_CXX__

#include "../include/ownFilters.hxx"

/*!
 * \brief Constructor
 *
 * Constructor which calculates lmabdaMin from kMax.
 * kMax is found from the Nyquist sampling theorem which basically states that
 * a given mode number k is properly resolved with 2*k points. In addition, in
 * order not to get aliasing from non-linear wave coupling, Orzsag 2/3 rule is
 * used.
 *
 * \note Do not confuse the mode number with the inverse wavelength. Instead,
 *       the relation \f$\lambda = \frac{C}{k}\f$ holds, where \f$C\f$ is the
 *       circumference
 */
OwnFilters::OwnFilters()
{
    TRACE("Halt in OwnFilters::ownFilters");

    ncz = mesh->ngz-1;

    // Get MXG
    int MXG;
    Options *options = Options::getRoot();
    options->get("MXG", MXG, 2);

    // GlobalNx includes the ghost points
    BoutReal outerRho = (mesh->GlobalNx - 2*MXG - 0.5)*mesh->dx(0,0);
    BoutReal outerCircumference = TWOPI*outerRho;

    /* FIXME: Write the Nyquist sampling theorem
     *        Found that can be expressed as
     *        dz = (1/(2k))*2*pi => nz = 2k
     *        but need to formalize it
     */
    // Calculate the kMax from the Nyquist sampling theorem
    int kMax = int(ncz/2.0);

    // Use the Orzag's 2/3 rule
    kMax = int (floor((2.0/3.0)*kMax));

    // Calculate the corresponding minimum wavelength
    lambdaMin = outerCircumference/kMax;

    // Allocate fourierArray
    fourierArray = new dcomplex[ncz/2 + 1];

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
const Field3D OwnFilters::radialLowPass(const Field3D &var)
{
    TRACE("Halt in OwnFilters::radialLowPass");

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

/*!
 * \brief Destructor
 *
 * Deallocates the fourierArray
 */
OwnFilters::~OwnFilters()
{
    TRACE("Halt in OwnFilters::~ownFilters");

    // Deallocate fourierArray
    delete[] fourierArray;
}

#endif
