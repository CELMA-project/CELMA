#ifndef __HELPERS_H__
#define __HELPERS_H__

#include <bout.hxx>

class PolAvg;
class VolumeIntegral;
class SurfaceIntegral;

/*!
 * \class Helpers
 *
 * \brief Base class containing helper functions used in monitoring and
 *        post-processing
 *
 * \author Michael Løiten
 * \date 2016.08.18
 */
class Helpers
{
    protected:
        int xInd; //!< x-looping index
        int yInd; //!< y-looping index
        int zInd; //!< z-looping index
};

/*!
 * \class PolAvg
 *
 * \brief Class which performs poloidal average
 *
 * \author Michael Løiten
 * \date 2016.08.19
 */
class PolAvg : private Helpers
{
    public:
        //! Function which returns the poloidal average of a field
        Field3D const polAvg(Field3D const &f);
};

/*!
 * \class Volume integral
 *
 * \brief Class containing volume integration
 *
 * \author Michael Løiten
 * \date 2016.08.19
 */
class VolumeIntegral : private Helpers
{
    public:
        //! Volume integral
        void volumeIntegral(Field3D const &f, BoutReal &result);
};

/*!
 * \class SurfaceIntegrals
 *
 * \brief Class containing helper functions used in monitoring and
 *        post-processing
 *
 * \author Michael Løiten
 * \date 2016.08.19
 */
class SurfaceIntegral : private Helpers
{
    private:
        Vector3D dSXin  ; //!< dS for inner x edge surface
        Vector3D dSXout ; //!< dS for outer x edge surface
        Vector3D dSYup  ; //!< dS for upper y edge surface
        Vector3D dSYdown; //!< dS for lower y edge surface

        Field3D vDotDSXin;   //!< v*dS for inner x edge surface
        Field3D vDotDSXout;  //!< v*dS for outer x edge surface
        Field3D vDotDSYup;   //!< v*dS for upper y edge surface
        Field3D vDotDSYdown; //!< v*dS for lower y edge surface

        /* \brief Results to be collected by MPI_Allreduce.
         *
         * localResult[0] - integral over v on the xin surface
         * localResult[1] - integral over v on the xout surface
         * localResult[2] - integral over v on the ydown surface
         * localResult[3] - integral over v on the yup surface
         */
        std::vector<BoutReal> localResults;

        bool useInner_; //!< If the inner surface should be used (pointing in -x)
        bool useOuter_; //!< If the inner surface should be used (pointing in +x)
        bool useLower_; //!< If the downer surface should be used (pointing in -y)
        bool useUpper_; //!< If the upper surface should be used (pointing in +y)

        int xIndInnerLoc; //!< Global xIndInner index on current processor
        int xIndOuterLoc; //!< Global xIndOuter index on current processor
        int yIndLowerLoc; //!< Global yIndLower index on current processor
        int yIndUpperLoc; //!< Global yIndUpper index on current processor

        int MXG; //!< Number of ghost points in x
        int MYG; //!< Number of ghost points in y

    public:
        // Constructor
        SurfaceIntegral();

        //! Set what surfaces which will be calculated
        void setSurfaces(bool useInner,
                         bool useOuter,
                         bool useLower,
                         bool useUpper);

        //! Surface integral
        void surfaceEdgeIntegral(Vector3D const &f,
                                 int const &xIndInner,
                                 int const &xIndOuter,
                                 int const &yIndLower,
                                 int const &yIndUpper,
                                 std::vector<BoutReal> &results);
};

#include "../src/helpers.cxx"

#endif
