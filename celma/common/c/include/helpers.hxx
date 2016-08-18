#ifndef __HELPERS_H__
#define __HELPERS_H__

#include <bout.hxx>

/*!
 * \class Helpers
 *
 * \brief Class containing helper functions used in monitoring and
 *        post-processing
 *
 * \author Michael Løiten
 * \date 2016.08.18
 */
class Helpers
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

        int xInd; //!< x-looping index
        int yInd; //!< y-looping index
        int zInd; //!< z-looping index

        /* \brief Results to be collected by MPI_Allreduce.
         *
         * localResult[0] - integral over v on the xin surface
         * localResult[1] - integral over v on the xout surface
         * localResult[2] - integral over v on the ydown surface
         * localResult[3] - integral over v on the yup surface
         */
        std::vector<BoutReal> localResults;

    public:
        // Constructor with list initialization
        Helpers ();

        //! Function which returns the poloidal average of a field
        Field3D const polAvg(Field3D const &f);

        // Integrals
        //! Surface integral
        void surfaceEdgeIntegral(Vector3D const &f,
                                 std::vector<BoutReal> &results);

        //! Volume integral
        void volumeIntegral(Field3D const &f, BoutReal &result);
};

#include "../src/helpers.cxx"

#endif
