#ifndef __HELPERFUNCTIONS_H__
#define __HELPERFUNCTIONS_H__

#include <bout.hxx>

//! Function which returns the poloidal average of a field
Field3D const polAvg(Field3D const &f, int const &xInd, int const &yInd);

// Integrals
//! Volume integral
void volumeIntegral(Field3D const &f, BoutReal &result);

//! Surface integral
void surfaceEdgeIntegral(Vector3D const &f             ,
                         std::vector<BoutReal> &results);

#include "../src/helperFunctions.cxx"

#endif
