#ifndef __HELPERFUNCTIONS_H__
#define __HELPERFUNCTIONS_H__

#include <bout.hxx>

//! Function which returns the poloidal average of a field
Field3D const polAvg(Field3D const &f, int const &xInd, int const &yInd);

// Monitors
//! Energy monitor
int energyIntMon(Solver *solver, BoutReal simtime, int iter, int NOUT)

#include "../src/helperFunctions.cxx"

#endif
