#ifndef COMMON_BOUTEXTENSIONS_INCLUDE_SIGNEDJACOBIAN_HXX_
#define COMMON_BOUTEXTENSIONS_INCLUDE_SIGNEDJACOBIAN_HXX_

#include <bout/coordinates.hxx>
#include <bout/mesh.hxx>
#include <field_factory.hxx>
#include <options.hxx>

/*!
 * \brief Restore the signed Jacobian J = xl of the cylinder
 *
 * The radial coordinate x runs through the axis, so the Jacobian J = xl is
 * negative in the guard cell across the axis. BOUT++ 5 refuses to start with
 * J <= 0, so the input file gives it J = abs(xl). This function puts the
 * signed value back, as BOUT++ 4 had it, and recalculates every metric term
 * derived from J (the connection coefficients) with it.
 *
 * Call it first in init(), before any operator is used.
 */
inline void restoreSignedJacobian() {
  Coordinates *coord = bout::globals::mesh->getCoordinates();
  coord->J = FieldFactory::get()->create2D("geom:xl", Options::getRoot(),
                                           bout::globals::mesh, CELL_CENTRE, 0);
  coord->geometry();
}

#endif // COMMON_BOUTEXTENSIONS_INCLUDE_SIGNEDJACOBIAN_HXX_
