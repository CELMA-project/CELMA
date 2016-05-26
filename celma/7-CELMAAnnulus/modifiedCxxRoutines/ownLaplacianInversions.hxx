#ifndef __OWNLAPLACIANINVERSIONS_H__
#define __OWNLAPLACIANINVERSIONS_H__

#include <bout.hxx>         // Includes all necessary classes and types
#include "../../common/c/include/ownOperators.hxx" // Must know about the operators class

/*!
 * \class OwnLaplacianInversions
 *
 * \brief Class which has member functions used to do a laplacian inversion
 *
 * This class is used to to do a laplacian inversion
 *
 * \author Michael Løiten
 * \date 2016.12.04
 */
class OwnLaplacianInversions
{
    private:
        // Data members
        BoutReal atol; //!< Absolute tolerance of the solver
        BoutReal rtol; //!< Relative tolerance of the solver
        int maxit;     //!< Max iteration before throwing an exception
        bool monitor;  //!< Enable/disable monitoring of the solver steps
        // Needed in NaulinSolver
        Field3D phiCur;         //!< \f$\phi\f$ for the old step
        Field3D phiNext;        //!< \f$\phi\f$ for the new step
        Vector3D GradPerp_ln_n; //!< \f$\nabla_\perp \ln(n)\f$
        int curCount;           //!< Counts number of steps in the iterator
        //! Absolute error between phiCur and phiNext in \f$L_\infty\f$ norm
        BoutReal EAbsLInf;
        //! Relative error between phiCur and phiNext in \f$L_\infty\f$ norm
        BoutReal ERelLInf;
        Laplacian *phiSolver; //!< Solver object for the FFT solver
        OwnOperators *ownOp;  //!< Object with own operators (used in create)

    public:
        // Constructors
        // Use a create function instead of constructor (makes it easier to
        // pass arguments

        // Destructors
        ~OwnLaplacianInversions(); //!< Destroys the laplace object

        // Member functions
        //! Alternative to a constructor
        void create(OwnOperators &opObj,
                    OwnBCs &BCObj,
                    const string &section = "phiSolver");
        //! The NaulinSolver
        Field3D NaulinSolver(const Vector3D &gradPerpLnN,
                             const Field3D &n,
                             const Field3D &vortD,
                             const Field3D &phiInit,
                             Field3D &vort);
};

// Function bodies of the non-inlined functions are located in the .cxx file
#include "ownLaplacianInversions.cxx"

#endif
