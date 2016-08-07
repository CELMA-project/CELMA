#ifndef __HELPERFUNCTIONS_H__
#define __HELPERFUNCTIONS_H__

#include <bout.hxx>

/* NOTE: Using a template class
 *       As we do not know the name of the physics class in advance, we
 *       implement this class as a template. This is because we need the
 *       derived class' "model" in order to access data members stored in that
 *       class.
 *       See
 *       http://stackoverflow.com/questions/1474416/c-passing-a-class-as-a-parameter
 */
template<class PhysicsClass>
class ownMonitors
{
private:
    //! Declaration of model (needed for accessing members of the PhysicsClass)
    static PhysicsClass *model;

public:
    // Constructors
    /* NOTE: Using default constructor
     *       The model is created without arguments in solver.hxx
     *       Thus,
     *
     *       classB m_input;
     *       classA classA(classB input) : m_input(input){};
     *
     *       would not work, as the physicsmodel is creating the object
     *       without arguments (and is not using list initialization)
     *       One could therefore either make an own main, or use "creators"
     *       as constructors (as used here)
     *
     *       http://stackoverflow.com/questions/7761676/calling-constructor-of-a-class-member-in-constructor
     */

    //! Alternative to a constructor
    void create(PhysicsClass *model);

    // Monitors
    //! Energy monitor
    static int energyIntMon(Solver *solver, BoutReal simtime, int iter, int NOUT);
};

/* NOTE: Reason for definition
 *       See
 *       http://stackoverflow.com/questions/7092765/what-does-it-mean-to-have-an-undefined-reference-to-a-static-member
 */
//! The model definition
template<class PhysicsClass> PhysicsClass* ownMonitors<PhysicsClass>::model;

#include "../src/ownMonitors.cxx"

#endif
