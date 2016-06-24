// *************** Simulation of DDXCylinderHighResCenter *********************
/* Geometry
 *  x - The radial coordinate (rho)
 *  y - The height of the cylinder (z)
 *  z - The azimuthal coordinate (theta)
 */

#include "DDXDDXCylinderHighResCenter.hxx"

Field3D DDXCylHighResCenter(Field3D const &f);


// Set xout ghost point of the resulting field of DDX
// *****************************************************************************
void xout_after_DDX(Field3D &out, Field3D const &in){
    /* Info:
     * Set xout ghost point of the resulting field of DDX
     *
     * Input:
     * in    - The field that was used to find out
     * out   - DDX(in), but without the ghost point in xout
     *
     * Output:
     * out- DDi(in), with the ghost point of xout set
     */
    TRACE("Halt in xout_after_DDX");

    if (mesh->lastX()){
        /* NOTE: xend
         *       xend = index value of last inner point on this processor
         *       xend+1 = first guard point
         */
        int x_ind = mesh->xend + 1;

        // Newton polynomial of fourth order (including boundary) evaluated at ghost
        for(int y_ind = mesh->ystart; y_ind <= mesh->yend; y_ind++){
            for(int z_ind = 0; z_ind < mesh->ngz -1; z_ind ++){
                out(x_ind, y_ind, z_ind) =
                    - (1.0/5.0)*out(x_ind-3, y_ind, z_ind)
                    +           out(x_ind-2, y_ind, z_ind)
                    - 3.0*      out(x_ind-1, y_ind, z_ind)
                    + (16.0/5.0)*
                      // Calculation of boundary value
                      ((in(x_ind, y_ind, z_ind)-in(x_ind-1, y_ind, z_ind))/
                       mesh->dx(x_ind, y_ind));
            }
        }
    }
}
// *****************************************************************************

// Initialization of the physics
// ############################################################################
int DDXCylinderHighResCenter::init(bool restarting) {
    TRACE("Halt in DDXCylinderHighResCenter::init");

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Load from the geometry
    // ************************************************************************
    Options *geom = options->getSection("geom");
    geom->get("Lx", Lx, 0.0);
    // ************************************************************************

    // Obtain the fields
    // ************************************************************************
    // f
    f = FieldFactory::get()
        ->create3D("f:function", Options::getRoot(), mesh, CELL_CENTRE, 0);

    // S
    S = FieldFactory::get()
        ->create3D("S:solution", Options::getRoot(), mesh, CELL_CENTRE, 0);
    // ************************************************************************

    // Add a FieldGroup to communicate
    // ************************************************************************
    // Only these fields will be taken derivatives of
    com_group.add(f);
    // ************************************************************************

    // Set boundaries manually
    // ************************************************************************
    f.setBoundary("f");
    f.applyBoundary();
    ownBC.innerRhoCylinder(f);
    // ************************************************************************

    // Communicate before taking derivatives
    mesh->communicate(com_group);

    output << "\n\n\n\n\n\n\nNow running test" << std::endl;

    // Calculate
    DDXf = DDXCylHighResCenter(f);
    ownBC.innerRhoCylinder(DDXf); // Set inner rho boundary
    xout_after_DDX(DDXf, f);      // Set outer rho boundary
    mesh->communicate(DDXf);
    S_num = DDXCylHighResCenter(DDXf);

    // Error in S
    e = S_num - S;

    // Save the variables
    SAVE_ONCE (Lx);
    SAVE_ONCE3(f, S, S_num);
    SAVE_ONCE (e);

    // Finalize
    dump.write();
    dump.close();

    output << "\nFinished running test, now quitting\n\n\n\n\n\n" << std::endl;

    // Wait for all processors to write data
    MPI_Barrier(BoutComm::get());

    return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int DDXCylinderHighResCenter::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(DDXCylinderHighResCenter);

// Destructor
DDXCylinderHighResCenter::~DDXCylinderHighResCenter(){
    TRACE("DDXCylinderHighResCenter::~DDXCylinderHighResCenter");
}

// New implementation
Field3D DDXCylHighResCenter(Field3D const &f)
{
    TRACE("DDXCylHighResCenter");

    Field3D result;

    // Calculate points away from the center
    result = DDX(f);

    /* NOTE: The index corresponding to pi
     *       Since z in [0, 2 pi[, the z index corresponding to pi is
     *       (mesh->ngz -1) / 2, where mesh->ngz - 1 is the last actual z point
     *       (in addition there is one extra z point never used)
     */
    int piIndex = (mesh->ngz -1)/2;

    if(mesh->firstX()) {
        // For all z indices corresponding to a theta angle below pi
        int xInd = mesh->xstart;
        for (int yInd = mesh->ystart; yInd <= mesh->yend; yInd ++) {
            for (int zInd = 0; zInd < piIndex; zInd ++){
                // Use a five point stencil
                result(xInd, yInd, zInd) =
                  (
                      // Second "ghost"
                      (1.0/12.0)*f(xInd+1, yInd, zInd + piIndex)
                      // First ghost (set ownBC.innerRhoCylinder(f))
                    - (2.0/3.0 )*f(xInd-1, yInd, zInd          )
                    + (2.0/3.0 )*f(xInd+1, yInd, zInd          )
                    - (1.0/12.0)*f(xInd+2, yInd, zInd          )
                  )/mesh->dx(xInd, yInd)
                    ;

            }
        }
        /* NOTE: Addressing "off by one" for the z points
         *        We loop up over the rest of the z points. Note however that
         *        ngz is a number that starts counting on 1. Thus we need to
         *        subtract by one since we count arrays starting from 0.
         */
        // For all z indices corresponding to a theta value including and above
        // pi
        for (int yInd = mesh->ystart; yInd <= mesh->yend; yInd ++) {
            for (int zInd = piIndex; zInd < mesh->ngz -1; zInd ++){
                // Use a five point stencil
                result(xInd, yInd, zInd) =
                  (
                      // Second "ghost"
                      (1.0/12.0)*f(xInd+1, yInd, zInd - piIndex)
                      // First ghost (set ownBC.innerRhoCylinder(f))
                    - (2.0/3.0 )*f(xInd-1, yInd, zInd          )
                    + (2.0/3.0 )*f(xInd+1, yInd, zInd          )
                    - (1.0/12.0)*f(xInd+2, yInd, zInd          )
                  )/mesh->dx(xInd, yInd)
                    ;
            }
        }
    }

    return result;
}
