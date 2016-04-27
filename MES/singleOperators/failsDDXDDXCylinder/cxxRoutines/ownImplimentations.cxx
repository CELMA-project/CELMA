#include "../DDXDDXCylinder.hxx"

// Own operators
// ############################################################################
// Get partial_rho phi
// ****************************************************************************
Field3D DDXDDXCylinder::calcPartialRhoPhi(const Field3D &phi){
    /* Explanation:
     * Calculate partial_rho phi
     * Set the ghost points in rho
     *  - Inner rho BC from circle condition
     *  - Outer rho BC from calculating the value of the BC, then 4th order
     *    extrapolate
     * Set the ghost points in the parallel direction
     * (needed for parAdvGradPerpPhi_n)
     *  - yup and down from using corner points
     */
    TRACE("Halt in DDXDDXCylinder::calcPartialRhoPhi");

    Field3D result;
    // Calculate the derivative
    result = DDX(phi);
    // Set boundaries
    /* Set upper/lower parallel boundary
     * We need this as
     * 1. We are taking the parallel derivative in u_par*partial_par
     * 2. DDZ( ,true) also calculates the y boundaries
     */
    DDX_yup_and_ydown(phi, result);
    innerRhoCylinder(result); // Set inner rho boundary
    xout_after_DDX(result, phi);     // Set outer rho boundary

    return result;
}
// ****************************************************************************
// ############################################################################


// Auxiliary functions
// ############################################################################
// ****************************************************************************
void DDXDDXCylinder::DDX_yup_and_ydown(const Field3D &in, Field3D &out){
    /* Explanation
     * Taking x-derivative of the ghost points in y, by calling
     * DDX_one_xz_plane_without_BC.
     *
     * Input
     * in  - The field to take x-derivative of
     * out - The field to apply the ghost point on
     *
     * Output
     * out - The field with the new ghost point
     *
     * NOTE: Make sure that in is communicated before using this routine
     * NOTE: No communication is needed until we need to take y-derivatives
     */
    TRACE("Halt in DDXDDXCylinder::DDX_yup_and_ydown");

    if(mesh->firstY()){
        /* NOTE: ystart
         *       ystart = index value of first inner point on this processor
         *       ystart-1 = first guard point
         */
        int y_ind = mesh->ystart-1;
        DDX_one_xz_plane_without_BC(in, y_ind, out);
    }
    if(mesh->lastY()){
        /* NOTE: yend
         *       yend = index value of last inner point on this processor
         *       yend+1 = first guard point
         */
        int y_ind = mesh->yend + 1;
        DDX_one_xz_plane_without_BC(in, y_ind, out);
    }
}
// ****************************************************************************

void DDXDDXCylinder::DDX_one_xz_plane_without_BC(Field3D const &in,
                                            int const &y_ind,
                                            Field3D &out){
    /* Explanation
     * Taking x-derivative of a specific y slice. We will use a one sided
     * scheme in the points just before the first ghost point. We will not use
     * information of the boundaries.
     *
     * Input
     * in    - The field to take x-derivative of
     * y_ind - The y-index to take the derivative in
     * out   - The field to set the result of the x-derivative in
     *
     * Output
     * out - The field after setting the result of the x-derivative
     *
     * NOTE: Make sure that in is communicated before using this routine
     * NOTE: No communication is needed until we need to take y-derivatives
     *
     * Source:
     * https://en.wikipedia.org/wiki/Finite_difference_coefficient
     */
    TRACE("Halt in DDXDDXCylinder::DDX_one_xz_plane_without_BC");

    // Forward
    // One sided stencil 3rd order accurate
     for(int z_ind = 0; z_ind < mesh->ngz -1; z_ind ++){
         out(mesh->xstart, y_ind, z_ind) =
             ( - (11.0/6.0)*in(mesh->xstart  , y_ind, z_ind)
               +      (3.0)*in(mesh->xstart+1, y_ind, z_ind)
               -  (3.0/2.0)*in(mesh->xstart+2, y_ind, z_ind)
               +  (1.0/3.0)*in(mesh->xstart+3, y_ind, z_ind)
             )/mesh->dx(0,0)
             ;
     }
    // Backward
    // One sided stencil 3rd order accurate
    for(int z_ind = 0; z_ind < mesh->ngz -1; z_ind ++){
        out(mesh->xend, y_ind, z_ind) =
            ( + (11.0/6.0)*in(mesh->xend  , y_ind, z_ind)
              -      (3.0)*in(mesh->xend-1, y_ind, z_ind)
              +  (3.0/2.0)*in(mesh->xend-2, y_ind, z_ind)
              -  (1.0/3.0)*in(mesh->xend-3, y_ind, z_ind)
            )/mesh->dx(0,0)
            ;
    }
    /* Centered
     * Second order accurate
     * Points closest to x guard cells are not looped over as this is taken
     * care of in the one sided stencils
     */
    for(int x_ind =  mesh->xstart + 1;
        x_ind <= mesh->xend   - 1;
        x_ind++){
        for(int z_ind = 0; z_ind < mesh->ngz -1; z_ind ++){
            out(x_ind, y_ind, z_ind) =
                ( - (1.0/2.0)*in(x_ind-1, y_ind, z_ind)
                  + (1.0/2.0)*in(x_ind+1, y_ind, z_ind)
                )/mesh->dx(0,0)
                ;
        }
    }
}
// ############################################################################
