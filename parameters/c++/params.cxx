// *************** Test of parameters *********************

#include "params.hxx"

// Initialization of the physics
// ############################################################################
int Params::init(bool restarting) {
    TRACE("Halt in Params::init");

    // Get the option (before any sections) in the BOUT.inp file
    Options *options = Options::getRoot();

    // Load from the geometry
    // ************************************************************************
    Options *geom = options->getSection("geom");
    geom->get("radius", radius, 0.0);
    geom->get("len"   , len   , 0.0);
    // ************************************************************************

    // Load the input
    // ************************************************************************
    Options *input = options->getSection("input");
    input->get("n0"             , n0             , 0.0);
    input->get("Te0"            , Te0            , 0.0);
    input->get("Ti0"            , Ti0            , 0.0);
    input->get("B0"             , B0             , 0.0);
    input->get("S"              , S              , 0.0);
    input->get("nuENUnormalized", nuENUnormalized, 0.0);
    input->get("nuINUnormalized", nuINUnormalized, 0.0);
    // ************************************************************************

    Parameters params(radius, len, n0, Te0, Ti0, B0, S,
                      nuENUnormalized, nuINUnormalized);

    // Get the variables
    // ************************************************************************
    Lx     = params.getLx();
    Ly     = params.getLy();
    nuEI   = params.getNuEINorm();
    nuEN   = params.getNuENNorm();
    nuIN   = params.getNuINNorm();
    SNorm  = params.getSNorm();
    mu     = params.getMu();
    Lambda = params.getLambda();
    beta   = params.getBeta();
    omCI   = params.getOmCI();
    rhoS   = params.getRhoS();
    // ************************************************************************

    output << "\n\n\n\n\n\n\nSaving" << std::endl;
    SAVE_ONCE2(Lx, Ly);
    SAVE_ONCE4(nuEI, nuEN, nuIN, SNorm);
    SAVE_ONCE3(mu, Lambda, beta);
    SAVE_ONCE2(omCI, rhoS);

    // Finalize
    dump.write();
    dump.close();

    output << "\nFinished saving, now quitting\n\n\n\n\n\n" << std::endl;

    // Wait for all processors to write data
    MPI_Barrier(BoutComm::get());

    return 0;
}
// ############################################################################

// Solving the equations
// ############################################################################
int Params::rhs(BoutReal t) {
    return 0;
}
// ############################################################################

// Create a simple main() function
BOUTMAIN(Params);

// Destructor
Params::~Params(){
    TRACE("Params::~Params");
}
