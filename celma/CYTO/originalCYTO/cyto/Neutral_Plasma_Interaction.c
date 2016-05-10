#include <mpi.h>
#include <utilities.h>
#include <par_helper.h>


// Calculates the ionisation,recombination and charge exchange rates in GB normlised units and 
// per electron density, that is Gamma_ion/n for example.

void Neutral_plasma_interaction(HDF_DS *d,PARA *p,double*** G_ion, double*** G_cx, double*** G_rec, 
                             	   double ***Momentum_source, double***H_slow_0, double***exp_n_0, double***exp_t_0)
{
	int ip,ir, iz;
	int nx,ny,nz;
	
	double IonEnergy;
    double n_e, n_n;
    double tphi;
	double norm;
    

    nx = d->lnx; 
    ny = d->lny;  
    nz = d->lnz;


	// Ionisation potentials:
    // H = 13.5984 eV
    // Ar = 15.7596 eV according to periodic table of the elements

	if(p->Mi > 10.)
		IonEnergy = p->Te/15.7596; //Te in ev/ Ionisation for Argon
	else
		IonEnergy = p->Te/13.5984; // Te in EV / Hydrogen (and everything else, really

	
	// Ionisation rate Te/Phi_ion = T_ir
	// results all in meters and seconds. Rates in per cubic meter per second
	// Vi normalise to density p->n_0
	// and gyro bohm units -> multiplication of rates by unit-time (in u seconds).
   
	norm = p->n0*1.e19*p->unit_time*1e-6;
	
    // ORNL-6086 V1 and R.J. Goldston, P.H. Rutherford, Introduction to Plasma Physics, IOP, 1995
    // normalising density p->n0 in 10^19/m^3
	// assuming n_i = n_e

	FORALL_BD
    {
        tphi = exp_t_0[iz][ip][ir]*IonEnergy;
       
        n_e = exp_n_0[iz][ip][ir]*norm;
        n_n = exp(H_slow_0[iz][ip][ir])*norm;

		// Ionisation rate
      	// 2e-13/(6+T_ir ) * sqrt(T_ir) (e^-T_ir) n_e n_n 
		G_ion[iz][ip][ir] = 2.e-13/(6.+ tphi )*sqrt(tphi) *exp(-1./tphi)*n_n;
	

        // Recombination rate
        // 7e-20 sqrt(T_ir) n_e^2
		G_rec[iz][ip][ir] = .7e-20*n_e*sqrt(tphi);

        // Charge exchange rate 
		// we approximate for sub sound speeds of the neutrals (v_i >> v_n)
		// and small ion velocities
		// sqrt[ (v_i-v_n)² + 2/MT_i)
		// as C_s only
        // 3e-19 c_s n_n  n_e
		G_cx[iz][ip][ir] = 3.e-19*n_n*p->c_s;
	}
    
	
}
