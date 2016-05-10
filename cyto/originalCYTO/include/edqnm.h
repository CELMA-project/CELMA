#include <integrator.h>

void step_edqnm(HDF_DS *ds, PARA* p,long nx,long ny,double xmin,double ymin,double dx,double dy,double xmax,double*,double*,char*);

#ifndef linux
void ieee_trap(int ,int ,struct sigcontext  ,char *);
void mysignalhandler(int ,int ,struct sigcontext ,char *);
#endif
