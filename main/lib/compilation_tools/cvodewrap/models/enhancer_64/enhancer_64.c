#include "enhancer_64.h"
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_bandpre.h>
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_sparse.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <symbolic_functions.c>
#include <udata.h>
#include <math.h>
#include <mex.h>

#define pi 3.141592653589793


 int xdot_enhancer_64(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int ix;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  memset(xdot_tmp,0,sizeof(double)*14);
xdot_tmp[0] = (2.0*(pow(k[0],2))*p[2]*x_tmp[13] - 2.0*(pow(k[0],2))*p[3]*x_tmp[0] + k[0]*p[3]*x_tmp[1] + k[0]*p[2]*x_tmp[7])/(pow(k[0],2));
xdot_tmp[1] = -(1.0*(k[0]*p[3]*x_tmp[1] - 1.0*k[0]*p[2]*x_tmp[7]))/k[0];
xdot_tmp[2] = ((pow(k[0],2))*p[0]*x_tmp[4] + 2.0*(pow(k[0],2))*p[1]*x_tmp[9] + k[0]*p[1]*x_tmp[7] - 2.0*(pow(k[0],3))*p[0]*x_tmp[2]*x_tmp[8] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4]*x_tmp[6] + (pow(k[0],2))*p[0]*x_tmp[6]*x_tmp[8])/(pow(k[0],2));
xdot_tmp[3] = (144.0*p[4] + (pow(k[0],2))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[5]*x_tmp[3] + 2.0*(pow(k[0],2))*p[1]*x_tmp[11] + k[0]*p[1]*x_tmp[7] + k[0]*p[5]*x_tmp[8] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4]*x_tmp[8] + (pow(k[0],2))*p[0]*x_tmp[6]*x_tmp[8])/(pow(k[0],2));
xdot_tmp[4] = ((pow(k[0],2))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[5]*x_tmp[4] + (pow(k[0],2))*p[1]*x_tmp[9] + (pow(k[0],2))*p[1]*x_tmp[11] + k[0]*p[1]*x_tmp[7] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[2]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[4]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[4]*x_tmp[8] + (pow(k[0],2))*p[0]*x_tmp[6]*x_tmp[8])/(pow(k[0],2));
xdot_tmp[5] = ((pow(k[0],2))*p[0]*x_tmp[4] - 2.0*p[0]*((pow(k[0],3))*x_tmp[4]*x_tmp[7] + (pow(k[0],3))*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*x_tmp[8]*x_tmp[9] - 1.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) - 1.0*k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + 3.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]) - 2.0*(pow(k[0],2))*p[1]*x_tmp[5] + k[0]*p[1]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[6]*x_tmp[8] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[9])/(pow(k[0],2));
xdot_tmp[6] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[4] - 1.0*k[0]*p[1]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[6]*x_tmp[8]))/k[0];
xdot_tmp[7] = ((pow(k[0],2))*p[0]*x_tmp[4] - 1.0*k[0]*p[1]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[6]*x_tmp[8])/k[0];
xdot_tmp[8] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[4] - 12.0*p[4] - 1.0*k[0]*p[1]*x_tmp[7] + k[0]*p[5]*x_tmp[8] + (pow(k[0],2))*p[0]*x_tmp[6]*x_tmp[8]))/k[0];
xdot_tmp[9] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[4]*x_tmp[7] + (pow(k[0],3))*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*x_tmp[8]*x_tmp[9] - 1.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) - 1.0*k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + 3.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]) - 1.0*(pow(k[0],2))*p[1]*x_tmp[5] + (pow(k[0],2))*p[1]*x_tmp[9] + k[0]*p[1]*x_tmp[7] - 1.0*(pow(k[0],3))*p[0]*x_tmp[2]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[4]*x_tmp[6] + (pow(k[0],2))*p[0]*x_tmp[6]*x_tmp[8] + (pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[9]))/(pow(k[0],2));
xdot_tmp[10] = (p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[4] + (pow(k[0],3))*x_tmp[6]*x_tmp[12] + (pow(k[0],3))*x_tmp[8]*x_tmp[10] - 1.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) - 1.0*k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]) + (pow(k[0],2))*p[2]*x_tmp[9] - 1.0*(pow(k[0],2))*p[3]*x_tmp[10] + (pow(k[0],2))*p[1]*x_tmp[13] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[10])/(pow(k[0],2));
xdot_tmp[11] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[4]*x_tmp[7] + (pow(k[0],3))*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*x_tmp[8]*x_tmp[9] - 1.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) - 1.0*k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + 3.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]) - 1.0*(pow(k[0],2))*p[1]*x_tmp[5] + (pow(k[0],2))*p[1]*x_tmp[11] + (pow(k[0],2))*p[5]*x_tmp[11] + k[0]*p[1]*x_tmp[7] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[4]*x_tmp[8] + (pow(k[0],2))*p[0]*x_tmp[6]*x_tmp[8] + (pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[9]))/(pow(k[0],2));
xdot_tmp[12] = -(1.0*((pow(k[0],2))*p[3]*x_tmp[12] - 1.0*(pow(k[0],2))*p[2]*x_tmp[11] - 1.0*(pow(k[0],2))*p[1]*x_tmp[13] - 1.0*p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[4] + (pow(k[0],3))*x_tmp[6]*x_tmp[12] + (pow(k[0],3))*x_tmp[8]*x_tmp[10] - 1.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) - 1.0*k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]) + (pow(k[0],2))*p[5]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[10]))/(pow(k[0],2));
xdot_tmp[13] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[4] + (pow(k[0],3))*x_tmp[6]*x_tmp[12] + (pow(k[0],3))*x_tmp[8]*x_tmp[10] - 1.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) - 1.0*k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],2))*p[2]*x_tmp[5] + (pow(k[0],2))*p[1]*x_tmp[13] + (pow(k[0],2))*p[3]*x_tmp[13] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[10]))/(pow(k[0],2));

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_enhancer_64(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*28);
xBdot_tmp[0] = 2.0*p[3]*xB_tmp[0];
xBdot_tmp[1] = p[3]*xB_tmp[1] - (1.0*p[3]*xB_tmp[0])/k[0] - (1.0*p[0]*xB_tmp[10]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[12]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + (p[0]*xB_tmp[13]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
xBdot_tmp[2] = 2.0*k[0]*p[0]*x_tmp[8]*xB_tmp[2] + k[0]*p[0]*x_tmp[8]*xB_tmp[4] - 1.0*k[0]*p[0]*x_tmp[8]*xB_tmp[9];
xBdot_tmp[3] = ((2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*xB_tmp[3])/(pow(k[0],2)) + k[0]*p[0]*x_tmp[6]*xB_tmp[4] - 1.0*k[0]*p[0]*x_tmp[6]*xB_tmp[11];
xBdot_tmp[4] = (xB_tmp[4]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[5] - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*xB_tmp[2])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*xB_tmp[3])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*xB_tmp[9])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*xB_tmp[11])/(pow(k[0],2)) + k[0]*p[0]*xB_tmp[6] - 1.0*k[0]*p[0]*xB_tmp[7] + k[0]*p[0]*xB_tmp[8];
xBdot_tmp[5] = 2.0*p[1]*xB_tmp[5] - 1.0*p[1]*xB_tmp[9] - 1.0*p[1]*xB_tmp[11] - 1.0*p[2]*xB_tmp[13];
xBdot_tmp[6] = (xB_tmp[4]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*xB_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*xB_tmp[5]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*xB_tmp[3])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*xB_tmp[2])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[10])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[13])/(pow(k[0],2)) + k[0]*p[0]*x_tmp[8]*xB_tmp[6] - 1.0*k[0]*p[0]*x_tmp[8]*xB_tmp[7] + k[0]*p[0]*x_tmp[8]*xB_tmp[8];
xBdot_tmp[7] = p[1]*xB_tmp[7] - 1.0*p[1]*xB_tmp[6] - 1.0*p[2]*xB_tmp[1] - 1.0*p[1]*xB_tmp[8] - (1.0*p[2]*xB_tmp[0])/k[0] - (1.0*p[1]*xB_tmp[2])/k[0] - (1.0*p[1]*xB_tmp[3])/k[0] - (1.0*p[1]*xB_tmp[4])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*xB_tmp[9])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*xB_tmp[11])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*xB_tmp[5])/(pow(k[0],2));
xBdot_tmp[8] = (xB_tmp[4]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*xB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*xB_tmp[5]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*xB_tmp[2])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*xB_tmp[10])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*xB_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*xB_tmp[13])/(pow(k[0],2)) + ((k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*xB_tmp[8])/k[0] + k[0]*p[0]*x_tmp[6]*xB_tmp[6] - 1.0*k[0]*p[0]*x_tmp[6]*xB_tmp[7];
xBdot_tmp[9] = (xB_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[4] - 1.0*p[2]*xB_tmp[10] - 2.0*p[1]*xB_tmp[2] - 2.0*k[0]*p[0]*x_tmp[8]*xB_tmp[5] + k[0]*p[0]*x_tmp[8]*xB_tmp[11];
xBdot_tmp[10] = (xB_tmp[10]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[8]*xB_tmp[12] - 1.0*k[0]*p[0]*x_tmp[8]*xB_tmp[13];
xBdot_tmp[11] = (xB_tmp[11]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[4] - 1.0*p[2]*xB_tmp[12] - 2.0*p[1]*xB_tmp[3] - 2.0*k[0]*p[0]*x_tmp[6]*xB_tmp[5] + k[0]*p[0]*x_tmp[6]*xB_tmp[9];
xBdot_tmp[12] = (xB_tmp[12]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[6]*xB_tmp[10] - 1.0*k[0]*p[0]*x_tmp[6]*xB_tmp[13];
xBdot_tmp[13] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[13])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[10] - 1.0*p[1]*xB_tmp[12] - 2.0*p[2]*xB_tmp[0];
xBdot_tmp[14] = 2.0*p[3]*xB_tmp[14];
xBdot_tmp[15] = p[3]*xB_tmp[15] - (1.0*p[3]*xB_tmp[14])/k[0] - (1.0*p[0]*xB_tmp[24]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[26]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + (p[0]*xB_tmp[27]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
xBdot_tmp[16] = 2.0*k[0]*p[0]*x_tmp[8]*xB_tmp[16] + k[0]*p[0]*x_tmp[8]*xB_tmp[18] - 1.0*k[0]*p[0]*x_tmp[8]*xB_tmp[23];
xBdot_tmp[17] = ((2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*xB_tmp[17])/(pow(k[0],2)) + k[0]*p[0]*x_tmp[6]*xB_tmp[18] - 1.0*k[0]*p[0]*x_tmp[6]*xB_tmp[25];
xBdot_tmp[18] = (xB_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[19] - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*xB_tmp[16])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*xB_tmp[17])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*xB_tmp[23])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*xB_tmp[25])/(pow(k[0],2)) + k[0]*p[0]*xB_tmp[20] - 1.0*k[0]*p[0]*xB_tmp[21] + k[0]*p[0]*xB_tmp[22];
xBdot_tmp[19] = 2.0*p[1]*xB_tmp[19] - 1.0*p[1]*xB_tmp[23] - 1.0*p[1]*xB_tmp[25] - 1.0*p[2]*xB_tmp[27];
xBdot_tmp[20] = (xB_tmp[18]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*xB_tmp[25]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*xB_tmp[19]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*xB_tmp[17])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*xB_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[26])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[27])/(pow(k[0],2)) + k[0]*p[0]*x_tmp[8]*xB_tmp[20] - 1.0*k[0]*p[0]*x_tmp[8]*xB_tmp[21] + k[0]*p[0]*x_tmp[8]*xB_tmp[22];
xBdot_tmp[21] = p[1]*xB_tmp[21] - 1.0*p[1]*xB_tmp[20] - 1.0*p[2]*xB_tmp[15] - 1.0*p[1]*xB_tmp[22] - (1.0*p[2]*xB_tmp[14])/k[0] - (1.0*p[1]*xB_tmp[16])/k[0] - (1.0*p[1]*xB_tmp[17])/k[0] - (1.0*p[1]*xB_tmp[18])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*xB_tmp[23])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*xB_tmp[25])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*xB_tmp[19])/(pow(k[0],2));
xBdot_tmp[22] = (xB_tmp[18]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[25]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[17]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*xB_tmp[19]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*xB_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*xB_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*xB_tmp[26])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*xB_tmp[27])/(pow(k[0],2)) + ((k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*xB_tmp[22])/k[0] + k[0]*p[0]*x_tmp[6]*xB_tmp[20] - 1.0*k[0]*p[0]*x_tmp[6]*xB_tmp[21];
xBdot_tmp[23] = (xB_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[18] - 1.0*p[2]*xB_tmp[24] - 2.0*p[1]*xB_tmp[16] - 2.0*k[0]*p[0]*x_tmp[8]*xB_tmp[19] + k[0]*p[0]*x_tmp[8]*xB_tmp[25];
xBdot_tmp[24] = (xB_tmp[24]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[8]*xB_tmp[26] - 1.0*k[0]*p[0]*x_tmp[8]*xB_tmp[27];
xBdot_tmp[25] = (xB_tmp[25]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[18] - 1.0*p[2]*xB_tmp[26] - 2.0*p[1]*xB_tmp[17] - 2.0*k[0]*p[0]*x_tmp[6]*xB_tmp[19] + k[0]*p[0]*x_tmp[6]*xB_tmp[23];
xBdot_tmp[26] = (xB_tmp[26]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[6]*xB_tmp[24] - 1.0*k[0]*p[0]*x_tmp[6]*xB_tmp[27];
xBdot_tmp[27] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[27])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[24] - 1.0*p[1]*xB_tmp[26] - 2.0*p[2]*xB_tmp[14];

  for (ixB=0; ixB<28; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_enhancer_64(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
{
  int iyp;
  int ip;
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  int np = *data->np;
  int ny = *data->ny;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *qBdot_tmp = N_VGetArrayPointer(qBdot);
  memset(qBdot_tmp,0,sizeof(double)*2*np);
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
qBdot_tmp[0+ip*ny] = (xB_tmp[9]*((pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],3))*x_tmp[2]*x_tmp[8] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[6] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[7] + (pow(k[0],2))*x_tmp[6]*x_tmp[8] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[11]*((pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],3))*x_tmp[3]*x_tmp[6] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[7] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[8] + (pow(k[0],2))*x_tmp[6]*x_tmp[8] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[2]*((pow(k[0],2))*x_tmp[4] - 2.0*(pow(k[0],3))*x_tmp[2]*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[6] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*((pow(k[0],2))*x_tmp[4] - 2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[6] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[8] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[6]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/k[0] - (1.0*xB_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/k[0] + (xB_tmp[8]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/k[0] + (xB_tmp[10]*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[4] + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[12]*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[4] + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[13]*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[4] + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[4]*((pow(k[0],3))*x_tmp[3]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],3))*x_tmp[2]*x_tmp[8] + (pow(k[0],3))*x_tmp[4]*x_tmp[6] + (pow(k[0],3))*x_tmp[4]*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[5]*((pow(k[0],2))*x_tmp[4] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[7] + (pow(k[0],2))*x_tmp[6]*x_tmp[8] + 2.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + 2.0*k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + 2.0*k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) - 6.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = (xB_tmp[23]*((pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],3))*x_tmp[2]*x_tmp[8] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[6] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[7] + (pow(k[0],2))*x_tmp[6]*x_tmp[8] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[25]*((pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],3))*x_tmp[3]*x_tmp[6] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[7] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[8] + (pow(k[0],2))*x_tmp[6]*x_tmp[8] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[16]*((pow(k[0],2))*x_tmp[4] - 2.0*(pow(k[0],3))*x_tmp[2]*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[6] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[17]*((pow(k[0],2))*x_tmp[4] - 2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[6] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[8] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[20]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/k[0] - (1.0*xB_tmp[21]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/k[0] + (xB_tmp[22]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/k[0] + (xB_tmp[24]*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[4] + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[26]*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[4] + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[27]*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[4] + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[18]*((pow(k[0],3))*x_tmp[3]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],3))*x_tmp[2]*x_tmp[8] + (pow(k[0],3))*x_tmp[4]*x_tmp[6] + (pow(k[0],3))*x_tmp[4]*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[19]*((pow(k[0],2))*x_tmp[4] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[7] + (pow(k[0],2))*x_tmp[6]*x_tmp[8] + 2.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + 2.0*k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + 2.0*k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) - 6.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]))/(pow(k[0],2));

  } break;

  case 1: {
qBdot_tmp[0+ip*ny] = x_tmp[7]*xB_tmp[7] - 1.0*x_tmp[7]*xB_tmp[6] - 1.0*x_tmp[7]*xB_tmp[8] - 1.0*x_tmp[13]*xB_tmp[10] - 1.0*x_tmp[13]*xB_tmp[12] + x_tmp[13]*xB_tmp[13] - (1.0*(k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[5])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[7] + 2.0*(pow(k[0],2))*x_tmp[9])*xB_tmp[2])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[7] + 2.0*(pow(k[0],2))*x_tmp[11])*xB_tmp[3])/(pow(k[0],2)) - (1.0*xB_tmp[4]*(k[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) + (xB_tmp[9]*(k[0]*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) + (xB_tmp[11]*(k[0]*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[11]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = x_tmp[7]*xB_tmp[21] - 1.0*x_tmp[7]*xB_tmp[20] - 1.0*x_tmp[7]*xB_tmp[22] - 1.0*x_tmp[13]*xB_tmp[24] - 1.0*x_tmp[13]*xB_tmp[26] + x_tmp[13]*xB_tmp[27] - (1.0*(k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[19])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[7] + 2.0*(pow(k[0],2))*x_tmp[9])*xB_tmp[16])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[7] + 2.0*(pow(k[0],2))*x_tmp[11])*xB_tmp[17])/(pow(k[0],2)) - (1.0*xB_tmp[18]*(k[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) + (xB_tmp[23]*(k[0]*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) + (xB_tmp[25]*(k[0]*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[11]))/(pow(k[0],2));

  } break;

  case 2: {
qBdot_tmp[0+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[1] - 1.0*x_tmp[5]*xB_tmp[13] - 1.0*x_tmp[9]*xB_tmp[10] - 1.0*x_tmp[11]*xB_tmp[12] - (1.0*(k[0]*x_tmp[7] + 2.0*(pow(k[0],2))*x_tmp[13])*xB_tmp[0])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = - 1.0*x_tmp[7]*xB_tmp[15] - 1.0*x_tmp[5]*xB_tmp[27] - 1.0*x_tmp[9]*xB_tmp[24] - 1.0*x_tmp[11]*xB_tmp[26] - (1.0*(k[0]*x_tmp[7] + 2.0*(pow(k[0],2))*x_tmp[13])*xB_tmp[14])/(pow(k[0],2));

  } break;

  case 3: {
qBdot_tmp[0+ip*ny] = x_tmp[1]*xB_tmp[1] + x_tmp[10]*xB_tmp[10] + x_tmp[12]*xB_tmp[12] + x_tmp[13]*xB_tmp[13] - (1.0*(k[0]*x_tmp[1] - 2.0*(pow(k[0],2))*x_tmp[0])*xB_tmp[0])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = x_tmp[1]*xB_tmp[15] + x_tmp[10]*xB_tmp[24] + x_tmp[12]*xB_tmp[26] + x_tmp[13]*xB_tmp[27] - (1.0*(k[0]*x_tmp[1] - 2.0*(pow(k[0],2))*x_tmp[0])*xB_tmp[14])/(pow(k[0],2));

  } break;

  case 4: {
qBdot_tmp[0+ip*ny] = - (144.0*xB_tmp[3])/(pow(k[0],2)) - (12.0*xB_tmp[8])/k[0];
qBdot_tmp[1+ip*ny] = - (144.0*xB_tmp[17])/(pow(k[0],2)) - (12.0*xB_tmp[22])/k[0];

  } break;

  case 5: {
qBdot_tmp[0+ip*ny] = x_tmp[4]*xB_tmp[4] + x_tmp[8]*xB_tmp[8] + x_tmp[11]*xB_tmp[11] + x_tmp[12]*xB_tmp[12] - (1.0*(k[0]*x_tmp[8] - 2.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[3])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = x_tmp[4]*xB_tmp[18] + x_tmp[8]*xB_tmp[22] + x_tmp[11]*xB_tmp[25] + x_tmp[12]*xB_tmp[26] - (1.0*(k[0]*x_tmp[8] - 2.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[17])/(pow(k[0],2));

  } break;

  }
  }

  for (iyp=0; iyp<2*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_enhancer_64(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*14);
x0_tmp[0] = (k[14]*k[28])/(pow(k[0],2));
x0_tmp[1] = (k[4]*k[18])/k[0];
x0_tmp[2] = (k[5]*k[19])/(pow(k[0],2));
x0_tmp[3] = (k[12]*k[26])/(pow(k[0],2));
x0_tmp[4] = (k[7]*k[21])/(pow(k[0],2));
x0_tmp[5] = (k[9]*k[23])/(pow(k[0],2));
x0_tmp[6] = (k[1]*k[15] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[7] = (k[2]*k[16])/k[0];
x0_tmp[8] = (k[3]*k[17])/k[0];
x0_tmp[9] = (k[6]*k[20])/(pow(k[0],2));
x0_tmp[10] = (k[8]*k[22])/(pow(k[0],2));
x0_tmp[11] = (k[10]*k[24])/(pow(k[0],2));
x0_tmp[12] = (k[13]*k[27])/(pow(k[0],2));
x0_tmp[13] = (k[11]*k[25])/(pow(k[0],2));
  
  
  return;
}


 int Jv_enhancer_64(N_Vector v, N_Vector Jv, realtype t,
  	N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp)
{
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *v_tmp = N_VGetArrayPointer(v);
  double *Jv_tmp = N_VGetArrayPointer(Jv);
  memset(Jv_tmp,0,sizeof(double)*14);
Jv_tmp[0] = 2.0*p[2]*v_tmp[13] - 2.0*p[3]*v_tmp[0] + (p[3]*v_tmp[1])/k[0] + (p[2]*v_tmp[7])/k[0];
Jv_tmp[1] = p[2]*v_tmp[7] - 1.0*p[3]*v_tmp[1];
Jv_tmp[2] = 2.0*p[1]*v_tmp[9] + (p[1]*v_tmp[7])/k[0] + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*v_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*v_tmp[8])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*v_tmp[6])/(pow(k[0],2)) - 2.0*k[0]*p[0]*v_tmp[2]*x_tmp[8];
Jv_tmp[3] = 2.0*p[1]*v_tmp[11] + (p[1]*v_tmp[7])/k[0] + (v_tmp[8]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*v_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*v_tmp[6])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*v_tmp[3])/(pow(k[0],2));
Jv_tmp[4] = p[1]*v_tmp[9] + p[1]*v_tmp[11] - (1.0*v_tmp[4]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (p[1]*v_tmp[7])/k[0] - (1.0*v_tmp[8]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*v_tmp[6]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*v_tmp[2]*x_tmp[8];
Jv_tmp[5] = p[0]*v_tmp[4] - 2.0*p[1]*v_tmp[5] + (v_tmp[8]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (v_tmp[6]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*v_tmp[7])/(pow(k[0],2)) + 2.0*k[0]*p[0]*v_tmp[9]*x_tmp[8] + 2.0*k[0]*p[0]*v_tmp[11]*x_tmp[6];
Jv_tmp[6] = p[1]*v_tmp[7] - 1.0*k[0]*p[0]*v_tmp[4] - 1.0*k[0]*p[0]*v_tmp[6]*x_tmp[8] - 1.0*k[0]*p[0]*v_tmp[8]*x_tmp[6];
Jv_tmp[7] = k[0]*p[0]*v_tmp[4] - 1.0*p[1]*v_tmp[7] + k[0]*p[0]*v_tmp[6]*x_tmp[8] + k[0]*p[0]*v_tmp[8]*x_tmp[6];
Jv_tmp[8] = p[1]*v_tmp[7] - 1.0*k[0]*p[0]*v_tmp[4] - (1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*v_tmp[8])/k[0] - 1.0*k[0]*p[0]*v_tmp[6]*x_tmp[8];
Jv_tmp[9] = p[1]*v_tmp[5] - (1.0*v_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (v_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (v_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*v_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*v_tmp[4])/(pow(k[0],2)) + k[0]*p[0]*v_tmp[2]*x_tmp[8] - 1.0*k[0]*p[0]*v_tmp[11]*x_tmp[6];
Jv_tmp[10] = p[2]*v_tmp[9] + p[1]*v_tmp[13] - (1.0*v_tmp[10]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*v_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*v_tmp[6])/(pow(k[0],2)) + (p[0]*v_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[12]*x_tmp[6];
Jv_tmp[11] = p[1]*v_tmp[5] + (v_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (v_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*v_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*v_tmp[4])/(pow(k[0],2)) - (1.0*v_tmp[11]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*v_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*v_tmp[9]*x_tmp[8];
Jv_tmp[12] = p[2]*v_tmp[11] + p[1]*v_tmp[13] - (1.0*v_tmp[12]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*v_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*v_tmp[6])/(pow(k[0],2)) + (p[0]*v_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[10]*x_tmp[8];
Jv_tmp[13] = p[2]*v_tmp[5] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*v_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*v_tmp[8])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*v_tmp[6])/(pow(k[0],2)) - (1.0*p[0]*v_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*v_tmp[10]*x_tmp[8] + k[0]*p[0]*v_tmp[12]*x_tmp[6];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_enhancer_64(N_Vector vB, N_Vector JvB, realtype t,
  	N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data, N_Vector tmpB)
{
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *vB_tmp = N_VGetArrayPointer(vB);
  double *JvB_tmp = N_VGetArrayPointer(JvB);
  memset(JvB_tmp,0,sizeof(double)*14);
JvB_tmp[0] = 2.0*p[3]*vB_tmp[0];
JvB_tmp[1] = p[3]*vB_tmp[1] - (1.0*p[3]*vB_tmp[0])/k[0] - (1.0*p[0]*vB_tmp[10]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*p[0]*vB_tmp[12]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + (p[0]*vB_tmp[13]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
JvB_tmp[2] = 2.0*k[0]*p[0]*vB_tmp[2]*x_tmp[8] + k[0]*p[0]*vB_tmp[4]*x_tmp[8] - 1.0*k[0]*p[0]*vB_tmp[9]*x_tmp[8];
JvB_tmp[3] = ((2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*vB_tmp[3])/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[4]*x_tmp[6] - 1.0*k[0]*p[0]*vB_tmp[11]*x_tmp[6];
JvB_tmp[4] = (vB_tmp[4]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[5] - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*vB_tmp[2])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*vB_tmp[3])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*vB_tmp[9])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*vB_tmp[11])/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[6] - 1.0*k[0]*p[0]*vB_tmp[7] + k[0]*p[0]*vB_tmp[8];
JvB_tmp[5] = 2.0*p[1]*vB_tmp[5] - 1.0*p[1]*vB_tmp[9] - 1.0*p[1]*vB_tmp[11] - 1.0*p[2]*vB_tmp[13];
JvB_tmp[6] = (vB_tmp[4]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - (1.0*vB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*vB_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*vB_tmp[5]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*vB_tmp[3])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*vB_tmp[2])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*vB_tmp[10])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*vB_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*vB_tmp[13])/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[6]*x_tmp[8] - 1.0*k[0]*p[0]*vB_tmp[7]*x_tmp[8] + k[0]*p[0]*vB_tmp[8]*x_tmp[8];
JvB_tmp[7] = p[1]*vB_tmp[7] - 1.0*p[1]*vB_tmp[6] - 1.0*p[2]*vB_tmp[1] - 1.0*p[1]*vB_tmp[8] - (1.0*p[2]*vB_tmp[0])/k[0] - (1.0*p[1]*vB_tmp[2])/k[0] - (1.0*p[1]*vB_tmp[3])/k[0] - (1.0*p[1]*vB_tmp[4])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*vB_tmp[9])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*vB_tmp[11])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*vB_tmp[5])/(pow(k[0],2));
JvB_tmp[8] = (vB_tmp[4]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*vB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*vB_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*vB_tmp[3]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*vB_tmp[5]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*vB_tmp[2])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*vB_tmp[10])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*vB_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*vB_tmp[13])/(pow(k[0],2)) + ((k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*vB_tmp[8])/k[0] + k[0]*p[0]*vB_tmp[6]*x_tmp[6] - 1.0*k[0]*p[0]*vB_tmp[7]*x_tmp[6];
JvB_tmp[9] = (vB_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[4] - 1.0*p[2]*vB_tmp[10] - 2.0*p[1]*vB_tmp[2] - 2.0*k[0]*p[0]*vB_tmp[5]*x_tmp[8] + k[0]*p[0]*vB_tmp[11]*x_tmp[8];
JvB_tmp[10] = (vB_tmp[10]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[12]*x_tmp[8] - 1.0*k[0]*p[0]*vB_tmp[13]*x_tmp[8];
JvB_tmp[11] = (vB_tmp[11]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[4] - 1.0*p[2]*vB_tmp[12] - 2.0*p[1]*vB_tmp[3] - 2.0*k[0]*p[0]*vB_tmp[5]*x_tmp[6] + k[0]*p[0]*vB_tmp[9]*x_tmp[6];
JvB_tmp[12] = (vB_tmp[12]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[10]*x_tmp[6] - 1.0*k[0]*p[0]*vB_tmp[13]*x_tmp[6];
JvB_tmp[13] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*vB_tmp[13])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[10] - 1.0*p[1]*vB_tmp[12] - 2.0*p[2]*vB_tmp[0];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_enhancer_64(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_enhancer_64(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_enhancer_64(long int N, realtype t, N_Vector x,
  	N_Vector xdot, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  memset(J->data,0,sizeof(double)*196);
J->data[0] = -2.0*p[3];
J->data[14] = p[3]/k[0];
J->data[15] = -1.0*p[3];
J->data[24] = (p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
J->data[26] = (p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
J->data[27] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
J->data[30] = -2.0*k[0]*p[0]*x_tmp[8];
J->data[32] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[37] = k[0]*p[0]*x_tmp[8];
J->data[45] = -(1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[46] = -1.0*k[0]*p[0]*x_tmp[6];
J->data[53] = k[0]*p[0]*x_tmp[6];
J->data[58] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
J->data[59] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])/(pow(k[0],2));
J->data[60] = -(1.0*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[61] = p[0];
J->data[62] = -1.0*k[0]*p[0];
J->data[63] = k[0]*p[0];
J->data[64] = -1.0*k[0]*p[0];
J->data[65] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[67] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[75] = -2.0*p[1];
J->data[79] = p[1];
J->data[81] = p[1];
J->data[83] = p[2];
J->data[86] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[87] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[88] = -(1.0*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[89] = ((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[90] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[91] = k[0]*p[0]*x_tmp[8];
J->data[92] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[93] = (p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[94] = (p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[95] = (p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[96] = (p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[97] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]))/(pow(k[0],2));
J->data[98] = p[2]/k[0];
J->data[99] = p[2];
J->data[100] = p[1]/k[0];
J->data[101] = p[1]/k[0];
J->data[102] = p[1]/k[0];
J->data[103] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
J->data[104] = p[1];
J->data[105] = -1.0*p[1];
J->data[106] = p[1];
J->data[107] = (p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[109] = (p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[114] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[115] = (k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6])/(pow(k[0],2));
J->data[116] = -(1.0*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[117] = ((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
J->data[118] = -1.0*k[0]*p[0]*x_tmp[6];
J->data[119] = k[0]*p[0]*x_tmp[6];
J->data[120] = -(1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6]))/k[0];
J->data[121] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
J->data[122] = (p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])/(pow(k[0],2));
J->data[123] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
J->data[124] = (p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])/(pow(k[0],2));
J->data[125] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10]))/(pow(k[0],2));
J->data[128] = 2.0*p[1];
J->data[130] = p[1];
J->data[131] = 2.0*k[0]*p[0]*x_tmp[8];
J->data[135] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[136] = p[2];
J->data[137] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[150] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[152] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[153] = k[0]*p[0]*x_tmp[8];
J->data[157] = 2.0*p[1];
J->data[158] = p[1];
J->data[159] = 2.0*k[0]*p[0]*x_tmp[6];
J->data[163] = -1.0*k[0]*p[0]*x_tmp[6];
J->data[165] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[166] = p[2];
J->data[178] = -1.0*k[0]*p[0]*x_tmp[6];
J->data[180] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[181] = k[0]*p[0]*x_tmp[6];
J->data[182] = 2.0*p[2];
J->data[192] = p[1];
J->data[194] = p[1];
J->data[195] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));

  for (iJ=0; iJ<196; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_enhancer_64(realtype t, N_Vector x,
  	N_Vector xdot, SlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  SlsSetToZero(J);
  J->rowvals[0] = 0;
  J->rowvals[1] = 0;
  J->rowvals[2] = 1;
  J->rowvals[3] = 10;
  J->rowvals[4] = 12;
  J->rowvals[5] = 13;
  J->rowvals[6] = 2;
  J->rowvals[7] = 4;
  J->rowvals[8] = 9;
  J->rowvals[9] = 3;
  J->rowvals[10] = 4;
  J->rowvals[11] = 11;
  J->rowvals[12] = 2;
  J->rowvals[13] = 3;
  J->rowvals[14] = 4;
  J->rowvals[15] = 5;
  J->rowvals[16] = 6;
  J->rowvals[17] = 7;
  J->rowvals[18] = 8;
  J->rowvals[19] = 9;
  J->rowvals[20] = 11;
  J->rowvals[21] = 5;
  J->rowvals[22] = 9;
  J->rowvals[23] = 11;
  J->rowvals[24] = 13;
  J->rowvals[25] = 2;
  J->rowvals[26] = 3;
  J->rowvals[27] = 4;
  J->rowvals[28] = 5;
  J->rowvals[29] = 6;
  J->rowvals[30] = 7;
  J->rowvals[31] = 8;
  J->rowvals[32] = 9;
  J->rowvals[33] = 10;
  J->rowvals[34] = 11;
  J->rowvals[35] = 12;
  J->rowvals[36] = 13;
  J->rowvals[37] = 0;
  J->rowvals[38] = 1;
  J->rowvals[39] = 2;
  J->rowvals[40] = 3;
  J->rowvals[41] = 4;
  J->rowvals[42] = 5;
  J->rowvals[43] = 6;
  J->rowvals[44] = 7;
  J->rowvals[45] = 8;
  J->rowvals[46] = 9;
  J->rowvals[47] = 11;
  J->rowvals[48] = 2;
  J->rowvals[49] = 3;
  J->rowvals[50] = 4;
  J->rowvals[51] = 5;
  J->rowvals[52] = 6;
  J->rowvals[53] = 7;
  J->rowvals[54] = 8;
  J->rowvals[55] = 9;
  J->rowvals[56] = 10;
  J->rowvals[57] = 11;
  J->rowvals[58] = 12;
  J->rowvals[59] = 13;
  J->rowvals[60] = 2;
  J->rowvals[61] = 4;
  J->rowvals[62] = 5;
  J->rowvals[63] = 9;
  J->rowvals[64] = 10;
  J->rowvals[65] = 11;
  J->rowvals[66] = 10;
  J->rowvals[67] = 12;
  J->rowvals[68] = 13;
  J->rowvals[69] = 3;
  J->rowvals[70] = 4;
  J->rowvals[71] = 5;
  J->rowvals[72] = 9;
  J->rowvals[73] = 11;
  J->rowvals[74] = 12;
  J->rowvals[75] = 10;
  J->rowvals[76] = 12;
  J->rowvals[77] = 13;
  J->rowvals[78] = 0;
  J->rowvals[79] = 10;
  J->rowvals[80] = 12;
  J->rowvals[81] = 13;
  J->colptrs[0] = 0;
  J->colptrs[1] = 1;
  J->colptrs[2] = 6;
  J->colptrs[3] = 9;
  J->colptrs[4] = 12;
  J->colptrs[5] = 21;
  J->colptrs[6] = 25;
  J->colptrs[7] = 37;
  J->colptrs[8] = 48;
  J->colptrs[9] = 60;
  J->colptrs[10] = 66;
  J->colptrs[11] = 69;
  J->colptrs[12] = 75;
  J->colptrs[13] = 78;
  J->colptrs[14] = 82;
J->data[0] = -2.0*p[3];
J->data[1] = p[3]/k[0];
J->data[2] = -1.0*p[3];
J->data[3] = (p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
J->data[4] = (p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
J->data[5] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
J->data[6] = -2.0*k[0]*p[0]*x_tmp[8];
J->data[7] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[8] = k[0]*p[0]*x_tmp[8];
J->data[9] = -(1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[10] = -1.0*k[0]*p[0]*x_tmp[6];
J->data[11] = k[0]*p[0]*x_tmp[6];
J->data[12] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
J->data[13] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])/(pow(k[0],2));
J->data[14] = -(1.0*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[15] = p[0];
J->data[16] = -1.0*k[0]*p[0];
J->data[17] = k[0]*p[0];
J->data[18] = -1.0*k[0]*p[0];
J->data[19] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[20] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[21] = -2.0*p[1];
J->data[22] = p[1];
J->data[23] = p[1];
J->data[24] = p[2];
J->data[25] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[26] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[27] = -(1.0*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[28] = ((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[29] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[30] = k[0]*p[0]*x_tmp[8];
J->data[31] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[32] = (p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[33] = (p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[34] = (p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[35] = (p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[36] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]))/(pow(k[0],2));
J->data[37] = p[2]/k[0];
J->data[38] = p[2];
J->data[39] = p[1]/k[0];
J->data[40] = p[1]/k[0];
J->data[41] = p[1]/k[0];
J->data[42] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
J->data[43] = p[1];
J->data[44] = -1.0*p[1];
J->data[45] = p[1];
J->data[46] = (p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[47] = (p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[48] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[49] = (k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6])/(pow(k[0],2));
J->data[50] = -(1.0*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[51] = ((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
J->data[52] = -1.0*k[0]*p[0]*x_tmp[6];
J->data[53] = k[0]*p[0]*x_tmp[6];
J->data[54] = -(1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6]))/k[0];
J->data[55] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
J->data[56] = (p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])/(pow(k[0],2));
J->data[57] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
J->data[58] = (p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])/(pow(k[0],2));
J->data[59] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10]))/(pow(k[0],2));
J->data[60] = 2.0*p[1];
J->data[61] = p[1];
J->data[62] = 2.0*k[0]*p[0]*x_tmp[8];
J->data[63] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[64] = p[2];
J->data[65] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[66] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2));
J->data[67] = -1.0*k[0]*p[0]*x_tmp[8];
J->data[68] = k[0]*p[0]*x_tmp[8];
J->data[69] = 2.0*p[1];
J->data[70] = p[1];
J->data[71] = 2.0*k[0]*p[0]*x_tmp[6];
J->data[72] = -1.0*k[0]*p[0]*x_tmp[6];
J->data[73] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[74] = p[2];
J->data[75] = -1.0*k[0]*p[0]*x_tmp[6];
J->data[76] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[77] = k[0]*p[0]*x_tmp[6];
J->data[78] = 2.0*p[2];
J->data[79] = p[1];
J->data[80] = p[1];
J->data[81] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
  return(0);
}


 int JBBand_enhancer_64(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_enhancer_64(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_enhancer_64(long int N, realtype t, N_Vector x,
  	N_Vector xB, N_Vector xdotB, DlsMat JB, void *user_data, 
  	N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
JB->data[0] = 2.0*p[3];
JB->data[1] = -(1.0*p[3])/k[0];
JB->data[7] = -(1.0*p[2])/k[0];
JB->data[13] = -2.0*p[2];
JB->data[15] = p[3];
JB->data[21] = -1.0*p[2];
JB->data[30] = 2.0*k[0]*p[0]*x_tmp[8];
JB->data[32] = -(1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
JB->data[34] = (2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])/(pow(k[0],2));
JB->data[35] = -(1.0*p[1])/k[0];
JB->data[36] = (2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])/(pow(k[0],2));
JB->data[37] = -2.0*p[1];
JB->data[45] = (2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
JB->data[46] = -(1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2));
JB->data[48] = (2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])/(pow(k[0],2));
JB->data[49] = -(1.0*p[1])/k[0];
JB->data[50] = -(1.0*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2));
JB->data[53] = -2.0*p[1];
JB->data[58] = k[0]*p[0]*x_tmp[8];
JB->data[59] = k[0]*p[0]*x_tmp[6];
JB->data[60] = ((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8])/(pow(k[0],2));
JB->data[62] = ((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])/(pow(k[0],2));
JB->data[63] = -(1.0*p[1])/k[0];
JB->data[64] = ((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])/(pow(k[0],2));
JB->data[65] = -1.0*p[1];
JB->data[67] = -1.0*p[1];
JB->data[74] = -1.0*p[0];
JB->data[75] = 2.0*p[1];
JB->data[76] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
JB->data[77] = (2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])/(pow(k[0],2));
JB->data[78] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2));
JB->data[79] = -2.0*k[0]*p[0]*x_tmp[8];
JB->data[81] = -2.0*k[0]*p[0]*x_tmp[6];
JB->data[88] = k[0]*p[0];
JB->data[90] = k[0]*p[0]*x_tmp[8];
JB->data[91] = -1.0*p[1];
JB->data[92] = k[0]*p[0]*x_tmp[6];
JB->data[102] = -1.0*k[0]*p[0];
JB->data[104] = -1.0*k[0]*p[0]*x_tmp[8];
JB->data[105] = p[1];
JB->data[106] = -1.0*k[0]*p[0]*x_tmp[6];
JB->data[116] = k[0]*p[0];
JB->data[118] = k[0]*p[0]*x_tmp[8];
JB->data[119] = -1.0*p[1];
JB->data[120] = (k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])/k[0];
JB->data[128] = -1.0*k[0]*p[0]*x_tmp[8];
JB->data[130] = ((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
JB->data[131] = -1.0*p[1];
JB->data[132] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
JB->data[133] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
JB->data[134] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2));
JB->data[135] = ((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8])/(pow(k[0],2));
JB->data[137] = k[0]*p[0]*x_tmp[6];
JB->data[141] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
JB->data[146] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]))/(pow(k[0],2));
JB->data[148] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10]))/(pow(k[0],2));
JB->data[149] = -1.0*p[2];
JB->data[150] = ((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8])/(pow(k[0],2));
JB->data[152] = k[0]*p[0]*x_tmp[6];
JB->data[153] = -1.0*p[1];
JB->data[157] = -1.0*k[0]*p[0]*x_tmp[6];
JB->data[158] = ((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])/(pow(k[0],2));
JB->data[159] = -1.0*p[1];
JB->data[160] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
JB->data[161] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
JB->data[162] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2));
JB->data[163] = k[0]*p[0]*x_tmp[8];
JB->data[165] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
JB->data[169] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
JB->data[174] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]))/(pow(k[0],2));
JB->data[176] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10]))/(pow(k[0],2));
JB->data[178] = k[0]*p[0]*x_tmp[8];
JB->data[179] = -1.0*p[2];
JB->data[180] = ((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
JB->data[181] = -1.0*p[1];
JB->data[183] = (p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2));
JB->data[187] = -1.0*p[2];
JB->data[188] = (p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])/(pow(k[0],2));
JB->data[190] = (p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])/(pow(k[0],2));
JB->data[192] = -1.0*k[0]*p[0]*x_tmp[8];
JB->data[194] = -1.0*k[0]*p[0]*x_tmp[6];
JB->data[195] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])/(pow(k[0],2));

  for (iJ=0; iJ<196; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_enhancer_64(realtype t, N_Vector x,
  	N_Vector xB, N_Vector xdotB, SlsMat JB, void *user_data, 
  	N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  SlsSetToZero(JB);
  JB->rowvals[0] = 0;
  JB->rowvals[1] = 1;
  JB->rowvals[2] = 7;
  JB->rowvals[3] = 13;
  JB->rowvals[4] = 1;
  JB->rowvals[5] = 7;
  JB->rowvals[6] = 2;
  JB->rowvals[7] = 4;
  JB->rowvals[8] = 6;
  JB->rowvals[9] = 7;
  JB->rowvals[10] = 8;
  JB->rowvals[11] = 9;
  JB->rowvals[12] = 3;
  JB->rowvals[13] = 4;
  JB->rowvals[14] = 6;
  JB->rowvals[15] = 7;
  JB->rowvals[16] = 8;
  JB->rowvals[17] = 11;
  JB->rowvals[18] = 2;
  JB->rowvals[19] = 3;
  JB->rowvals[20] = 4;
  JB->rowvals[21] = 6;
  JB->rowvals[22] = 7;
  JB->rowvals[23] = 8;
  JB->rowvals[24] = 9;
  JB->rowvals[25] = 11;
  JB->rowvals[26] = 4;
  JB->rowvals[27] = 5;
  JB->rowvals[28] = 6;
  JB->rowvals[29] = 7;
  JB->rowvals[30] = 8;
  JB->rowvals[31] = 9;
  JB->rowvals[32] = 11;
  JB->rowvals[33] = 4;
  JB->rowvals[34] = 6;
  JB->rowvals[35] = 7;
  JB->rowvals[36] = 8;
  JB->rowvals[37] = 4;
  JB->rowvals[38] = 6;
  JB->rowvals[39] = 7;
  JB->rowvals[40] = 8;
  JB->rowvals[41] = 4;
  JB->rowvals[42] = 6;
  JB->rowvals[43] = 7;
  JB->rowvals[44] = 8;
  JB->rowvals[45] = 2;
  JB->rowvals[46] = 4;
  JB->rowvals[47] = 5;
  JB->rowvals[48] = 6;
  JB->rowvals[49] = 7;
  JB->rowvals[50] = 8;
  JB->rowvals[51] = 9;
  JB->rowvals[52] = 11;
  JB->rowvals[53] = 1;
  JB->rowvals[54] = 6;
  JB->rowvals[55] = 8;
  JB->rowvals[56] = 9;
  JB->rowvals[57] = 10;
  JB->rowvals[58] = 12;
  JB->rowvals[59] = 13;
  JB->rowvals[60] = 3;
  JB->rowvals[61] = 4;
  JB->rowvals[62] = 5;
  JB->rowvals[63] = 6;
  JB->rowvals[64] = 7;
  JB->rowvals[65] = 8;
  JB->rowvals[66] = 9;
  JB->rowvals[67] = 11;
  JB->rowvals[68] = 1;
  JB->rowvals[69] = 6;
  JB->rowvals[70] = 8;
  JB->rowvals[71] = 10;
  JB->rowvals[72] = 11;
  JB->rowvals[73] = 12;
  JB->rowvals[74] = 13;
  JB->rowvals[75] = 1;
  JB->rowvals[76] = 5;
  JB->rowvals[77] = 6;
  JB->rowvals[78] = 8;
  JB->rowvals[79] = 10;
  JB->rowvals[80] = 12;
  JB->rowvals[81] = 13;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 4;
  JB->colptrs[2] = 6;
  JB->colptrs[3] = 12;
  JB->colptrs[4] = 18;
  JB->colptrs[5] = 26;
  JB->colptrs[6] = 33;
  JB->colptrs[7] = 37;
  JB->colptrs[8] = 41;
  JB->colptrs[9] = 45;
  JB->colptrs[10] = 53;
  JB->colptrs[11] = 60;
  JB->colptrs[12] = 68;
  JB->colptrs[13] = 75;
  JB->colptrs[14] = 82;
  return(0);
}


 int sx_enhancer_64(int Ns, realtype t, N_Vector x, N_Vector xdot,
  	int ip, N_Vector sx, N_Vector sxdot, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2)
{
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double *sxdot_tmp = N_VGetArrayPointer(sxdot);
  memset(sxdot_tmp,0,sizeof(double)*14);
  switch (plist[ip]) {
  case 0: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[13] - 2.0*p[3]*sx_tmp[0] + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[9] + ((pow(k[0],2))*x_tmp[4] - 2.0*(pow(k[0],3))*x_tmp[2]*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[6] + (pow(k[0],2))*x_tmp[6]*x_tmp[8])/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[3] = 2.0*p[1]*sx_tmp[11] + ((pow(k[0],2))*x_tmp[4] - 2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[6] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[8] + (pow(k[0],2))*x_tmp[6]*x_tmp[8])/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[8]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[1]*sx_tmp[11] - (1.0*((pow(k[0],3))*x_tmp[3]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],3))*x_tmp[2]*x_tmp[8] + (pow(k[0],3))*x_tmp[4]*x_tmp[6] + (pow(k[0],3))*x_tmp[4]*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*sx_tmp[4]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] - (1.0*sx_tmp[8]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[5] = p[0]*sx_tmp[4] - 2.0*p[1]*sx_tmp[5] + ((pow(k[0],2))*x_tmp[4] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[7] + (pow(k[0],2))*x_tmp[6]*x_tmp[8] + 2.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + 2.0*k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + 2.0*k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) - 6.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8] + 2.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[6] = p[1]*sx_tmp[7] - (1.0*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[4] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[7] = ((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8])/k[0] - 1.0*p[1]*sx_tmp[7] + k[0]*p[0]*sx_tmp[4] + k[0]*p[0]*sx_tmp[6]*x_tmp[8] + k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[7] - (1.0*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[4] - (1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/k[0] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8];
sxdot_tmp[9] = p[1]*sx_tmp[5] - (1.0*((pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],3))*x_tmp[2]*x_tmp[8] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[6] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[7] + (pow(k[0],2))*x_tmp[6]*x_tmp[8] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[10] = p[2]*sx_tmp[9] + p[1]*sx_tmp[13] - (1.0*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[4] + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*sx_tmp[10]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[11] = p[1]*sx_tmp[5] - (1.0*((pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],3))*x_tmp[3]*x_tmp[6] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[7] - 1.0*(pow(k[0],3))*x_tmp[4]*x_tmp[8] + (pow(k[0],2))*x_tmp[6]*x_tmp[8] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[6]*x_tmp[7]*x_tmp[8]))/(pow(k[0],2)) + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*sx_tmp[11]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8];
sxdot_tmp[12] = p[2]*sx_tmp[11] + p[1]*sx_tmp[13] - (1.0*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[4] + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - (1.0*sx_tmp[12]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[10]*x_tmp[8];
sxdot_tmp[13] = p[2]*sx_tmp[5] + (k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[4] + k[0]*x_tmp[8]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + k[0]*x_tmp[6]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6]*x_tmp[8])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[10]*x_tmp[8] + k[0]*p[0]*sx_tmp[12]*x_tmp[6];

  } break;

  case 1: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[13] - 2.0*p[3]*sx_tmp[0] + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[9] + (k[0]*x_tmp[7] + 2.0*(pow(k[0],2))*x_tmp[9])/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[3] = 2.0*p[1]*sx_tmp[11] + (k[0]*x_tmp[7] + 2.0*(pow(k[0],2))*x_tmp[11])/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[8]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[1]*sx_tmp[11] + (k[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[11])/(pow(k[0],2)) - (1.0*sx_tmp[4]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] - (1.0*sx_tmp[8]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[5] = p[0]*sx_tmp[4] - 2.0*p[1]*sx_tmp[5] + (k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[5])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8] + 2.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[6] = p[1]*sx_tmp[7] + x_tmp[7] - 1.0*k[0]*p[0]*sx_tmp[4] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[7] = k[0]*p[0]*sx_tmp[4] - 1.0*x_tmp[7] - 1.0*p[1]*sx_tmp[7] + k[0]*p[0]*sx_tmp[6]*x_tmp[8] + k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[7] + x_tmp[7] - 1.0*k[0]*p[0]*sx_tmp[4] - (1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/k[0] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8];
sxdot_tmp[9] = p[1]*sx_tmp[5] - (1.0*(k[0]*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[10] = p[2]*sx_tmp[9] + p[1]*sx_tmp[13] + x_tmp[13] - (1.0*sx_tmp[10]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[11] = p[1]*sx_tmp[5] - (1.0*(k[0]*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*sx_tmp[11]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8];
sxdot_tmp[12] = p[2]*sx_tmp[11] + p[1]*sx_tmp[13] + x_tmp[13] - (1.0*sx_tmp[12]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[10]*x_tmp[8];
sxdot_tmp[13] = p[2]*sx_tmp[5] - 1.0*x_tmp[13] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[10]*x_tmp[8] + k[0]*p[0]*sx_tmp[12]*x_tmp[6];

  } break;

  case 2: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[13] - 2.0*p[3]*sx_tmp[0] + (k[0]*x_tmp[7] + 2.0*(pow(k[0],2))*x_tmp[13])/(pow(k[0],2)) + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[1] + x_tmp[7];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[9] + (p[1]*sx_tmp[7])/k[0] + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[3] = 2.0*p[1]*sx_tmp[11] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[8]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[1]*sx_tmp[11] - (1.0*sx_tmp[4]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] - (1.0*sx_tmp[8]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[5] = p[0]*sx_tmp[4] - 2.0*p[1]*sx_tmp[5] + (sx_tmp[8]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8] + 2.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*k[0]*p[0]*sx_tmp[4] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[7] = k[0]*p[0]*sx_tmp[4] - 1.0*p[1]*sx_tmp[7] + k[0]*p[0]*sx_tmp[6]*x_tmp[8] + k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[7] - 1.0*k[0]*p[0]*sx_tmp[4] - (1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/k[0] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8];
sxdot_tmp[9] = p[1]*sx_tmp[5] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[10] = p[2]*sx_tmp[9] + p[1]*sx_tmp[13] + x_tmp[9] - (1.0*sx_tmp[10]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[11] = p[1]*sx_tmp[5] + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*sx_tmp[11]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8];
sxdot_tmp[12] = p[2]*sx_tmp[11] + p[1]*sx_tmp[13] + x_tmp[11] - (1.0*sx_tmp[12]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[10]*x_tmp[8];
sxdot_tmp[13] = p[2]*sx_tmp[5] + x_tmp[5] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[10]*x_tmp[8] + k[0]*p[0]*sx_tmp[12]*x_tmp[6];

  } break;

  case 3: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[13] - 2.0*p[3]*sx_tmp[0] + (k[0]*x_tmp[1] - 2.0*(pow(k[0],2))*x_tmp[0])/(pow(k[0],2)) + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[1] - 1.0*x_tmp[1];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[9] + (p[1]*sx_tmp[7])/k[0] + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[3] = 2.0*p[1]*sx_tmp[11] + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[8]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[1]*sx_tmp[11] - (1.0*sx_tmp[4]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] - (1.0*sx_tmp[8]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[5] = p[0]*sx_tmp[4] - 2.0*p[1]*sx_tmp[5] + (sx_tmp[8]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8] + 2.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*k[0]*p[0]*sx_tmp[4] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[7] = k[0]*p[0]*sx_tmp[4] - 1.0*p[1]*sx_tmp[7] + k[0]*p[0]*sx_tmp[6]*x_tmp[8] + k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[7] - 1.0*k[0]*p[0]*sx_tmp[4] - (1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/k[0] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8];
sxdot_tmp[9] = p[1]*sx_tmp[5] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[10] = p[2]*sx_tmp[9] + p[1]*sx_tmp[13] - 1.0*x_tmp[10] - (1.0*sx_tmp[10]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[11] = p[1]*sx_tmp[5] + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*sx_tmp[11]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8];
sxdot_tmp[12] = p[2]*sx_tmp[11] + p[1]*sx_tmp[13] - 1.0*x_tmp[12] - (1.0*sx_tmp[12]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[10]*x_tmp[8];
sxdot_tmp[13] = p[2]*sx_tmp[5] - 1.0*x_tmp[13] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[10]*x_tmp[8] + k[0]*p[0]*sx_tmp[12]*x_tmp[6];

  } break;

  case 4: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[13] - 2.0*p[3]*sx_tmp[0] + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[9] + (p[1]*sx_tmp[7])/k[0] + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[3] = 2.0*p[1]*sx_tmp[11] + 144.0/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[8]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[1]*sx_tmp[11] - (1.0*sx_tmp[4]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] - (1.0*sx_tmp[8]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[5] = p[0]*sx_tmp[4] - 2.0*p[1]*sx_tmp[5] + (sx_tmp[8]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8] + 2.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*k[0]*p[0]*sx_tmp[4] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[7] = k[0]*p[0]*sx_tmp[4] - 1.0*p[1]*sx_tmp[7] + k[0]*p[0]*sx_tmp[6]*x_tmp[8] + k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[7] + 12.0/k[0] - 1.0*k[0]*p[0]*sx_tmp[4] - (1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/k[0] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8];
sxdot_tmp[9] = p[1]*sx_tmp[5] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[10] = p[2]*sx_tmp[9] + p[1]*sx_tmp[13] - (1.0*sx_tmp[10]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[11] = p[1]*sx_tmp[5] + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*sx_tmp[11]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8];
sxdot_tmp[12] = p[2]*sx_tmp[11] + p[1]*sx_tmp[13] - (1.0*sx_tmp[12]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[10]*x_tmp[8];
sxdot_tmp[13] = p[2]*sx_tmp[5] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[10]*x_tmp[8] + k[0]*p[0]*sx_tmp[12]*x_tmp[6];

  } break;

  case 5: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[13] - 2.0*p[3]*sx_tmp[0] + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[7])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[7] - 1.0*p[3]*sx_tmp[1];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[9] + (p[1]*sx_tmp[7])/k[0] + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[3] = 2.0*p[1]*sx_tmp[11] + (k[0]*x_tmp[8] - 2.0*(pow(k[0],2))*x_tmp[3])/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] + (sx_tmp[8]*(k[0]*p[5] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8])*sx_tmp[6])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[1]*sx_tmp[11] - 1.0*x_tmp[4] - (1.0*sx_tmp[4]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (p[1]*sx_tmp[7])/k[0] - (1.0*sx_tmp[8]*((pow(k[0],3))*p[0]*x_tmp[2] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[8];
sxdot_tmp[5] = p[0]*sx_tmp[4] - 2.0*p[1]*sx_tmp[5] + (sx_tmp[8]*((pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8] + 2.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 1.0*k[0]*p[0]*sx_tmp[4] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[7] = k[0]*p[0]*sx_tmp[4] - 1.0*p[1]*sx_tmp[7] + k[0]*p[0]*sx_tmp[6]*x_tmp[8] + k[0]*p[0]*sx_tmp[8]*x_tmp[6];
sxdot_tmp[8] = p[1]*sx_tmp[7] - 1.0*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[4] - (1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[6])*sx_tmp[8])/k[0] - 1.0*k[0]*p[0]*sx_tmp[6]*x_tmp[8];
sxdot_tmp[9] = p[1]*sx_tmp[5] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[4])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[8] - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[6];
sxdot_tmp[10] = p[2]*sx_tmp[9] + p[1]*sx_tmp[13] - (1.0*sx_tmp[10]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[8]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[6];
sxdot_tmp[11] = p[1]*sx_tmp[5] - 1.0*x_tmp[11] + (sx_tmp[8]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*x_tmp[6]*x_tmp[7]) + (pow(k[0],3))*p[0]*x_tmp[4] - 1.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) + (sx_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*x_tmp[7]*x_tmp[8]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]) - 1.0*k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8])*sx_tmp[4])/(pow(k[0],2)) - (1.0*sx_tmp[11]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[8];
sxdot_tmp[12] = p[2]*sx_tmp[11] + p[1]*sx_tmp[13] - 1.0*x_tmp[12] - (1.0*sx_tmp[12]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) + (p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[10]*x_tmp[8];
sxdot_tmp[13] = p[2]*sx_tmp[5] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[10] - 1.0*k[0]*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[1]*x_tmp[6]) + (pow(k[0],3))*x_tmp[1]*x_tmp[6]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[8])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[12] - 1.0*k[0]*((pow(k[0],2))*x_tmp[12] + (pow(k[0],2))*x_tmp[1]*x_tmp[8]) + (pow(k[0],3))*x_tmp[1]*x_tmp[8]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[6])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[4] - 1.0*k[0]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[6]*x_tmp[8]) + (pow(k[0],3))*x_tmp[6]*x_tmp[8]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[10]*x_tmp[8] + k[0]*p[0]*sx_tmp[12]*x_tmp[6];

  } break;

  }
 for (ix=0; ix<14; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_enhancer_64(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  memset(sx0_tmp,0,sizeof(double)*14);
  switch (ip) {
  }

  return;
}


void y_enhancer_64(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*1];
y[it+nt*1] = x[it+nt*0];
    
    return;
}


void dydp_enhancer_64(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }
  
  return;
}


void dydx_enhancer_64(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*28);
dydx[1] = 1.0;
dydx[2] = 1.0;
  
  return;
}


void sy_enhancer_64(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(0+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(0+np*nx)];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(0+np*nx)];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(0+np*nx)];

  } break;

  case 4: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(0+np*nx)];

  } break;

  case 5: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(1+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(0+np*nx)];

  } break;

  }
  
  return;
}
int root_enhancer_64(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_enhancer_64(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double dr_dp;
  switch (ip) {
  }
  return(dr_dp);
}
double s2root_enhancer_64(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double ddr_dpdp;
  switch (ip) {
  }
  return(ddr_dpdp);
}
double srootval_enhancer_64(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double dg_dp;
  switch (ip) {
  }
  return(dg_dp);
}
double s2rootval_enhancer_64(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp = N_VGetArrayPointer(sx);
  double ddg_dpdp;
  switch (ip) {
  }
  return(ddg_dpdp);
}
void deltadisc_enhancer_64(double t, int idisc, N_Vector x, void *user_data){
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double deltadisc[14];
  memset(deltadisc,0,sizeof(double)*14);
  for(ix = 0; ix<14;ix++){;
  x_tmp[ix] += deltadisc[ix];
  };
}
void sdeltadisc_enhancer_64(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data){
  int ix;
  int ip;
  UserData data = (UserData) user_data;
  int *plist = data->plist;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *sx_tmp;
  int np = *data->np;
  double deltadisc[14];
  double *sdeltadisc;
  memset(deltadisc,0,sizeof(double)*14);
  sdeltadisc = mxMalloc(sizeof(double)*14*np);
  memset(sdeltadisc,0,sizeof(double)*14*np);
  for (ip=0; ip<np; ip++) {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
     switch (plist[ip]) {
     }
  }
  for(ip = 0; ip<np;ip++){
      sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
      for(ix = 0; ix<14;ix++){
      sx_tmp[ix] += sdeltadisc[plist[ip]+np*ix];
     }
  }
  for(ix = 0; ix<14;ix++){
  x_tmp[ix] += deltadisc[ix];
  };
 mxFree(sdeltadisc);
}


void dxdotdp_enhancer_64(double t, int nt, int it, double *dxdotdp, double *p, double *k, double *u, double *x, int *plist, int np, int nx){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
dxdotdp[(2+ip*nx)] = ((pow(k[0],2))*x[it+nt*4] - 2.0*(pow(k[0],3))*x[it+nt*2]*x[it+nt*8] - 2.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*6] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8])/(pow(k[0],2));
dxdotdp[(3+ip*nx)] = ((pow(k[0],2))*x[it+nt*4] - 2.0*(pow(k[0],3))*x[it+nt*3]*x[it+nt*6] - 2.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*8] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8])/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = -(1.0*((pow(k[0],3))*x[it+nt*3]*x[it+nt*6] - 1.0*(pow(k[0],2))*x[it+nt*4] + (pow(k[0],3))*x[it+nt*2]*x[it+nt*8] + (pow(k[0],3))*x[it+nt*4]*x[it+nt*6] + (pow(k[0],3))*x[it+nt*4]*x[it+nt*8] - 1.0*(pow(k[0],2))*x[it+nt*6]*x[it+nt*8]))/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = ((pow(k[0],2))*x[it+nt*4] - 2.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*7] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8] + 2.0*k[0]*x[it+nt*7]*((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8]) + 2.0*k[0]*x[it+nt*8]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*7]) + 2.0*k[0]*x[it+nt*6]*((pow(k[0],2))*x[it+nt*11] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*8]) - 6.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*7]*x[it+nt*8])/(pow(k[0],2));
dxdotdp[(6+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8]))/k[0];
dxdotdp[(7+ip*nx)] = ((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8])/k[0];
dxdotdp[(8+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8]))/k[0];
dxdotdp[(9+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*4] - 1.0*(pow(k[0],3))*x[it+nt*2]*x[it+nt*8] - 1.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*6] - 1.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*7] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8] + k[0]*x[it+nt*7]*((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8]) + k[0]*x[it+nt*8]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*7]) + k[0]*x[it+nt*6]*((pow(k[0],2))*x[it+nt*11] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*8]) - 3.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*7]*x[it+nt*8]))/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = -(1.0*(k[0]*x[it+nt*1]*((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8]) - 1.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*4] + k[0]*x[it+nt*8]*((pow(k[0],2))*x[it+nt*10] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*6]) + k[0]*x[it+nt*6]*((pow(k[0],2))*x[it+nt*12] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*8]) - 3.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*6]*x[it+nt*8]))/(pow(k[0],2));
dxdotdp[(11+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*4] - 1.0*(pow(k[0],3))*x[it+nt*3]*x[it+nt*6] - 1.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*7] - 1.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*8] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8] + k[0]*x[it+nt*7]*((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8]) + k[0]*x[it+nt*8]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*7]) + k[0]*x[it+nt*6]*((pow(k[0],2))*x[it+nt*11] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*8]) - 3.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*7]*x[it+nt*8]))/(pow(k[0],2));
dxdotdp[(12+ip*nx)] = -(1.0*(k[0]*x[it+nt*1]*((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8]) - 1.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*4] + k[0]*x[it+nt*8]*((pow(k[0],2))*x[it+nt*10] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*6]) + k[0]*x[it+nt*6]*((pow(k[0],2))*x[it+nt*12] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*8]) - 3.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*6]*x[it+nt*8]))/(pow(k[0],2));
dxdotdp[(13+ip*nx)] = (k[0]*x[it+nt*1]*((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*6]*x[it+nt*8]) - 1.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*4] + k[0]*x[it+nt*8]*((pow(k[0],2))*x[it+nt*10] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*6]) + k[0]*x[it+nt*6]*((pow(k[0],2))*x[it+nt*12] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*8]) - 3.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*6]*x[it+nt*8])/(pow(k[0],2));

  } break;

  case 1: {
dxdotdp[(2+ip*nx)] = (k[0]*x[it+nt*7] + 2.0*(pow(k[0],2))*x[it+nt*9])/(pow(k[0],2));
dxdotdp[(3+ip*nx)] = (k[0]*x[it+nt*7] + 2.0*(pow(k[0],2))*x[it+nt*11])/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = (k[0]*x[it+nt*7] + (pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*11])/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = (k[0]*x[it+nt*7] - 2.0*(pow(k[0],2))*x[it+nt*5])/(pow(k[0],2));
dxdotdp[(6+ip*nx)] = x[it+nt*7];
dxdotdp[(7+ip*nx)] = -1.0*x[it+nt*7];
dxdotdp[(8+ip*nx)] = x[it+nt*7];
dxdotdp[(9+ip*nx)] = -(1.0*(k[0]*x[it+nt*7] - 1.0*(pow(k[0],2))*x[it+nt*5] + (pow(k[0],2))*x[it+nt*9]))/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = x[it+nt*13];
dxdotdp[(11+ip*nx)] = -(1.0*(k[0]*x[it+nt*7] - 1.0*(pow(k[0],2))*x[it+nt*5] + (pow(k[0],2))*x[it+nt*11]))/(pow(k[0],2));
dxdotdp[(12+ip*nx)] = x[it+nt*13];
dxdotdp[(13+ip*nx)] = -1.0*x[it+nt*13];

  } break;

  case 2: {
dxdotdp[(0+ip*nx)] = (k[0]*x[it+nt*7] + 2.0*(pow(k[0],2))*x[it+nt*13])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = x[it+nt*7];
dxdotdp[(10+ip*nx)] = x[it+nt*9];
dxdotdp[(12+ip*nx)] = x[it+nt*11];
dxdotdp[(13+ip*nx)] = x[it+nt*5];

  } break;

  case 3: {
dxdotdp[(0+ip*nx)] = (k[0]*x[it+nt*1] - 2.0*(pow(k[0],2))*x[it+nt*0])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = -1.0*x[it+nt*1];
dxdotdp[(10+ip*nx)] = -1.0*x[it+nt*10];
dxdotdp[(12+ip*nx)] = -1.0*x[it+nt*12];
dxdotdp[(13+ip*nx)] = -1.0*x[it+nt*13];

  } break;

  case 4: {
dxdotdp[(3+ip*nx)] = 144.0/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = 12.0/k[0];

  } break;

  case 5: {
dxdotdp[(3+ip*nx)] = (k[0]*x[it+nt*8] - 2.0*(pow(k[0],2))*x[it+nt*3])/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = -1.0*x[it+nt*4];
dxdotdp[(8+ip*nx)] = -1.0*x[it+nt*8];
dxdotdp[(11+ip*nx)] = -1.0*x[it+nt*11];
dxdotdp[(12+ip*nx)] = -1.0*x[it+nt*12];

  } break;

  }
  }
  
  return;
}
