#include "enhancer_321.h"
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


 int xdot_enhancer_321(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int ix;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  memset(xdot_tmp,0,sizeof(double)*20);
xdot_tmp[0] = (2.0*(pow(k[0],2))*p[2]*x_tmp[17] - 2.0*(pow(k[0],2))*p[3]*x_tmp[0] + 4.0*(pow(k[0],2))*p[2]*x_tmp[18] + k[0]*p[3]*x_tmp[1] + k[0]*p[2]*x_tmp[11] + 2.0*k[0]*p[2]*x_tmp[12])/(pow(k[0],2));
xdot_tmp[1] = (k[0]*p[2]*x_tmp[11] - 1.0*k[0]*p[3]*x_tmp[1] + 2.0*k[0]*p[2]*x_tmp[12])/k[0];
xdot_tmp[2] = ((pow(k[0],2))*p[0]*x_tmp[8] - 4.0*(pow(k[0],2))*p[1]*x_tmp[2] - 2.0*p[0]*((pow(k[0],3))*x_tmp[8]*x_tmp[12] + (pow(k[0],3))*x_tmp[9]*x_tmp[11] + (pow(k[0],3))*x_tmp[13]*x_tmp[16] - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + 3.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]) + 2.0*k[0]*p[1]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]*x_tmp[11] + (pow(k[0],2))*p[0]*x_tmp[11]*x_tmp[13] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[16])/(pow(k[0],2));
xdot_tmp[3] = (16.0*p[4] + 2.0*(pow(k[0],2))*p[0]*x_tmp[6] + (pow(k[0],2))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[5]*x_tmp[3] + 2.0*(pow(k[0],2))*p[1]*x_tmp[8] + 4.0*(pow(k[0],2))*p[1]*x_tmp[9] + k[0]*p[1]*x_tmp[11] + 2.0*k[0]*p[1]*x_tmp[12] + k[0]*p[5]*x_tmp[13] - 4.0*(pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[11] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[13] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[13] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*p[0]*x_tmp[11]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[4] = (2.0*(pow(k[0],2))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[1]*x_tmp[5] + k[0]*p[1]*x_tmp[11] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[10] - 4.0*(pow(k[0],3))*p[0]*x_tmp[4]*x_tmp[13] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[5] = -(1.0*(2.0*(pow(k[0],2))*p[0]*x_tmp[6] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*x_tmp[5]*x_tmp[13] + (pow(k[0],3))*x_tmp[8]*x_tmp[10] - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]) + (pow(k[0],2))*p[1]*x_tmp[5] - 2.0*(pow(k[0],2))*p[1]*x_tmp[7] - 1.0*(pow(k[0],2))*p[1]*x_tmp[15] + k[0]*p[1]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4]*x_tmp[13] + (pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[11] + 3.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[13] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[10] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
xdot_tmp[6] = (p[0]*((pow(k[0],3))*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*x_tmp[5]*x_tmp[13] + (pow(k[0],3))*x_tmp[8]*x_tmp[10] - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[6] + (pow(k[0],2))*p[1]*x_tmp[5] + 2.0*(pow(k[0],2))*p[1]*x_tmp[7] + (pow(k[0],2))*p[1]*x_tmp[8] - 1.0*(pow(k[0],2))*p[5]*x_tmp[6] + k[0]*p[1]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[4]*x_tmp[13] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[13] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[13] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[7] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*x_tmp[5]*x_tmp[13] + (pow(k[0],3))*x_tmp[8]*x_tmp[10] - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[6]*x_tmp[12] + (pow(k[0],3))*x_tmp[9]*x_tmp[10] + (pow(k[0],3))*x_tmp[7]*x_tmp[13] - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[1]*x_tmp[7] - 1.0*(pow(k[0],2))*p[1]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[13] + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]*x_tmp[13]))/(pow(k[0],2));
xdot_tmp[8] = -(1.0*(2.0*(pow(k[0],2))*p[0]*x_tmp[6] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*x_tmp[5]*x_tmp[13] + (pow(k[0],3))*x_tmp[8]*x_tmp[10] - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] + (pow(k[0],2))*p[1]*x_tmp[8] - 2.0*(pow(k[0],2))*p[1]*x_tmp[9] + (pow(k[0],2))*p[5]*x_tmp[8] - 1.0*(pow(k[0],2))*p[1]*x_tmp[15] - 2.0*(pow(k[0],2))*p[1]*x_tmp[16] + k[0]*p[1]*x_tmp[11] - 2.0*k[0]*p[1]*x_tmp[12] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[13] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6]*x_tmp[13] + (pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[13] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10]*x_tmp[13] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11]*x_tmp[13] + (pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[15]))/(pow(k[0],2));
xdot_tmp[9] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[8] - 1.0*p[0]*((pow(k[0],3))*x_tmp[8]*x_tmp[12] + (pow(k[0],3))*x_tmp[9]*x_tmp[11] + (pow(k[0],3))*x_tmp[13]*x_tmp[16] - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + 3.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],2))*p[1]*x_tmp[2] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6]*x_tmp[12] + (pow(k[0],3))*x_tmp[9]*x_tmp[10] + (pow(k[0],3))*x_tmp[7]*x_tmp[13] - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[1]*x_tmp[9] + (pow(k[0],2))*p[5]*x_tmp[9] - 1.0*(pow(k[0],2))*p[1]*x_tmp[16] + 2.0*k[0]*p[1]*x_tmp[12] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]*x_tmp[13] + (pow(k[0],3))*p[0]*x_tmp[9]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[13] + (pow(k[0],2))*p[0]*x_tmp[11]*x_tmp[13] + (pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[16]))/(pow(k[0],2));
xdot_tmp[10] = -(1.0*(2.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*k[0]*p[1]*x_tmp[11] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10]*x_tmp[13]))/k[0];
xdot_tmp[11] = (2.0*(pow(k[0],2))*p[0]*x_tmp[6] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8] - 1.0*k[0]*p[1]*x_tmp[11] + 2.0*k[0]*p[1]*x_tmp[12] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10]*x_tmp[13] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11]*x_tmp[13])/k[0];
xdot_tmp[12] = ((pow(k[0],2))*p[0]*x_tmp[8] - 2.0*k[0]*p[1]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[11]*x_tmp[13])/k[0];
xdot_tmp[13] = -(1.0*(2.0*(pow(k[0],2))*p[0]*x_tmp[6] - 4.0*p[4] + (pow(k[0],2))*p[0]*x_tmp[8] - 1.0*k[0]*p[1]*x_tmp[11] - 2.0*k[0]*p[1]*x_tmp[12] + k[0]*p[5]*x_tmp[13] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*p[0]*x_tmp[11]*x_tmp[13]))/k[0];
xdot_tmp[14] = (2.0*p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[6] + (pow(k[0],3))*x_tmp[13]*x_tmp[14] + (pow(k[0],3))*x_tmp[10]*x_tmp[19] - 1.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13]) + (pow(k[0],2))*p[2]*x_tmp[5] + 2.0*(pow(k[0],2))*p[2]*x_tmp[7] - 1.0*(pow(k[0],2))*p[3]*x_tmp[14] + (pow(k[0],2))*p[1]*x_tmp[17] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10]*x_tmp[19])/(pow(k[0],2));
xdot_tmp[15] = (2.0*(pow(k[0],2))*p[0]*x_tmp[6] - 4.0*p[0]*((pow(k[0],3))*x_tmp[6]*x_tmp[11] + (pow(k[0],3))*x_tmp[5]*x_tmp[13] + (pow(k[0],3))*x_tmp[8]*x_tmp[10] - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[1]*x_tmp[15] + 4.0*(pow(k[0],2))*p[1]*x_tmp[16] + k[0]*p[1]*x_tmp[11] + 2.0*k[0]*p[1]*x_tmp[12] + 4.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[13] + 4.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[11] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*p[0]*x_tmp[11]*x_tmp[13] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[15])/(pow(k[0],2));
xdot_tmp[16] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6]*x_tmp[12] + (pow(k[0],3))*x_tmp[9]*x_tmp[10] + (pow(k[0],3))*x_tmp[7]*x_tmp[13] - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8]*x_tmp[12] + (pow(k[0],3))*x_tmp[9]*x_tmp[11] + (pow(k[0],3))*x_tmp[13]*x_tmp[16] - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + 3.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],2))*p[1]*x_tmp[2] + (pow(k[0],2))*p[0]*x_tmp[8] + 3.0*(pow(k[0],2))*p[1]*x_tmp[16] + 2.0*k[0]*p[1]*x_tmp[12] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[9]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[7]*x_tmp[13] + (pow(k[0],3))*p[0]*x_tmp[9]*x_tmp[11] + (pow(k[0],2))*p[0]*x_tmp[11]*x_tmp[13] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[16]))/(pow(k[0],2));
xdot_tmp[17] = (p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[8] + (pow(k[0],3))*x_tmp[11]*x_tmp[19] + (pow(k[0],3))*x_tmp[13]*x_tmp[17] - 1.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[6] + (pow(k[0],3))*x_tmp[13]*x_tmp[14] + (pow(k[0],3))*x_tmp[10]*x_tmp[19] - 1.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13]) + (pow(k[0],2))*p[2]*x_tmp[15] - 1.0*(pow(k[0],2))*p[1]*x_tmp[17] + 2.0*(pow(k[0],2))*p[2]*x_tmp[16] + 2.0*(pow(k[0],2))*p[1]*x_tmp[18] - 1.0*(pow(k[0],2))*p[3]*x_tmp[17] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[14] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10]*x_tmp[19] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]*x_tmp[19] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[17])/(pow(k[0],2));
xdot_tmp[18] = (2.0*(pow(k[0],2))*p[2]*x_tmp[2] - 1.0*p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[8] + (pow(k[0],3))*x_tmp[11]*x_tmp[19] + (pow(k[0],3))*x_tmp[13]*x_tmp[17] - 1.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]) + (pow(k[0],2))*p[2]*x_tmp[16] - 2.0*(pow(k[0],2))*p[1]*x_tmp[18] - 1.0*(pow(k[0],2))*p[3]*x_tmp[18] + (pow(k[0],3))*p[0]*x_tmp[11]*x_tmp[19] + (pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[17])/(pow(k[0],2));
xdot_tmp[19] = (2.0*p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[6] + (pow(k[0],3))*x_tmp[13]*x_tmp[14] + (pow(k[0],3))*x_tmp[10]*x_tmp[19] - 1.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) - 1.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[8] + (pow(k[0],3))*x_tmp[11]*x_tmp[19] + (pow(k[0],3))*x_tmp[13]*x_tmp[17] - 1.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) - 1.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]) + (pow(k[0],2))*p[2]*x_tmp[8] + 2.0*(pow(k[0],2))*p[2]*x_tmp[9] + (pow(k[0],2))*p[1]*x_tmp[17] + 2.0*(pow(k[0],2))*p[1]*x_tmp[18] - 1.0*(pow(k[0],2))*p[3]*x_tmp[19] - 1.0*(pow(k[0],2))*p[5]*x_tmp[19] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10]*x_tmp[19] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]*x_tmp[19] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13]*x_tmp[17])/(pow(k[0],2));

  for (ix=0; ix<20; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_enhancer_321(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*40);
xBdot_tmp[0] = 2.0*p[3]*xB_tmp[0];
xBdot_tmp[1] = p[3]*xB_tmp[1] - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*xB_tmp[19])/(pow(k[0],2)) - (1.0*p[3]*xB_tmp[0])/k[0] + (xB_tmp[17]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) - (2.0*p[0]*xB_tmp[14]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + (p[0]*xB_tmp[18]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
xBdot_tmp[2] = 4.0*p[1]*xB_tmp[2] - 2.0*p[1]*xB_tmp[9] - 2.0*p[1]*xB_tmp[16] - 2.0*p[2]*xB_tmp[18];
xBdot_tmp[3] = (xB_tmp[3]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[6] - 1.0*k[0]*p[0]*x_tmp[11]*xB_tmp[9];
xBdot_tmp[4] = 4.0*k[0]*p[0]*x_tmp[13]*xB_tmp[4] - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[5] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[6];
xBdot_tmp[5] = (((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[5])/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[14] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[6])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[4] - 1.0*k[0]*p[0]*x_tmp[13]*xB_tmp[7] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[8] - 4.0*k[0]*p[0]*x_tmp[13]*xB_tmp[15];
xBdot_tmp[6] = (xB_tmp[6]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*p[0]*xB_tmp[15] - (1.0*(2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*xB_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[3])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[0]*xB_tmp[10] - 2.0*k[0]*p[0]*xB_tmp[11] + 2.0*k[0]*p[0]*xB_tmp[13] + (xB_tmp[5]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*x_tmp[11]*xB_tmp[7];
xBdot_tmp[7] = ((2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[7])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[6] - 2.0*p[2]*xB_tmp[14] - 2.0*p[1]*xB_tmp[5] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[9] - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[16];
xBdot_tmp[8] = (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[9])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[6] - 1.0*p[2]*xB_tmp[19] - (1.0*xB_tmp[3]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[2] + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[16])/(pow(k[0],2)) + k[0]*p[0]*xB_tmp[11] - 1.0*k[0]*p[0]*xB_tmp[12] + k[0]*p[0]*xB_tmp[13] + (xB_tmp[8]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[15]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[5];
xBdot_tmp[9] = (xB_tmp[9]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[8] - 2.0*p[2]*xB_tmp[19] - 4.0*p[1]*xB_tmp[3] - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*x_tmp[11]*xB_tmp[2] + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[7];
xBdot_tmp[10] = (xB_tmp[7]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[15]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*xB_tmp[3])/(pow(k[0],2)) + ((4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*xB_tmp[4])/(pow(k[0],2)) - (1.0*xB_tmp[8]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*xB_tmp[9])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*xB_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[17])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[19])/(pow(k[0],2)) - (1.0*xB_tmp[5]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[10] - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[11] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[13];
xBdot_tmp[11] = (xB_tmp[8]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[10] - (1.0*xB_tmp[2]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*xB_tmp[13])/k[0] - (1.0*p[2]*xB_tmp[0])/k[0] - (1.0*p[1]*xB_tmp[4])/k[0] - 1.0*p[2]*xB_tmp[1] - (1.0*xB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[15]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*xB_tmp[7])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[17])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[19])/(pow(k[0],2)) + (xB_tmp[5]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*xB_tmp[11])/k[0] - 1.0*k[0]*p[0]*x_tmp[13]*xB_tmp[12];
xBdot_tmp[12] = 2.0*p[1]*xB_tmp[12] - 2.0*p[1]*xB_tmp[11] - 2.0*p[2]*xB_tmp[1] - 2.0*p[1]*xB_tmp[13] - (2.0*p[2]*xB_tmp[0])/k[0] - (2.0*p[1]*xB_tmp[3])/k[0] - (2.0*p[1]*xB_tmp[8])/k[0] - (2.0*p[1]*xB_tmp[15])/k[0] - (1.0*xB_tmp[9]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) + (xB_tmp[16]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*xB_tmp[2])/(pow(k[0],2)) - (2.0*p[0]*xB_tmp[7]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
xBdot_tmp[13] = (xB_tmp[6]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[15]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[5]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[2]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (xB_tmp[13]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - (1.0*xB_tmp[9]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[16]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*xB_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*xB_tmp[11])/k[0] + (xB_tmp[8]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*xB_tmp[18])/(pow(k[0],2)) + (xB_tmp[7]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (xB_tmp[17]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[19]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*xB_tmp[14])/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[10] - 1.0*k[0]*p[0]*x_tmp[11]*xB_tmp[12];
xBdot_tmp[14] = (((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[14])/(pow(k[0],2)) - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[17] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[19];
xBdot_tmp[15] = ((2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[15])/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[17] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[8])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[5] - 1.0*k[0]*p[0]*x_tmp[13]*xB_tmp[16];
xBdot_tmp[16] = ((3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[16])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[8] - 4.0*p[1]*xB_tmp[15] - 2.0*p[2]*xB_tmp[17] - 1.0*p[2]*xB_tmp[18] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[9])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[7] - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[2];
xBdot_tmp[17] = (xB_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[14] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[19])/(pow(k[0],2)) - 2.0*p[2]*xB_tmp[0] - 1.0*k[0]*p[0]*x_tmp[13]*xB_tmp[18];
xBdot_tmp[18] = ((2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[18])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[17] - 2.0*p[1]*xB_tmp[19] - 4.0*p[2]*xB_tmp[0];
xBdot_tmp[19] = (xB_tmp[19]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[17])/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[14] - 1.0*k[0]*p[0]*x_tmp[11]*xB_tmp[18];
xBdot_tmp[20] = 2.0*p[3]*xB_tmp[20];
xBdot_tmp[21] = p[3]*xB_tmp[21] - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*xB_tmp[39])/(pow(k[0],2)) - (1.0*p[3]*xB_tmp[20])/k[0] + (xB_tmp[37]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) - (2.0*p[0]*xB_tmp[34]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + (p[0]*xB_tmp[38]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
xBdot_tmp[22] = 4.0*p[1]*xB_tmp[22] - 2.0*p[1]*xB_tmp[29] - 2.0*p[1]*xB_tmp[36] - 2.0*p[2]*xB_tmp[38];
xBdot_tmp[23] = (xB_tmp[23]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[28])/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[26] - 1.0*k[0]*p[0]*x_tmp[11]*xB_tmp[29];
xBdot_tmp[24] = 4.0*k[0]*p[0]*x_tmp[13]*xB_tmp[24] - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[25] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[26];
xBdot_tmp[25] = (((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[25])/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[34] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[26])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[24] - 1.0*k[0]*p[0]*x_tmp[13]*xB_tmp[27] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[28] - 4.0*k[0]*p[0]*x_tmp[13]*xB_tmp[35];
xBdot_tmp[26] = (xB_tmp[26]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*p[0]*xB_tmp[35] - (1.0*(2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*xB_tmp[24])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[23])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[28])/(pow(k[0],2)) + 2.0*k[0]*p[0]*xB_tmp[30] - 2.0*k[0]*p[0]*xB_tmp[31] + 2.0*k[0]*p[0]*xB_tmp[33] + (xB_tmp[25]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*x_tmp[11]*xB_tmp[27];
xBdot_tmp[27] = ((2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[27])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[26] - 2.0*p[2]*xB_tmp[34] - 2.0*p[1]*xB_tmp[25] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[29] - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[36];
xBdot_tmp[28] = (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[29])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[26] - 1.0*p[2]*xB_tmp[39] - (1.0*xB_tmp[23]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[22] + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[36])/(pow(k[0],2)) + k[0]*p[0]*xB_tmp[31] - 1.0*k[0]*p[0]*xB_tmp[32] + k[0]*p[0]*xB_tmp[33] + (xB_tmp[28]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[35]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[25];
xBdot_tmp[29] = (xB_tmp[29]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[28] - 2.0*p[2]*xB_tmp[39] - 4.0*p[1]*xB_tmp[23] - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[36])/(pow(k[0],2)) - 2.0*k[0]*p[0]*x_tmp[11]*xB_tmp[22] + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[27];
xBdot_tmp[30] = (xB_tmp[27]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[35]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[26]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*xB_tmp[23])/(pow(k[0],2)) + ((4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*xB_tmp[24])/(pow(k[0],2)) - (1.0*xB_tmp[28]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*xB_tmp[29])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*xB_tmp[36])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[34])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[37])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[39])/(pow(k[0],2)) - (1.0*xB_tmp[25]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[30] - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[31] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[33];
xBdot_tmp[31] = (xB_tmp[28]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[30] - (1.0*xB_tmp[22]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*xB_tmp[33])/k[0] - (1.0*p[2]*xB_tmp[20])/k[0] - (1.0*p[1]*xB_tmp[24])/k[0] - 1.0*p[2]*xB_tmp[21] - (1.0*xB_tmp[29]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[36]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[35]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[26]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*xB_tmp[27])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[37])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[38])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*xB_tmp[39])/(pow(k[0],2)) + (xB_tmp[25]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*xB_tmp[31])/k[0] - 1.0*k[0]*p[0]*x_tmp[13]*xB_tmp[32];
xBdot_tmp[32] = 2.0*p[1]*xB_tmp[32] - 2.0*p[1]*xB_tmp[31] - 2.0*p[2]*xB_tmp[21] - 2.0*p[1]*xB_tmp[33] - (2.0*p[2]*xB_tmp[20])/k[0] - (2.0*p[1]*xB_tmp[23])/k[0] - (2.0*p[1]*xB_tmp[28])/k[0] - (2.0*p[1]*xB_tmp[35])/k[0] - (1.0*xB_tmp[29]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) + (xB_tmp[36]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*xB_tmp[22])/(pow(k[0],2)) - (2.0*p[0]*xB_tmp[27]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
xBdot_tmp[33] = (xB_tmp[26]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[35]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[25]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[22]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (xB_tmp[33]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - (1.0*xB_tmp[29]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[36]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*xB_tmp[24])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*xB_tmp[31])/k[0] + (xB_tmp[28]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*xB_tmp[38])/(pow(k[0],2)) + (xB_tmp[27]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (xB_tmp[37]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[39]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*xB_tmp[34])/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[30] - 1.0*k[0]*p[0]*x_tmp[11]*xB_tmp[32];
xBdot_tmp[34] = (((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[34])/(pow(k[0],2)) - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[37] + 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[39];
xBdot_tmp[35] = ((2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[35])/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[37] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[28])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[25] - 1.0*k[0]*p[0]*x_tmp[13]*xB_tmp[36];
xBdot_tmp[36] = ((3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[36])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[28] - 4.0*p[1]*xB_tmp[35] - 2.0*p[2]*xB_tmp[37] - 1.0*p[2]*xB_tmp[38] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[29])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[27] - 2.0*k[0]*p[0]*x_tmp[13]*xB_tmp[22];
xBdot_tmp[37] = (xB_tmp[37]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[34] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*xB_tmp[39])/(pow(k[0],2)) - 2.0*p[2]*xB_tmp[20] - 1.0*k[0]*p[0]*x_tmp[13]*xB_tmp[38];
xBdot_tmp[38] = ((2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[38])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[37] - 2.0*p[1]*xB_tmp[39] - 4.0*p[2]*xB_tmp[20];
xBdot_tmp[39] = (xB_tmp[39]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[37])/(pow(k[0],2)) + 2.0*k[0]*p[0]*x_tmp[10]*xB_tmp[34] - 1.0*k[0]*p[0]*x_tmp[11]*xB_tmp[38];

  for (ixB=0; ixB<40; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_enhancer_321(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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
qBdot_tmp[0+ip*ny] = (xB_tmp[5]*(2.0*(pow(k[0],2))*x_tmp[6] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[10] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[10] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + 3.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + 3.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 3.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 9.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[17]*(2.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6] - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[8] - 2.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) - 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + 6.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[19]*(2.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6] + k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) - 6.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[7]*(2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[10] + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) - 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) - 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13] + 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*(2.0*(pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[8] - 4.0*(pow(k[0],3))*x_tmp[3]*x_tmp[10] - 2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[11] - 4.0*(pow(k[0],3))*x_tmp[6]*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[14]*(2.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 2.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6] + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) - 6.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[8]*(2.0*(pow(k[0],2))*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[10] + (pow(k[0],3))*x_tmp[3]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[13] + (pow(k[0],3))*x_tmp[8]*x_tmp[11] + (pow(k[0],3))*x_tmp[8]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[11]*x_tmp[13] + (pow(k[0],3))*x_tmp[13]*x_tmp[15] + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + 2.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]))/k[0] + (xB_tmp[13]*(2.0*(pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]))/k[0] + (xB_tmp[6]*(2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[10] - 2.0*(pow(k[0],2))*x_tmp[6] + 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[10] + 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[13] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[10] + 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[18]*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[8] + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[9]*((pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],3))*x_tmp[3]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] + 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) - 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[16]*((pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[11] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[12] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[15] - 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) - 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[4]*(2.0*(pow(k[0],2))*x_tmp[6] - 4.0*(pow(k[0],3))*x_tmp[6]*x_tmp[10] - 4.0*(pow(k[0],3))*x_tmp[4]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*x_tmp[6] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13])*xB_tmp[10])/k[0] - (1.0*xB_tmp[2]*((pow(k[0],2))*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[12] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] + 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) - 6.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[11]*(2.0*(pow(k[0],2))*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[11]*x_tmp[13]))/k[0] - (1.0*xB_tmp[15]*(2.0*(pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[8] - 4.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[13]*x_tmp[15] + 4.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + 4.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 4.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 12.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = (xB_tmp[25]*(2.0*(pow(k[0],2))*x_tmp[6] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[10] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[10] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + 3.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + 3.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 3.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 9.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[37]*(2.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6] - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[8] - 2.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) - 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + 6.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[39]*(2.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6] + k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) - 6.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[27]*(2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[10] + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) - 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) - 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13] + 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*(2.0*(pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[8] - 4.0*(pow(k[0],3))*x_tmp[3]*x_tmp[10] - 2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[11] - 4.0*(pow(k[0],3))*x_tmp[6]*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[34]*(2.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 2.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6] + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) - 6.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[28]*(2.0*(pow(k[0],2))*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[10] + (pow(k[0],3))*x_tmp[3]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[13] + (pow(k[0],3))*x_tmp[8]*x_tmp[11] + (pow(k[0],3))*x_tmp[8]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[11]*x_tmp[13] + (pow(k[0],3))*x_tmp[13]*x_tmp[15] + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + 2.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[32]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]))/k[0] + (xB_tmp[33]*(2.0*(pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]))/k[0] + (xB_tmp[26]*(2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[10] - 2.0*(pow(k[0],2))*x_tmp[6] + 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[10] + 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[13] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[10] + 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[38]*(k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[8] + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[29]*((pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],3))*x_tmp[3]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] + 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) - 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) + (xB_tmp[36]*((pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[11] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[12] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[15] - 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) - 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[24]*(2.0*(pow(k[0],2))*x_tmp[6] - 4.0*(pow(k[0],3))*x_tmp[6]*x_tmp[10] - 4.0*(pow(k[0],3))*x_tmp[4]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*x_tmp[6] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13])*xB_tmp[30])/k[0] - (1.0*xB_tmp[22]*((pow(k[0],2))*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[12] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] + 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) - 6.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) - (1.0*xB_tmp[31]*(2.0*(pow(k[0],2))*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[11]*x_tmp[13]))/k[0] - (1.0*xB_tmp[35]*(2.0*(pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[8] - 4.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[13]*x_tmp[15] + 4.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + 4.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 4.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 12.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));

  } break;

  case 1: {
qBdot_tmp[0+ip*ny] = 2.0*x_tmp[12]*xB_tmp[12] - 1.0*x_tmp[11]*xB_tmp[10] - 1.0*x_tmp[17]*xB_tmp[14] + 2.0*x_tmp[18]*xB_tmp[18] + ((2.0*(pow(k[0],2))*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[16])*xB_tmp[7])/(pow(k[0],2)) + ((k[0]*x_tmp[11] - 2.0*k[0]*x_tmp[12])*xB_tmp[11])/k[0] - (1.0*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12])*xB_tmp[13])/k[0] - (1.0*xB_tmp[6]*(k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[8] + 4.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[15]*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12] - 2.0*(pow(k[0],2))*x_tmp[15] + 4.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[4])/(pow(k[0],2)) + (xB_tmp[16]*(2.0*k[0]*x_tmp[12] - 2.0*(pow(k[0],2))*x_tmp[2] + 3.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[5]*(k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[8]*(2.0*k[0]*x_tmp[12] - 1.0*k[0]*x_tmp[11] - 1.0*(pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[12] - 4.0*(pow(k[0],2))*x_tmp[2])*xB_tmp[2])/(pow(k[0],2)) + (((pow(k[0],2))*x_tmp[17] - 2.0*(pow(k[0],2))*x_tmp[18])*xB_tmp[17])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[17] + 2.0*(pow(k[0],2))*x_tmp[18])*xB_tmp[19])/(pow(k[0],2)) + (xB_tmp[9]*(2.0*k[0]*x_tmp[12] - 2.0*(pow(k[0],2))*x_tmp[2] + 2.0*(pow(k[0],2))*x_tmp[9] - 1.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = 2.0*x_tmp[12]*xB_tmp[32] - 1.0*x_tmp[11]*xB_tmp[30] - 1.0*x_tmp[17]*xB_tmp[34] + 2.0*x_tmp[18]*xB_tmp[38] + ((2.0*(pow(k[0],2))*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[16])*xB_tmp[27])/(pow(k[0],2)) + ((k[0]*x_tmp[11] - 2.0*k[0]*x_tmp[12])*xB_tmp[31])/k[0] - (1.0*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12])*xB_tmp[33])/k[0] - (1.0*xB_tmp[26]*(k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[8] + 4.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[35]*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12] - 2.0*(pow(k[0],2))*x_tmp[15] + 4.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[24])/(pow(k[0],2)) + (xB_tmp[36]*(2.0*k[0]*x_tmp[12] - 2.0*(pow(k[0],2))*x_tmp[2] + 3.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[25]*(k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[28]*(2.0*k[0]*x_tmp[12] - 1.0*k[0]*x_tmp[11] - 1.0*(pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[12] - 4.0*(pow(k[0],2))*x_tmp[2])*xB_tmp[22])/(pow(k[0],2)) + (((pow(k[0],2))*x_tmp[17] - 2.0*(pow(k[0],2))*x_tmp[18])*xB_tmp[37])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[17] + 2.0*(pow(k[0],2))*x_tmp[18])*xB_tmp[39])/(pow(k[0],2)) + (xB_tmp[29]*(2.0*k[0]*x_tmp[12] - 2.0*(pow(k[0],2))*x_tmp[2] + 2.0*(pow(k[0],2))*x_tmp[9] - 1.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2));

  } break;

  case 2: {
qBdot_tmp[0+ip*ny] = - (1.0*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12])*xB_tmp[1])/k[0] - (1.0*xB_tmp[0]*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[17] + 4.0*(pow(k[0],2))*x_tmp[18]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[7])*xB_tmp[14])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[9])*xB_tmp[19])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[16])*xB_tmp[18])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[16])*xB_tmp[17])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = - (1.0*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12])*xB_tmp[21])/k[0] - (1.0*xB_tmp[20]*(k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[17] + 4.0*(pow(k[0],2))*x_tmp[18]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[7])*xB_tmp[34])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[9])*xB_tmp[39])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[16])*xB_tmp[38])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[16])*xB_tmp[37])/(pow(k[0],2));

  } break;

  case 3: {
qBdot_tmp[0+ip*ny] = x_tmp[1]*xB_tmp[1] + x_tmp[14]*xB_tmp[14] + x_tmp[17]*xB_tmp[17] + x_tmp[18]*xB_tmp[18] + x_tmp[19]*xB_tmp[19] - (1.0*(k[0]*x_tmp[1] - 2.0*(pow(k[0],2))*x_tmp[0])*xB_tmp[0])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = x_tmp[1]*xB_tmp[21] + x_tmp[14]*xB_tmp[34] + x_tmp[17]*xB_tmp[37] + x_tmp[18]*xB_tmp[38] + x_tmp[19]*xB_tmp[39] - (1.0*(k[0]*x_tmp[1] - 2.0*(pow(k[0],2))*x_tmp[0])*xB_tmp[20])/(pow(k[0],2));

  } break;

  case 4: {
qBdot_tmp[0+ip*ny] = - (16.0*xB_tmp[3])/(pow(k[0],2)) - (4.0*xB_tmp[13])/k[0];
qBdot_tmp[1+ip*ny] = - (16.0*xB_tmp[23])/(pow(k[0],2)) - (4.0*xB_tmp[33])/k[0];

  } break;

  case 5: {
qBdot_tmp[0+ip*ny] = x_tmp[6]*xB_tmp[6] + x_tmp[8]*xB_tmp[8] + x_tmp[9]*xB_tmp[9] + x_tmp[13]*xB_tmp[13] + x_tmp[19]*xB_tmp[19] - (1.0*(k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[3])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = x_tmp[6]*xB_tmp[26] + x_tmp[8]*xB_tmp[28] + x_tmp[9]*xB_tmp[29] + x_tmp[13]*xB_tmp[33] + x_tmp[19]*xB_tmp[39] - (1.0*(k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[23])/(pow(k[0],2));

  } break;

  }
  }

  for (iyp=0; iyp<2*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_enhancer_321(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*20);
x0_tmp[0] = (k[20]*k[40])/(pow(k[0],2));
x0_tmp[1] = (k[5]*k[25])/k[0];
x0_tmp[2] = (k[15]*k[35])/(pow(k[0],2));
x0_tmp[3] = (k[18]*k[38])/(pow(k[0],2));
x0_tmp[4] = (k[6]*k[26])/(pow(k[0],2));
x0_tmp[5] = (k[7]*k[27])/(pow(k[0],2));
x0_tmp[6] = (k[9]*k[29])/(pow(k[0],2));
x0_tmp[7] = (k[8]*k[28])/(pow(k[0],2));
x0_tmp[8] = (k[13]*k[33])/(pow(k[0],2));
x0_tmp[9] = (k[16]*k[36])/(pow(k[0],2));
x0_tmp[10] = (k[1]*k[21] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[11] = (k[2]*k[22])/k[0];
x0_tmp[12] = (k[3]*k[23])/k[0];
x0_tmp[13] = (k[4]*k[24])/k[0];
x0_tmp[14] = (k[10]*k[30])/(pow(k[0],2));
x0_tmp[15] = (k[11]*k[31])/(pow(k[0],2));
x0_tmp[16] = (k[12]*k[32])/(pow(k[0],2));
x0_tmp[17] = (k[14]*k[34])/(pow(k[0],2));
x0_tmp[18] = (k[17]*k[37])/(pow(k[0],2));
x0_tmp[19] = (k[19]*k[39])/(pow(k[0],2));
  
  
  return;
}


 int Jv_enhancer_321(N_Vector v, N_Vector Jv, realtype t,
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
  memset(Jv_tmp,0,sizeof(double)*20);
Jv_tmp[0] = 2.0*p[2]*v_tmp[17] - 2.0*p[3]*v_tmp[0] + 4.0*p[2]*v_tmp[18] + (p[3]*v_tmp[1])/k[0] + (p[2]*v_tmp[11])/k[0] + (2.0*p[2]*v_tmp[12])/k[0];
Jv_tmp[1] = p[2]*v_tmp[11] - 1.0*p[3]*v_tmp[1] + 2.0*p[2]*v_tmp[12];
Jv_tmp[2] = p[0]*v_tmp[8] - 4.0*p[1]*v_tmp[2] + (v_tmp[11]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (v_tmp[13]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*v_tmp[12])/(pow(k[0],2)) + 2.0*k[0]*p[0]*v_tmp[9]*x_tmp[11] + 2.0*k[0]*p[0]*v_tmp[16]*x_tmp[13];
Jv_tmp[3] = 4.0*p[1]*v_tmp[9] + (v_tmp[8]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (2.0*p[1]*v_tmp[12])/k[0] - (1.0*v_tmp[3]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (v_tmp[13]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (v_tmp[11]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*v_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[6])/(pow(k[0],2));
Jv_tmp[4] = 2.0*p[1]*v_tmp[5] + (p[1]*v_tmp[11])/k[0] - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*v_tmp[13])/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*v_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*v_tmp[6])/(pow(k[0],2)) - 4.0*k[0]*p[0]*v_tmp[4]*x_tmp[13];
Jv_tmp[5] = 2.0*p[1]*v_tmp[7] + p[1]*v_tmp[15] + (v_tmp[13]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[5])/(pow(k[0],2)) - (1.0*v_tmp[11]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*v_tmp[6]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (v_tmp[10]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*v_tmp[4]*x_tmp[13] - 2.0*k[0]*p[0]*v_tmp[8]*x_tmp[10];
Jv_tmp[6] = 2.0*p[1]*v_tmp[7] + p[1]*v_tmp[8] - (1.0*v_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[5])/(pow(k[0],2)) + (v_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (v_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*v_tmp[6]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*v_tmp[3]*x_tmp[10] - 2.0*k[0]*p[0]*v_tmp[4]*x_tmp[13];
Jv_tmp[7] = p[1]*v_tmp[16] - (1.0*v_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*v_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[7])/(pow(k[0],2)) - (1.0*v_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (2.0*p[0]*v_tmp[12]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*v_tmp[6]*x_tmp[11] + k[0]*p[0]*v_tmp[5]*x_tmp[13] - 2.0*k[0]*p[0]*v_tmp[9]*x_tmp[10];
Jv_tmp[8] = 2.0*p[1]*v_tmp[9] + 2.0*p[1]*v_tmp[16] + (2.0*p[1]*v_tmp[12])/k[0] - (1.0*v_tmp[11]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[15])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*v_tmp[3])/(pow(k[0],2)) - (1.0*v_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[6])/(pow(k[0],2)) + (v_tmp[10]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*v_tmp[8]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*v_tmp[5]*x_tmp[13];
Jv_tmp[9] = 2.0*p[1]*v_tmp[2] + (v_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*v_tmp[9]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (v_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[8])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[16])/(pow(k[0],2)) + (v_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*v_tmp[10])/(pow(k[0],2)) + k[0]*p[0]*v_tmp[3]*x_tmp[11] - 2.0*k[0]*p[0]*v_tmp[7]*x_tmp[13];
Jv_tmp[10] = p[1]*v_tmp[11] - 2.0*k[0]*p[0]*v_tmp[6] - 2.0*k[0]*p[0]*v_tmp[10]*x_tmp[13] - 2.0*k[0]*p[0]*v_tmp[13]*x_tmp[10];
Jv_tmp[11] = 2.0*p[1]*v_tmp[12] + ((2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*v_tmp[13])/k[0] + 2.0*k[0]*p[0]*v_tmp[6] - 1.0*k[0]*p[0]*v_tmp[8] - (1.0*(k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*v_tmp[11])/k[0] + 2.0*k[0]*p[0]*v_tmp[10]*x_tmp[13];
Jv_tmp[12] = k[0]*p[0]*v_tmp[8] - 2.0*p[1]*v_tmp[12] + k[0]*p[0]*v_tmp[11]*x_tmp[13] + k[0]*p[0]*v_tmp[13]*x_tmp[11];
Jv_tmp[13] = 2.0*p[1]*v_tmp[12] + ((k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*v_tmp[11])/k[0] - (1.0*v_tmp[13]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - 2.0*k[0]*p[0]*v_tmp[6] - 1.0*k[0]*p[0]*v_tmp[8] - 2.0*k[0]*p[0]*v_tmp[10]*x_tmp[13];
Jv_tmp[14] = p[2]*v_tmp[5] + 2.0*p[2]*v_tmp[7] + p[1]*v_tmp[17] - (1.0*((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*v_tmp[13])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*v_tmp[10])/(pow(k[0],2)) + (2.0*p[0]*v_tmp[1]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*v_tmp[19]*x_tmp[10];
Jv_tmp[15] = 2.0*p[0]*v_tmp[6] + 4.0*p[1]*v_tmp[16] + (v_tmp[13]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + (2.0*p[1]*v_tmp[12])/k[0] + (v_tmp[10]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*v_tmp[11]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[15])/(pow(k[0],2)) + (v_tmp[8]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 4.0*k[0]*p[0]*v_tmp[5]*x_tmp[13];
Jv_tmp[16] = 2.0*p[1]*v_tmp[2] + (v_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*v_tmp[8])/(pow(k[0],2)) - (1.0*v_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*v_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*v_tmp[9])/(pow(k[0],2)) - (1.0*(3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*v_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*v_tmp[7]*x_tmp[13] + k[0]*p[0]*v_tmp[15]*x_tmp[13];
Jv_tmp[17] = p[2]*v_tmp[15] + 2.0*p[2]*v_tmp[16] + 2.0*p[1]*v_tmp[18] - (1.0*v_tmp[1]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*v_tmp[19])/(pow(k[0],2)) - (1.0*v_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*v_tmp[11])/(pow(k[0],2)) - (1.0*v_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*v_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*v_tmp[14]*x_tmp[13];
Jv_tmp[18] = 2.0*p[2]*v_tmp[2] + p[2]*v_tmp[16] - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*v_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*v_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*v_tmp[11])/(pow(k[0],2)) - (1.0*p[0]*v_tmp[1]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*v_tmp[17]*x_tmp[13] + k[0]*p[0]*v_tmp[19]*x_tmp[11];
Jv_tmp[19] = p[2]*v_tmp[8] + 2.0*p[2]*v_tmp[9] + 2.0*p[1]*v_tmp[18] - (1.0*v_tmp[19]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*v_tmp[1])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*v_tmp[17])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*v_tmp[11])/(pow(k[0],2)) + (v_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*v_tmp[10])/(pow(k[0],2)) - 2.0*k[0]*p[0]*v_tmp[14]*x_tmp[13];

  for (ix=0; ix<20; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_enhancer_321(N_Vector vB, N_Vector JvB, realtype t,
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
  memset(JvB_tmp,0,sizeof(double)*20);
JvB_tmp[0] = 2.0*p[3]*vB_tmp[0];
JvB_tmp[1] = p[3]*vB_tmp[1] - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*vB_tmp[19])/(pow(k[0],2)) - (1.0*p[3]*vB_tmp[0])/k[0] + (vB_tmp[17]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) - (2.0*p[0]*vB_tmp[14]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + (p[0]*vB_tmp[18]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
JvB_tmp[2] = 4.0*p[1]*vB_tmp[2] - 2.0*p[1]*vB_tmp[9] - 2.0*p[1]*vB_tmp[16] - 2.0*p[2]*vB_tmp[18];
JvB_tmp[3] = (vB_tmp[3]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*vB_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[0]*vB_tmp[6]*x_tmp[10] - 1.0*k[0]*p[0]*vB_tmp[9]*x_tmp[11];
JvB_tmp[4] = 4.0*k[0]*p[0]*vB_tmp[4]*x_tmp[13] - 2.0*k[0]*p[0]*vB_tmp[5]*x_tmp[13] + 2.0*k[0]*p[0]*vB_tmp[6]*x_tmp[13];
JvB_tmp[5] = (((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[5])/(pow(k[0],2)) - 1.0*p[2]*vB_tmp[14] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[6])/(pow(k[0],2)) - 2.0*p[1]*vB_tmp[4] - 1.0*k[0]*p[0]*vB_tmp[7]*x_tmp[13] + 2.0*k[0]*p[0]*vB_tmp[8]*x_tmp[13] - 4.0*k[0]*p[0]*vB_tmp[15]*x_tmp[13];
JvB_tmp[6] = (vB_tmp[6]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*p[0]*vB_tmp[15] - (1.0*(2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*vB_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[3])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[8])/(pow(k[0],2)) + 2.0*k[0]*p[0]*vB_tmp[10] - 2.0*k[0]*p[0]*vB_tmp[11] + 2.0*k[0]*p[0]*vB_tmp[13] + (vB_tmp[5]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*vB_tmp[7]*x_tmp[11];
JvB_tmp[7] = ((2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[7])/(pow(k[0],2)) - 2.0*p[1]*vB_tmp[6] - 2.0*p[2]*vB_tmp[14] - 2.0*p[1]*vB_tmp[5] + 2.0*k[0]*p[0]*vB_tmp[9]*x_tmp[13] - 2.0*k[0]*p[0]*vB_tmp[16]*x_tmp[13];
JvB_tmp[8] = (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[9])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[6] - 1.0*p[2]*vB_tmp[19] - (1.0*vB_tmp[3]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[2] + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*vB_tmp[16])/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[11] - 1.0*k[0]*p[0]*vB_tmp[12] + k[0]*p[0]*vB_tmp[13] + (vB_tmp[8]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*vB_tmp[15]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*vB_tmp[5]*x_tmp[10];
JvB_tmp[9] = (vB_tmp[9]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - 2.0*p[1]*vB_tmp[8] - 2.0*p[2]*vB_tmp[19] - 4.0*p[1]*vB_tmp[3] - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*vB_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*vB_tmp[2]*x_tmp[11] + 2.0*k[0]*p[0]*vB_tmp[7]*x_tmp[10];
JvB_tmp[10] = (vB_tmp[7]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*vB_tmp[15]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*vB_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*vB_tmp[3])/(pow(k[0],2)) + ((4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*vB_tmp[4])/(pow(k[0],2)) - (1.0*vB_tmp[8]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*vB_tmp[9])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*vB_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*vB_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*vB_tmp[17])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*vB_tmp[19])/(pow(k[0],2)) - (1.0*vB_tmp[5]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*vB_tmp[10]*x_tmp[13] - 2.0*k[0]*p[0]*vB_tmp[11]*x_tmp[13] + 2.0*k[0]*p[0]*vB_tmp[13]*x_tmp[13];
JvB_tmp[11] = (vB_tmp[8]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[10] - (1.0*vB_tmp[2]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*vB_tmp[13])/k[0] - (1.0*p[2]*vB_tmp[0])/k[0] - (1.0*p[1]*vB_tmp[4])/k[0] - 1.0*p[2]*vB_tmp[1] - (1.0*vB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*vB_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*vB_tmp[3]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (vB_tmp[15]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*vB_tmp[6]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*vB_tmp[7])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*vB_tmp[17])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*vB_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*vB_tmp[19])/(pow(k[0],2)) + (vB_tmp[5]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) + ((k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*vB_tmp[11])/k[0] - 1.0*k[0]*p[0]*vB_tmp[12]*x_tmp[13];
JvB_tmp[12] = 2.0*p[1]*vB_tmp[12] - 2.0*p[1]*vB_tmp[11] - 2.0*p[2]*vB_tmp[1] - 2.0*p[1]*vB_tmp[13] - (2.0*p[2]*vB_tmp[0])/k[0] - (2.0*p[1]*vB_tmp[3])/k[0] - (2.0*p[1]*vB_tmp[8])/k[0] - (2.0*p[1]*vB_tmp[15])/k[0] - (1.0*vB_tmp[9]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) + (vB_tmp[16]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*vB_tmp[2])/(pow(k[0],2)) - (2.0*p[0]*vB_tmp[7]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
JvB_tmp[13] = (vB_tmp[6]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*vB_tmp[15]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*vB_tmp[5]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*vB_tmp[2]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*vB_tmp[3]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (vB_tmp[13]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - (1.0*vB_tmp[9]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (vB_tmp[16]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*vB_tmp[4])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*vB_tmp[11])/k[0] + (vB_tmp[8]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*vB_tmp[18])/(pow(k[0],2)) + (vB_tmp[7]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (vB_tmp[17]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*vB_tmp[19]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*vB_tmp[14])/(pow(k[0],2)) + 2.0*k[0]*p[0]*vB_tmp[10]*x_tmp[10] - 1.0*k[0]*p[0]*vB_tmp[12]*x_tmp[11];
JvB_tmp[14] = (((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[14])/(pow(k[0],2)) - 2.0*k[0]*p[0]*vB_tmp[17]*x_tmp[13] + 2.0*k[0]*p[0]*vB_tmp[19]*x_tmp[13];
JvB_tmp[15] = ((2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[15])/(pow(k[0],2)) - 1.0*p[2]*vB_tmp[17] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[8])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[5] - 1.0*k[0]*p[0]*vB_tmp[16]*x_tmp[13];
JvB_tmp[16] = ((3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[16])/(pow(k[0],2)) - 2.0*p[1]*vB_tmp[8] - 4.0*p[1]*vB_tmp[15] - 2.0*p[2]*vB_tmp[17] - 1.0*p[2]*vB_tmp[18] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[9])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[7] - 2.0*k[0]*p[0]*vB_tmp[2]*x_tmp[13];
JvB_tmp[17] = (vB_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[14] - (1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*vB_tmp[19])/(pow(k[0],2)) - 2.0*p[2]*vB_tmp[0] - 1.0*k[0]*p[0]*vB_tmp[18]*x_tmp[13];
JvB_tmp[18] = ((2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*vB_tmp[18])/(pow(k[0],2)) - 2.0*p[1]*vB_tmp[17] - 2.0*p[1]*vB_tmp[19] - 4.0*p[2]*vB_tmp[0];
JvB_tmp[19] = (vB_tmp[19]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*vB_tmp[17])/(pow(k[0],2)) + 2.0*k[0]*p[0]*vB_tmp[14]*x_tmp[10] - 1.0*k[0]*p[0]*vB_tmp[18]*x_tmp[11];

  for (ix=0; ix<20; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_enhancer_321(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_enhancer_321(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_enhancer_321(long int N, realtype t, N_Vector x,
  	N_Vector xdot, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  memset(J->data,0,sizeof(double)*400);
J->data[0] = -2.0*p[3];
J->data[20] = p[3]/k[0];
J->data[21] = -1.0*p[3];
J->data[34] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
J->data[37] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2));
J->data[38] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
J->data[39] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
J->data[42] = -4.0*p[1];
J->data[49] = 2.0*p[1];
J->data[56] = 2.0*p[1];
J->data[58] = 2.0*p[2];
J->data[63] = -(1.0*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[66] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[68] = (2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[69] = k[0]*p[0]*x_tmp[11];
J->data[84] = -4.0*k[0]*p[0]*x_tmp[13];
J->data[85] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[86] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[104] = 2.0*p[1];
J->data[105] = -(1.0*((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[106] = ((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[107] = k[0]*p[0]*x_tmp[13];
J->data[108] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[114] = p[2];
J->data[115] = 4.0*k[0]*p[0]*x_tmp[13];
J->data[123] = (2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[124] = (2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])/(pow(k[0],2));
J->data[125] = -(1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[126] = -(1.0*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[127] = k[0]*p[0]*x_tmp[11];
J->data[128] = -(1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[130] = -2.0*k[0]*p[0];
J->data[131] = 2.0*k[0]*p[0];
J->data[133] = -2.0*k[0]*p[0];
J->data[135] = 2.0*p[0];
J->data[145] = 2.0*p[1];
J->data[146] = 2.0*p[1];
J->data[147] = -(1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[149] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[154] = 2.0*p[2];
J->data[156] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[162] = p[0];
J->data[163] = ((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[165] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[166] = p[1];
J->data[168] = -(1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[169] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[171] = -1.0*k[0]*p[0];
J->data[172] = k[0]*p[0];
J->data[173] = -1.0*k[0]*p[0];
J->data[175] = ((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[176] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[179] = p[2];
J->data[182] = 2.0*k[0]*p[0]*x_tmp[11];
J->data[183] = 4.0*p[1];
J->data[187] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[188] = 2.0*p[1];
J->data[189] = -(1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[196] = (2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[199] = 2.0*p[2];
J->data[203] = -(1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[204] = -(1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[205] = (3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[206] = (p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[207] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2));
J->data[208] = (2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[209] = (2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
J->data[210] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[211] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[213] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[214] = (2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
J->data[215] = (4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[216] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2));
J->data[217] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19]))/(pow(k[0],2));
J->data[219] = (2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
J->data[220] = p[2]/k[0];
J->data[221] = p[2];
J->data[222] = (2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[223] = (k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[224] = p[1]/k[0];
J->data[225] = -(1.0*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[226] = (p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
J->data[227] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[228] = -(1.0*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[229] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[230] = p[1];
J->data[231] = -(1.0*(k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13]))/k[0];
J->data[232] = k[0]*p[0]*x_tmp[13];
J->data[233] = (k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])/k[0];
J->data[235] = -(1.0*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[236] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[237] = (p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
J->data[238] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19]))/(pow(k[0],2));
J->data[239] = (p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
J->data[240] = (2.0*p[2])/k[0];
J->data[241] = 2.0*p[2];
J->data[242] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2));
J->data[243] = (2.0*p[1])/k[0];
J->data[247] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
J->data[248] = (2.0*p[1])/k[0];
J->data[249] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])/(pow(k[0],2));
J->data[251] = 2.0*p[1];
J->data[252] = -2.0*p[1];
J->data[253] = 2.0*p[1];
J->data[255] = (2.0*p[1])/k[0];
J->data[256] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2));
J->data[262] = ((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
J->data[263] = (k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[264] = -(1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2));
J->data[265] = (3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])/(pow(k[0],2));
J->data[266] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2));
J->data[267] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2));
J->data[268] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2));
J->data[269] = (2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
J->data[270] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[271] = (2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])/k[0];
J->data[272] = k[0]*p[0]*x_tmp[11];
J->data[273] = -(1.0*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0];
J->data[274] = (2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])/(pow(k[0],2));
J->data[275] = (4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15])/(pow(k[0],2));
J->data[276] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[277] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2));
J->data[278] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2));
J->data[279] = (2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])/(pow(k[0],2));
J->data[294] = -(1.0*((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[297] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[299] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[305] = p[1];
J->data[308] = ((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[315] = -(1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[316] = k[0]*p[0]*x_tmp[13];
J->data[317] = p[2];
J->data[322] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[327] = p[1];
J->data[328] = 2.0*p[1];
J->data[329] = ((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[335] = 4.0*p[1];
J->data[336] = -(1.0*(3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[337] = 2.0*p[2];
J->data[338] = p[2];
J->data[340] = 2.0*p[2];
J->data[354] = p[1];
J->data[357] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[358] = k[0]*p[0]*x_tmp[13];
J->data[359] = ((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[360] = 4.0*p[2];
J->data[377] = 2.0*p[1];
J->data[378] = -(1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[379] = 2.0*p[1];
J->data[394] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[397] = (2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[398] = k[0]*p[0]*x_tmp[11];
J->data[399] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));

  for (iJ=0; iJ<400; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_enhancer_321(realtype t, N_Vector x,
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
  J->rowvals[3] = 14;
  J->rowvals[4] = 17;
  J->rowvals[5] = 18;
  J->rowvals[6] = 19;
  J->rowvals[7] = 2;
  J->rowvals[8] = 9;
  J->rowvals[9] = 16;
  J->rowvals[10] = 18;
  J->rowvals[11] = 3;
  J->rowvals[12] = 6;
  J->rowvals[13] = 8;
  J->rowvals[14] = 9;
  J->rowvals[15] = 4;
  J->rowvals[16] = 5;
  J->rowvals[17] = 6;
  J->rowvals[18] = 4;
  J->rowvals[19] = 5;
  J->rowvals[20] = 6;
  J->rowvals[21] = 7;
  J->rowvals[22] = 8;
  J->rowvals[23] = 14;
  J->rowvals[24] = 15;
  J->rowvals[25] = 3;
  J->rowvals[26] = 4;
  J->rowvals[27] = 5;
  J->rowvals[28] = 6;
  J->rowvals[29] = 7;
  J->rowvals[30] = 8;
  J->rowvals[31] = 10;
  J->rowvals[32] = 11;
  J->rowvals[33] = 13;
  J->rowvals[34] = 15;
  J->rowvals[35] = 5;
  J->rowvals[36] = 6;
  J->rowvals[37] = 7;
  J->rowvals[38] = 9;
  J->rowvals[39] = 14;
  J->rowvals[40] = 16;
  J->rowvals[41] = 2;
  J->rowvals[42] = 3;
  J->rowvals[43] = 5;
  J->rowvals[44] = 6;
  J->rowvals[45] = 8;
  J->rowvals[46] = 9;
  J->rowvals[47] = 11;
  J->rowvals[48] = 12;
  J->rowvals[49] = 13;
  J->rowvals[50] = 15;
  J->rowvals[51] = 16;
  J->rowvals[52] = 19;
  J->rowvals[53] = 2;
  J->rowvals[54] = 3;
  J->rowvals[55] = 7;
  J->rowvals[56] = 8;
  J->rowvals[57] = 9;
  J->rowvals[58] = 16;
  J->rowvals[59] = 19;
  J->rowvals[60] = 3;
  J->rowvals[61] = 4;
  J->rowvals[62] = 5;
  J->rowvals[63] = 6;
  J->rowvals[64] = 7;
  J->rowvals[65] = 8;
  J->rowvals[66] = 9;
  J->rowvals[67] = 10;
  J->rowvals[68] = 11;
  J->rowvals[69] = 13;
  J->rowvals[70] = 14;
  J->rowvals[71] = 15;
  J->rowvals[72] = 16;
  J->rowvals[73] = 17;
  J->rowvals[74] = 19;
  J->rowvals[75] = 0;
  J->rowvals[76] = 1;
  J->rowvals[77] = 2;
  J->rowvals[78] = 3;
  J->rowvals[79] = 4;
  J->rowvals[80] = 5;
  J->rowvals[81] = 6;
  J->rowvals[82] = 7;
  J->rowvals[83] = 8;
  J->rowvals[84] = 9;
  J->rowvals[85] = 10;
  J->rowvals[86] = 11;
  J->rowvals[87] = 12;
  J->rowvals[88] = 13;
  J->rowvals[89] = 15;
  J->rowvals[90] = 16;
  J->rowvals[91] = 17;
  J->rowvals[92] = 18;
  J->rowvals[93] = 19;
  J->rowvals[94] = 0;
  J->rowvals[95] = 1;
  J->rowvals[96] = 2;
  J->rowvals[97] = 3;
  J->rowvals[98] = 7;
  J->rowvals[99] = 8;
  J->rowvals[100] = 9;
  J->rowvals[101] = 11;
  J->rowvals[102] = 12;
  J->rowvals[103] = 13;
  J->rowvals[104] = 15;
  J->rowvals[105] = 16;
  J->rowvals[106] = 2;
  J->rowvals[107] = 3;
  J->rowvals[108] = 4;
  J->rowvals[109] = 5;
  J->rowvals[110] = 6;
  J->rowvals[111] = 7;
  J->rowvals[112] = 8;
  J->rowvals[113] = 9;
  J->rowvals[114] = 10;
  J->rowvals[115] = 11;
  J->rowvals[116] = 12;
  J->rowvals[117] = 13;
  J->rowvals[118] = 14;
  J->rowvals[119] = 15;
  J->rowvals[120] = 16;
  J->rowvals[121] = 17;
  J->rowvals[122] = 18;
  J->rowvals[123] = 19;
  J->rowvals[124] = 14;
  J->rowvals[125] = 17;
  J->rowvals[126] = 19;
  J->rowvals[127] = 5;
  J->rowvals[128] = 8;
  J->rowvals[129] = 15;
  J->rowvals[130] = 16;
  J->rowvals[131] = 17;
  J->rowvals[132] = 2;
  J->rowvals[133] = 7;
  J->rowvals[134] = 8;
  J->rowvals[135] = 9;
  J->rowvals[136] = 15;
  J->rowvals[137] = 16;
  J->rowvals[138] = 17;
  J->rowvals[139] = 18;
  J->rowvals[140] = 0;
  J->rowvals[141] = 14;
  J->rowvals[142] = 17;
  J->rowvals[143] = 18;
  J->rowvals[144] = 19;
  J->rowvals[145] = 0;
  J->rowvals[146] = 17;
  J->rowvals[147] = 18;
  J->rowvals[148] = 19;
  J->rowvals[149] = 14;
  J->rowvals[150] = 17;
  J->rowvals[151] = 18;
  J->rowvals[152] = 19;
  J->colptrs[0] = 0;
  J->colptrs[1] = 1;
  J->colptrs[2] = 7;
  J->colptrs[3] = 11;
  J->colptrs[4] = 15;
  J->colptrs[5] = 18;
  J->colptrs[6] = 25;
  J->colptrs[7] = 35;
  J->colptrs[8] = 41;
  J->colptrs[9] = 53;
  J->colptrs[10] = 60;
  J->colptrs[11] = 75;
  J->colptrs[12] = 94;
  J->colptrs[13] = 106;
  J->colptrs[14] = 124;
  J->colptrs[15] = 127;
  J->colptrs[16] = 132;
  J->colptrs[17] = 140;
  J->colptrs[18] = 145;
  J->colptrs[19] = 149;
  J->colptrs[20] = 153;
J->data[0] = -2.0*p[3];
J->data[1] = p[3]/k[0];
J->data[2] = -1.0*p[3];
J->data[3] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
J->data[4] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2));
J->data[5] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
J->data[6] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
J->data[7] = -4.0*p[1];
J->data[8] = 2.0*p[1];
J->data[9] = 2.0*p[1];
J->data[10] = 2.0*p[2];
J->data[11] = -(1.0*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[12] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[13] = (2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[14] = k[0]*p[0]*x_tmp[11];
J->data[15] = -4.0*k[0]*p[0]*x_tmp[13];
J->data[16] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[17] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[18] = 2.0*p[1];
J->data[19] = -(1.0*((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[20] = ((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[21] = k[0]*p[0]*x_tmp[13];
J->data[22] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[23] = p[2];
J->data[24] = 4.0*k[0]*p[0]*x_tmp[13];
J->data[25] = (2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[26] = (2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])/(pow(k[0],2));
J->data[27] = -(1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[28] = -(1.0*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[29] = k[0]*p[0]*x_tmp[11];
J->data[30] = -(1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[31] = -2.0*k[0]*p[0];
J->data[32] = 2.0*k[0]*p[0];
J->data[33] = -2.0*k[0]*p[0];
J->data[34] = 2.0*p[0];
J->data[35] = 2.0*p[1];
J->data[36] = 2.0*p[1];
J->data[37] = -(1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[38] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[39] = 2.0*p[2];
J->data[40] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[41] = p[0];
J->data[42] = ((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[43] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[44] = p[1];
J->data[45] = -(1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[46] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[47] = -1.0*k[0]*p[0];
J->data[48] = k[0]*p[0];
J->data[49] = -1.0*k[0]*p[0];
J->data[50] = ((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[51] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[52] = p[2];
J->data[53] = 2.0*k[0]*p[0]*x_tmp[11];
J->data[54] = 4.0*p[1];
J->data[55] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[56] = 2.0*p[1];
J->data[57] = -(1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[58] = (2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[59] = 2.0*p[2];
J->data[60] = -(1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[61] = -(1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[62] = (3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[63] = (p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[64] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2));
J->data[65] = (2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[66] = (2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
J->data[67] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[68] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[69] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[70] = (2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
J->data[71] = (4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[72] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2));
J->data[73] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19]))/(pow(k[0],2));
J->data[74] = (2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
J->data[75] = p[2]/k[0];
J->data[76] = p[2];
J->data[77] = (2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[78] = (k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[79] = p[1]/k[0];
J->data[80] = -(1.0*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[81] = (p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
J->data[82] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
J->data[83] = -(1.0*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[84] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[85] = p[1];
J->data[86] = -(1.0*(k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13]))/k[0];
J->data[87] = k[0]*p[0]*x_tmp[13];
J->data[88] = (k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])/k[0];
J->data[89] = -(1.0*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[90] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[91] = (p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
J->data[92] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19]))/(pow(k[0],2));
J->data[93] = (p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
J->data[94] = (2.0*p[2])/k[0];
J->data[95] = 2.0*p[2];
J->data[96] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2));
J->data[97] = (2.0*p[1])/k[0];
J->data[98] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
J->data[99] = (2.0*p[1])/k[0];
J->data[100] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])/(pow(k[0],2));
J->data[101] = 2.0*p[1];
J->data[102] = -2.0*p[1];
J->data[103] = 2.0*p[1];
J->data[104] = (2.0*p[1])/k[0];
J->data[105] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2));
J->data[106] = ((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
J->data[107] = (k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[108] = -(1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2));
J->data[109] = (3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])/(pow(k[0],2));
J->data[110] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2));
J->data[111] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2));
J->data[112] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2));
J->data[113] = (2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
J->data[114] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[115] = (2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])/k[0];
J->data[116] = k[0]*p[0]*x_tmp[11];
J->data[117] = -(1.0*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0];
J->data[118] = (2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])/(pow(k[0],2));
J->data[119] = (4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15])/(pow(k[0],2));
J->data[120] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[121] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2));
J->data[122] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2));
J->data[123] = (2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])/(pow(k[0],2));
J->data[124] = -(1.0*((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[125] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[126] = -2.0*k[0]*p[0]*x_tmp[13];
J->data[127] = p[1];
J->data[128] = ((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[129] = -(1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[130] = k[0]*p[0]*x_tmp[13];
J->data[131] = p[2];
J->data[132] = 2.0*k[0]*p[0]*x_tmp[13];
J->data[133] = p[1];
J->data[134] = 2.0*p[1];
J->data[135] = ((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[136] = 4.0*p[1];
J->data[137] = -(1.0*(3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[138] = 2.0*p[2];
J->data[139] = p[2];
J->data[140] = 2.0*p[2];
J->data[141] = p[1];
J->data[142] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
J->data[143] = k[0]*p[0]*x_tmp[13];
J->data[144] = ((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
J->data[145] = 4.0*p[2];
J->data[146] = 2.0*p[1];
J->data[147] = -(1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[148] = 2.0*p[1];
J->data[149] = -2.0*k[0]*p[0]*x_tmp[10];
J->data[150] = (2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[151] = k[0]*p[0]*x_tmp[11];
J->data[152] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
  return(0);
}


 int JBBand_enhancer_321(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_enhancer_321(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_enhancer_321(long int N, realtype t, N_Vector x,
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
JB->data[11] = -(1.0*p[2])/k[0];
JB->data[12] = -(2.0*p[2])/k[0];
JB->data[17] = -2.0*p[2];
JB->data[18] = -4.0*p[2];
JB->data[21] = p[3];
JB->data[31] = -1.0*p[2];
JB->data[32] = -2.0*p[2];
JB->data[42] = 4.0*p[1];
JB->data[48] = -1.0*p[0];
JB->data[49] = -2.0*k[0]*p[0]*x_tmp[11];
JB->data[51] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[52] = (2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])/(pow(k[0],2));
JB->data[53] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
JB->data[56] = -2.0*k[0]*p[0]*x_tmp[13];
JB->data[63] = (2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
JB->data[66] = -(1.0*(2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[68] = -(1.0*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[69] = -4.0*p[1];
JB->data[70] = (4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[71] = -(1.0*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[72] = -(2.0*p[1])/k[0];
JB->data[73] = -(1.0*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2));
JB->data[84] = 4.0*k[0]*p[0]*x_tmp[13];
JB->data[85] = -2.0*p[1];
JB->data[86] = -(1.0*(2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10]))/(pow(k[0],2));
JB->data[90] = (4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[91] = -(1.0*p[1])/k[0];
JB->data[93] = (4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])/(pow(k[0],2));
JB->data[104] = -2.0*k[0]*p[0]*x_tmp[13];
JB->data[105] = ((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[106] = (2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
JB->data[107] = -2.0*p[1];
JB->data[108] = 2.0*k[0]*p[0]*x_tmp[10];
JB->data[110] = -(1.0*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[111] = (k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
JB->data[113] = -(1.0*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2));
JB->data[115] = -1.0*p[1];
JB->data[123] = 2.0*k[0]*p[0]*x_tmp[10];
JB->data[124] = 2.0*k[0]*p[0]*x_tmp[13];
JB->data[125] = -(1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[126] = ((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[127] = -2.0*p[1];
JB->data[128] = -1.0*p[1];
JB->data[130] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[131] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2));
JB->data[133] = (2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])/(pow(k[0],2));
JB->data[145] = -1.0*k[0]*p[0]*x_tmp[13];
JB->data[146] = -1.0*k[0]*p[0]*x_tmp[11];
JB->data[147] = (2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[149] = 2.0*k[0]*p[0]*x_tmp[10];
JB->data[150] = (p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
JB->data[151] = (p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])/(pow(k[0],2));
JB->data[152] = -(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
JB->data[153] = (p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7])/(pow(k[0],2));
JB->data[156] = -1.0*p[1];
JB->data[163] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
JB->data[165] = 2.0*k[0]*p[0]*x_tmp[13];
JB->data[166] = (2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[168] = ((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[169] = -2.0*p[1];
JB->data[170] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[171] = (k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[172] = -(2.0*p[1])/k[0];
JB->data[173] = (2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15])/(pow(k[0],2));
JB->data[175] = -(1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[176] = -2.0*p[1];
JB->data[182] = -2.0*p[1];
JB->data[183] = -1.0*k[0]*p[0]*x_tmp[11];
JB->data[187] = 2.0*k[0]*p[0]*x_tmp[13];
JB->data[188] = ((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[189] = (2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
JB->data[190] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2));
JB->data[191] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[192] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2));
JB->data[193] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
JB->data[196] = -(1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[206] = 2.0*k[0]*p[0];
JB->data[210] = 2.0*k[0]*p[0]*x_tmp[13];
JB->data[211] = -1.0*p[1];
JB->data[213] = 2.0*k[0]*p[0]*x_tmp[10];
JB->data[226] = -2.0*k[0]*p[0];
JB->data[228] = k[0]*p[0];
JB->data[230] = -2.0*k[0]*p[0]*x_tmp[13];
JB->data[231] = (k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])/k[0];
JB->data[232] = -2.0*p[1];
JB->data[233] = -(1.0*(2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11]))/k[0];
JB->data[248] = -1.0*k[0]*p[0];
JB->data[251] = -1.0*k[0]*p[0]*x_tmp[13];
JB->data[252] = 2.0*p[1];
JB->data[253] = -1.0*k[0]*p[0]*x_tmp[11];
JB->data[266] = 2.0*k[0]*p[0];
JB->data[268] = k[0]*p[0];
JB->data[270] = 2.0*k[0]*p[0]*x_tmp[13];
JB->data[271] = -(1.0*(k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/k[0];
JB->data[272] = -2.0*p[1];
JB->data[273] = (k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11])/k[0];
JB->data[281] = -(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2));
JB->data[285] = -1.0*p[2];
JB->data[287] = -2.0*p[2];
JB->data[290] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19]))/(pow(k[0],2));
JB->data[293] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
JB->data[294] = ((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[297] = -1.0*p[1];
JB->data[299] = 2.0*k[0]*p[0]*x_tmp[10];
JB->data[305] = -4.0*k[0]*p[0]*x_tmp[13];
JB->data[306] = -2.0*p[0];
JB->data[308] = -(1.0*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
JB->data[310] = -(1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[311] = (4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[312] = -(2.0*p[1])/k[0];
JB->data[313] = -(1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2));
JB->data[315] = (2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[316] = -4.0*p[1];
JB->data[322] = -2.0*p[1];
JB->data[327] = -2.0*k[0]*p[0]*x_tmp[13];
JB->data[328] = ((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
JB->data[329] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
JB->data[330] = (2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])/(pow(k[0],2));
JB->data[331] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[332] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1])/(pow(k[0],2));
JB->data[333] = (2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[335] = -1.0*k[0]*p[0]*x_tmp[13];
JB->data[336] = (3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[341] = (2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
JB->data[350] = (2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
JB->data[351] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19]))/(pow(k[0],2));
JB->data[353] = (2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17])/(pow(k[0],2));
JB->data[354] = -2.0*k[0]*p[0]*x_tmp[13];
JB->data[355] = -1.0*p[2];
JB->data[356] = -2.0*p[2];
JB->data[357] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13])/(pow(k[0],2));
JB->data[358] = -2.0*p[1];
JB->data[359] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
JB->data[361] = (p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2));
JB->data[362] = -2.0*p[2];
JB->data[371] = (p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])/(pow(k[0],2));
JB->data[373] = (p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])/(pow(k[0],2));
JB->data[376] = -1.0*p[2];
JB->data[377] = -1.0*k[0]*p[0]*x_tmp[13];
JB->data[378] = (2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])/(pow(k[0],2));
JB->data[379] = -1.0*k[0]*p[0]*x_tmp[11];
JB->data[381] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2));
JB->data[388] = -1.0*p[2];
JB->data[389] = -2.0*p[2];
JB->data[390] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19]))/(pow(k[0],2));
JB->data[391] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19]))/(pow(k[0],2));
JB->data[393] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2));
JB->data[394] = 2.0*k[0]*p[0]*x_tmp[13];
JB->data[397] = -(1.0*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2));
JB->data[398] = -2.0*p[1];
JB->data[399] = ((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));

  for (iJ=0; iJ<400; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_enhancer_321(realtype t, N_Vector x,
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
  JB->rowvals[2] = 11;
  JB->rowvals[3] = 12;
  JB->rowvals[4] = 17;
  JB->rowvals[5] = 18;
  JB->rowvals[6] = 1;
  JB->rowvals[7] = 11;
  JB->rowvals[8] = 12;
  JB->rowvals[9] = 2;
  JB->rowvals[10] = 8;
  JB->rowvals[11] = 9;
  JB->rowvals[12] = 11;
  JB->rowvals[13] = 12;
  JB->rowvals[14] = 13;
  JB->rowvals[15] = 16;
  JB->rowvals[16] = 3;
  JB->rowvals[17] = 6;
  JB->rowvals[18] = 8;
  JB->rowvals[19] = 9;
  JB->rowvals[20] = 10;
  JB->rowvals[21] = 11;
  JB->rowvals[22] = 12;
  JB->rowvals[23] = 13;
  JB->rowvals[24] = 4;
  JB->rowvals[25] = 5;
  JB->rowvals[26] = 6;
  JB->rowvals[27] = 10;
  JB->rowvals[28] = 11;
  JB->rowvals[29] = 13;
  JB->rowvals[30] = 4;
  JB->rowvals[31] = 5;
  JB->rowvals[32] = 6;
  JB->rowvals[33] = 7;
  JB->rowvals[34] = 8;
  JB->rowvals[35] = 10;
  JB->rowvals[36] = 11;
  JB->rowvals[37] = 13;
  JB->rowvals[38] = 15;
  JB->rowvals[39] = 3;
  JB->rowvals[40] = 4;
  JB->rowvals[41] = 5;
  JB->rowvals[42] = 6;
  JB->rowvals[43] = 7;
  JB->rowvals[44] = 8;
  JB->rowvals[45] = 10;
  JB->rowvals[46] = 11;
  JB->rowvals[47] = 13;
  JB->rowvals[48] = 5;
  JB->rowvals[49] = 6;
  JB->rowvals[50] = 7;
  JB->rowvals[51] = 9;
  JB->rowvals[52] = 10;
  JB->rowvals[53] = 11;
  JB->rowvals[54] = 12;
  JB->rowvals[55] = 13;
  JB->rowvals[56] = 16;
  JB->rowvals[57] = 3;
  JB->rowvals[58] = 5;
  JB->rowvals[59] = 6;
  JB->rowvals[60] = 8;
  JB->rowvals[61] = 9;
  JB->rowvals[62] = 10;
  JB->rowvals[63] = 11;
  JB->rowvals[64] = 12;
  JB->rowvals[65] = 13;
  JB->rowvals[66] = 15;
  JB->rowvals[67] = 16;
  JB->rowvals[68] = 2;
  JB->rowvals[69] = 3;
  JB->rowvals[70] = 7;
  JB->rowvals[71] = 8;
  JB->rowvals[72] = 9;
  JB->rowvals[73] = 10;
  JB->rowvals[74] = 11;
  JB->rowvals[75] = 12;
  JB->rowvals[76] = 13;
  JB->rowvals[77] = 16;
  JB->rowvals[78] = 6;
  JB->rowvals[79] = 10;
  JB->rowvals[80] = 11;
  JB->rowvals[81] = 13;
  JB->rowvals[82] = 6;
  JB->rowvals[83] = 8;
  JB->rowvals[84] = 10;
  JB->rowvals[85] = 11;
  JB->rowvals[86] = 12;
  JB->rowvals[87] = 13;
  JB->rowvals[88] = 8;
  JB->rowvals[89] = 11;
  JB->rowvals[90] = 12;
  JB->rowvals[91] = 13;
  JB->rowvals[92] = 6;
  JB->rowvals[93] = 8;
  JB->rowvals[94] = 10;
  JB->rowvals[95] = 11;
  JB->rowvals[96] = 12;
  JB->rowvals[97] = 13;
  JB->rowvals[98] = 1;
  JB->rowvals[99] = 5;
  JB->rowvals[100] = 7;
  JB->rowvals[101] = 10;
  JB->rowvals[102] = 13;
  JB->rowvals[103] = 14;
  JB->rowvals[104] = 17;
  JB->rowvals[105] = 19;
  JB->rowvals[106] = 5;
  JB->rowvals[107] = 6;
  JB->rowvals[108] = 8;
  JB->rowvals[109] = 10;
  JB->rowvals[110] = 11;
  JB->rowvals[111] = 12;
  JB->rowvals[112] = 13;
  JB->rowvals[113] = 15;
  JB->rowvals[114] = 16;
  JB->rowvals[115] = 2;
  JB->rowvals[116] = 7;
  JB->rowvals[117] = 8;
  JB->rowvals[118] = 9;
  JB->rowvals[119] = 10;
  JB->rowvals[120] = 11;
  JB->rowvals[121] = 12;
  JB->rowvals[122] = 13;
  JB->rowvals[123] = 15;
  JB->rowvals[124] = 16;
  JB->rowvals[125] = 1;
  JB->rowvals[126] = 10;
  JB->rowvals[127] = 11;
  JB->rowvals[128] = 13;
  JB->rowvals[129] = 14;
  JB->rowvals[130] = 15;
  JB->rowvals[131] = 16;
  JB->rowvals[132] = 17;
  JB->rowvals[133] = 18;
  JB->rowvals[134] = 19;
  JB->rowvals[135] = 1;
  JB->rowvals[136] = 2;
  JB->rowvals[137] = 11;
  JB->rowvals[138] = 13;
  JB->rowvals[139] = 16;
  JB->rowvals[140] = 17;
  JB->rowvals[141] = 18;
  JB->rowvals[142] = 19;
  JB->rowvals[143] = 1;
  JB->rowvals[144] = 8;
  JB->rowvals[145] = 9;
  JB->rowvals[146] = 10;
  JB->rowvals[147] = 11;
  JB->rowvals[148] = 13;
  JB->rowvals[149] = 14;
  JB->rowvals[150] = 17;
  JB->rowvals[151] = 18;
  JB->rowvals[152] = 19;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 6;
  JB->colptrs[2] = 9;
  JB->colptrs[3] = 16;
  JB->colptrs[4] = 24;
  JB->colptrs[5] = 30;
  JB->colptrs[6] = 39;
  JB->colptrs[7] = 48;
  JB->colptrs[8] = 57;
  JB->colptrs[9] = 68;
  JB->colptrs[10] = 78;
  JB->colptrs[11] = 82;
  JB->colptrs[12] = 88;
  JB->colptrs[13] = 92;
  JB->colptrs[14] = 98;
  JB->colptrs[15] = 106;
  JB->colptrs[16] = 115;
  JB->colptrs[17] = 125;
  JB->colptrs[18] = 135;
  JB->colptrs[19] = 143;
  JB->colptrs[20] = 153;
  return(0);
}


 int sx_enhancer_321(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
  memset(sxdot_tmp,0,sizeof(double)*20);
  switch (plist[ip]) {
  case 0: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[17] - 2.0*p[3]*sx_tmp[0] + 4.0*p[2]*sx_tmp[18] + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[11])/k[0] + (2.0*p[2]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[1] + 2.0*p[2]*sx_tmp[12];
sxdot_tmp[2] = p[0]*sx_tmp[8] - 4.0*p[1]*sx_tmp[2] + ((pow(k[0],2))*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[12] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] + 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) - 6.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13])/(pow(k[0],2)) + (sx_tmp[11]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[13]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[11] + 2.0*k[0]*p[0]*sx_tmp[16]*x_tmp[13];
sxdot_tmp[3] = 4.0*p[1]*sx_tmp[9] + (2.0*(pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[8] - 4.0*(pow(k[0],3))*x_tmp[3]*x_tmp[10] - 2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[11] - 4.0*(pow(k[0],3))*x_tmp[6]*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[3]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[13]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2));
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[5] + (2.0*(pow(k[0],2))*x_tmp[6] - 4.0*(pow(k[0],3))*x_tmp[6]*x_tmp[10] - 4.0*(pow(k[0],3))*x_tmp[4]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13])/(pow(k[0],2)) + (p[1]*sx_tmp[11])/k[0] - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[6])/(pow(k[0],2)) - 4.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[15] - (1.0*(2.0*(pow(k[0],2))*x_tmp[6] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[10] - 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[10] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + 3.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + 3.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 3.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 9.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[13]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) - (1.0*sx_tmp[11]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[10]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[8] - (1.0*(2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[10] - 2.0*(pow(k[0],2))*x_tmp[6] + 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[10] + 2.0*(pow(k[0],3))*x_tmp[4]*x_tmp[13] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[10] + 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) + (sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[10] - 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[7] = p[1]*sx_tmp[16] + (2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[10] + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) - 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) - 3.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13] + 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13])/(pow(k[0],2)) - (1.0*sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[7])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[12]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[6]*x_tmp[11] + k[0]*p[0]*sx_tmp[5]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[10];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + 2.0*p[1]*sx_tmp[16] - (1.0*(2.0*(pow(k[0],2))*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[3]*x_tmp[10] + (pow(k[0],3))*x_tmp[3]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[13] + (pow(k[0],3))*x_tmp[8]*x_tmp[11] + (pow(k[0],3))*x_tmp[8]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[11]*x_tmp[13] + (pow(k[0],3))*x_tmp[13]*x_tmp[15] + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + 2.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[11]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[3])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2)) + (sx_tmp[10]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[9] = 2.0*p[1]*sx_tmp[2] - (1.0*((pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],3))*x_tmp[3]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] + 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) - 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[8])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[10] = p[1]*sx_tmp[11] - (1.0*(2.0*(pow(k[0],2))*x_tmp[6] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13]))/k[0] - 2.0*k[0]*p[0]*sx_tmp[6] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[13]*x_tmp[10];
sxdot_tmp[11] = 2.0*p[1]*sx_tmp[12] + (2.0*(pow(k[0],2))*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[11]*x_tmp[13])/k[0] + ((2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*sx_tmp[13])/k[0] + 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - (1.0*(k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] + 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[12] = ((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13])/k[0] - 2.0*p[1]*sx_tmp[12] + k[0]*p[0]*sx_tmp[8] + k[0]*p[0]*sx_tmp[11]*x_tmp[13] + k[0]*p[0]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = 2.0*p[1]*sx_tmp[12] - (1.0*(2.0*(pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]))/k[0] + ((k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] - (1.0*sx_tmp[13]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] + 2.0*p[2]*sx_tmp[7] + p[1]*sx_tmp[17] - (1.0*(2.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 2.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6] + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) - 6.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[13])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[19]*x_tmp[10];
sxdot_tmp[15] = 2.0*p[0]*sx_tmp[6] + 4.0*p[1]*sx_tmp[16] + (2.0*(pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[8] - 4.0*(pow(k[0],3))*x_tmp[6]*x_tmp[11] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[10]*x_tmp[13] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] - 2.0*(pow(k[0],3))*x_tmp[13]*x_tmp[15] + 4.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + 4.0*k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + 4.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 12.0*(pow(k[0],3))*x_tmp[10]*x_tmp[11]*x_tmp[13])/(pow(k[0],2)) + (sx_tmp[13]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] + (sx_tmp[10]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[11]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 4.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[16] = 2.0*p[1]*sx_tmp[2] - (1.0*((pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],3))*x_tmp[6]*x_tmp[12] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[11] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[12] + (pow(k[0],2))*x_tmp[11]*x_tmp[13] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[15] - 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) - 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + 6.0*(pow(k[0],3))*x_tmp[10]*x_tmp[12]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[11]*x_tmp[12]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[8])/(pow(k[0],2)) - (1.0*sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[9])/(pow(k[0],2)) - (1.0*(3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13] + k[0]*p[0]*sx_tmp[15]*x_tmp[13];
sxdot_tmp[17] = p[2]*sx_tmp[15] + 2.0*p[2]*sx_tmp[16] + 2.0*p[1]*sx_tmp[18] - (1.0*(2.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6] - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[8] - 2.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) - 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + 6.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[1]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[19])/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];
sxdot_tmp[18] = 2.0*p[2]*sx_tmp[2] + p[2]*sx_tmp[16] + (k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[8] + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[17]*x_tmp[13] + k[0]*p[0]*sx_tmp[19]*x_tmp[11];
sxdot_tmp[19] = p[2]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + 2.0*p[1]*sx_tmp[18] - (1.0*(2.0*k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*x_tmp[1]*x_tmp[8] - 2.0*(pow(k[0],3))*x_tmp[1]*x_tmp[6] + k[0]*x_tmp[1]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + 2.0*k[0]*x_tmp[10]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + k[0]*x_tmp[11]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) - 6.0*(pow(k[0],3))*x_tmp[1]*x_tmp[10]*x_tmp[13] - 3.0*(pow(k[0],3))*x_tmp[1]*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[19]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*sx_tmp[1])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[17])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];

  } break;

  case 1: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[17] - 2.0*p[3]*sx_tmp[0] + 4.0*p[2]*sx_tmp[18] + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[11])/k[0] + (2.0*p[2]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[1] + 2.0*p[2]*sx_tmp[12];
sxdot_tmp[2] = p[0]*sx_tmp[8] - 4.0*p[1]*sx_tmp[2] + (2.0*k[0]*x_tmp[12] - 4.0*(pow(k[0],2))*x_tmp[2])/(pow(k[0],2)) + (sx_tmp[11]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[13]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[11] + 2.0*k[0]*p[0]*sx_tmp[16]*x_tmp[13];
sxdot_tmp[3] = 4.0*p[1]*sx_tmp[9] + (k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[8] + 4.0*(pow(k[0],2))*x_tmp[9])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[3]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[13]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2));
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[5] + (k[0]*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[5])/(pow(k[0],2)) + (p[1]*sx_tmp[11])/k[0] - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[6])/(pow(k[0],2)) - 4.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[15] - (1.0*(k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) + (sx_tmp[13]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) - (1.0*sx_tmp[11]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[10]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[8] + (k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[8])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) + (sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[10] - 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[7] = p[1]*sx_tmp[16] - (1.0*(2.0*(pow(k[0],2))*x_tmp[7] - 1.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[7])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[12]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[6]*x_tmp[11] + k[0]*p[0]*sx_tmp[5]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[10];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + 2.0*p[1]*sx_tmp[16] + (2.0*k[0]*x_tmp[12] - 1.0*k[0]*x_tmp[11] - 1.0*(pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[16])/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[11]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[3])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2)) + (sx_tmp[10]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[9] = 2.0*p[1]*sx_tmp[2] - (1.0*(2.0*k[0]*x_tmp[12] - 2.0*(pow(k[0],2))*x_tmp[2] + 2.0*(pow(k[0],2))*x_tmp[9] - 1.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[8])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[10] = p[1]*sx_tmp[11] + x_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[6] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[13]*x_tmp[10];
sxdot_tmp[11] = 2.0*p[1]*sx_tmp[12] - (1.0*(k[0]*x_tmp[11] - 2.0*k[0]*x_tmp[12]))/k[0] + ((2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*sx_tmp[13])/k[0] + 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - (1.0*(k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] + 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[12] = k[0]*p[0]*sx_tmp[8] - 2.0*x_tmp[12] - 2.0*p[1]*sx_tmp[12] + k[0]*p[0]*sx_tmp[11]*x_tmp[13] + k[0]*p[0]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = 2.0*p[1]*sx_tmp[12] + (k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12])/k[0] + ((k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] - (1.0*sx_tmp[13]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] + 2.0*p[2]*sx_tmp[7] + p[1]*sx_tmp[17] + x_tmp[17] - (1.0*((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[13])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[19]*x_tmp[10];
sxdot_tmp[15] = 2.0*p[0]*sx_tmp[6] + 4.0*p[1]*sx_tmp[16] + (k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12] - 2.0*(pow(k[0],2))*x_tmp[15] + 4.0*(pow(k[0],2))*x_tmp[16])/(pow(k[0],2)) + (sx_tmp[13]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] + (sx_tmp[10]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[11]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 4.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[16] = 2.0*p[1]*sx_tmp[2] - (1.0*(2.0*k[0]*x_tmp[12] - 2.0*(pow(k[0],2))*x_tmp[2] + 3.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[8])/(pow(k[0],2)) - (1.0*sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[9])/(pow(k[0],2)) - (1.0*(3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13] + k[0]*p[0]*sx_tmp[15]*x_tmp[13];
sxdot_tmp[17] = p[2]*sx_tmp[15] + 2.0*p[2]*sx_tmp[16] + 2.0*p[1]*sx_tmp[18] - (1.0*((pow(k[0],2))*x_tmp[17] - 2.0*(pow(k[0],2))*x_tmp[18]))/(pow(k[0],2)) - (1.0*sx_tmp[1]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[19])/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];
sxdot_tmp[18] = 2.0*p[2]*sx_tmp[2] + p[2]*sx_tmp[16] - 2.0*x_tmp[18] - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[17]*x_tmp[13] + k[0]*p[0]*sx_tmp[19]*x_tmp[11];
sxdot_tmp[19] = p[2]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + 2.0*p[1]*sx_tmp[18] + ((pow(k[0],2))*x_tmp[17] + 2.0*(pow(k[0],2))*x_tmp[18])/(pow(k[0],2)) - (1.0*sx_tmp[19]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*sx_tmp[1])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[17])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];

  } break;

  case 2: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[17] - 2.0*p[3]*sx_tmp[0] + 4.0*p[2]*sx_tmp[18] + (k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[17] + 4.0*(pow(k[0],2))*x_tmp[18])/(pow(k[0],2)) + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[11])/k[0] + (2.0*p[2]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[1] + 2.0*p[2]*sx_tmp[12] + (k[0]*x_tmp[11] + 2.0*k[0]*x_tmp[12])/k[0];
sxdot_tmp[2] = p[0]*sx_tmp[8] - 4.0*p[1]*sx_tmp[2] + (sx_tmp[11]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[13]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[11] + 2.0*k[0]*p[0]*sx_tmp[16]*x_tmp[13];
sxdot_tmp[3] = 4.0*p[1]*sx_tmp[9] + (sx_tmp[8]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[3]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[13]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2));
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[5] + (p[1]*sx_tmp[11])/k[0] - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[6])/(pow(k[0],2)) - 4.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[15] + (sx_tmp[13]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) - (1.0*sx_tmp[11]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[10]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[8] - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) + (sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[10] - 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[7] = p[1]*sx_tmp[16] - (1.0*sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[7])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[12]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[6]*x_tmp[11] + k[0]*p[0]*sx_tmp[5]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[10];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + 2.0*p[1]*sx_tmp[16] + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[11]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[3])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2)) + (sx_tmp[10]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[9] = 2.0*p[1]*sx_tmp[2] + (sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[8])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[10] = p[1]*sx_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[6] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[13]*x_tmp[10];
sxdot_tmp[11] = 2.0*p[1]*sx_tmp[12] + ((2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*sx_tmp[13])/k[0] + 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - (1.0*(k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] + 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[12] = k[0]*p[0]*sx_tmp[8] - 2.0*p[1]*sx_tmp[12] + k[0]*p[0]*sx_tmp[11]*x_tmp[13] + k[0]*p[0]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = 2.0*p[1]*sx_tmp[12] + ((k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] - (1.0*sx_tmp[13]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] + 2.0*p[2]*sx_tmp[7] + p[1]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[7])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[13])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[19]*x_tmp[10];
sxdot_tmp[15] = 2.0*p[0]*sx_tmp[6] + 4.0*p[1]*sx_tmp[16] + (sx_tmp[13]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] + (sx_tmp[10]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[11]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 4.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[16] = 2.0*p[1]*sx_tmp[2] + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[8])/(pow(k[0],2)) - (1.0*sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[9])/(pow(k[0],2)) - (1.0*(3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13] + k[0]*p[0]*sx_tmp[15]*x_tmp[13];
sxdot_tmp[17] = p[2]*sx_tmp[15] + 2.0*p[2]*sx_tmp[16] + 2.0*p[1]*sx_tmp[18] + ((pow(k[0],2))*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[16])/(pow(k[0],2)) - (1.0*sx_tmp[1]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[19])/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];
sxdot_tmp[18] = 2.0*p[2]*sx_tmp[2] + p[2]*sx_tmp[16] + (2.0*(pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[17]*x_tmp[13] + k[0]*p[0]*sx_tmp[19]*x_tmp[11];
sxdot_tmp[19] = p[2]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + 2.0*p[1]*sx_tmp[18] + ((pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[9])/(pow(k[0],2)) - (1.0*sx_tmp[19]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*sx_tmp[1])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[17])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];

  } break;

  case 3: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[17] - 2.0*p[3]*sx_tmp[0] + 4.0*p[2]*sx_tmp[18] + (k[0]*x_tmp[1] - 2.0*(pow(k[0],2))*x_tmp[0])/(pow(k[0],2)) + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[11])/k[0] + (2.0*p[2]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[1] + 2.0*p[2]*sx_tmp[12] - 1.0*x_tmp[1];
sxdot_tmp[2] = p[0]*sx_tmp[8] - 4.0*p[1]*sx_tmp[2] + (sx_tmp[11]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[13]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[11] + 2.0*k[0]*p[0]*sx_tmp[16]*x_tmp[13];
sxdot_tmp[3] = 4.0*p[1]*sx_tmp[9] + (sx_tmp[8]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[3]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[13]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2));
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[5] + (p[1]*sx_tmp[11])/k[0] - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[6])/(pow(k[0],2)) - 4.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[15] + (sx_tmp[13]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) - (1.0*sx_tmp[11]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[10]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[8] - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) + (sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[10] - 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[7] = p[1]*sx_tmp[16] - (1.0*sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[7])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[12]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[6]*x_tmp[11] + k[0]*p[0]*sx_tmp[5]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[10];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + 2.0*p[1]*sx_tmp[16] + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[11]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[3])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2)) + (sx_tmp[10]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[9] = 2.0*p[1]*sx_tmp[2] + (sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[8])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[10] = p[1]*sx_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[6] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[13]*x_tmp[10];
sxdot_tmp[11] = 2.0*p[1]*sx_tmp[12] + ((2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*sx_tmp[13])/k[0] + 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - (1.0*(k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] + 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[12] = k[0]*p[0]*sx_tmp[8] - 2.0*p[1]*sx_tmp[12] + k[0]*p[0]*sx_tmp[11]*x_tmp[13] + k[0]*p[0]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = 2.0*p[1]*sx_tmp[12] + ((k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] - (1.0*sx_tmp[13]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] + 2.0*p[2]*sx_tmp[7] + p[1]*sx_tmp[17] - 1.0*x_tmp[14] - (1.0*((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[13])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[19]*x_tmp[10];
sxdot_tmp[15] = 2.0*p[0]*sx_tmp[6] + 4.0*p[1]*sx_tmp[16] + (sx_tmp[13]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] + (sx_tmp[10]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[11]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 4.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[16] = 2.0*p[1]*sx_tmp[2] + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[8])/(pow(k[0],2)) - (1.0*sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[9])/(pow(k[0],2)) - (1.0*(3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13] + k[0]*p[0]*sx_tmp[15]*x_tmp[13];
sxdot_tmp[17] = p[2]*sx_tmp[15] + 2.0*p[2]*sx_tmp[16] + 2.0*p[1]*sx_tmp[18] - 1.0*x_tmp[17] - (1.0*sx_tmp[1]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[19])/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];
sxdot_tmp[18] = 2.0*p[2]*sx_tmp[2] + p[2]*sx_tmp[16] - 1.0*x_tmp[18] - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[17]*x_tmp[13] + k[0]*p[0]*sx_tmp[19]*x_tmp[11];
sxdot_tmp[19] = p[2]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + 2.0*p[1]*sx_tmp[18] - 1.0*x_tmp[19] - (1.0*sx_tmp[19]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*sx_tmp[1])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[17])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];

  } break;

  case 4: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[17] - 2.0*p[3]*sx_tmp[0] + 4.0*p[2]*sx_tmp[18] + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[11])/k[0] + (2.0*p[2]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[1] + 2.0*p[2]*sx_tmp[12];
sxdot_tmp[2] = p[0]*sx_tmp[8] - 4.0*p[1]*sx_tmp[2] + (sx_tmp[11]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[13]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[11] + 2.0*k[0]*p[0]*sx_tmp[16]*x_tmp[13];
sxdot_tmp[3] = 4.0*p[1]*sx_tmp[9] + 16.0/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[3]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[13]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2));
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[5] + (p[1]*sx_tmp[11])/k[0] - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[6])/(pow(k[0],2)) - 4.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[15] + (sx_tmp[13]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) - (1.0*sx_tmp[11]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[10]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[8] - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) + (sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[10] - 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[7] = p[1]*sx_tmp[16] - (1.0*sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[7])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[12]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[6]*x_tmp[11] + k[0]*p[0]*sx_tmp[5]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[10];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + 2.0*p[1]*sx_tmp[16] + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[11]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[3])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2)) + (sx_tmp[10]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[9] = 2.0*p[1]*sx_tmp[2] + (sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[8])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[10] = p[1]*sx_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[6] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[13]*x_tmp[10];
sxdot_tmp[11] = 2.0*p[1]*sx_tmp[12] + ((2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*sx_tmp[13])/k[0] + 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - (1.0*(k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] + 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[12] = k[0]*p[0]*sx_tmp[8] - 2.0*p[1]*sx_tmp[12] + k[0]*p[0]*sx_tmp[11]*x_tmp[13] + k[0]*p[0]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = 2.0*p[1]*sx_tmp[12] + 4.0/k[0] + ((k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] - (1.0*sx_tmp[13]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] + 2.0*p[2]*sx_tmp[7] + p[1]*sx_tmp[17] - (1.0*((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[13])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[19]*x_tmp[10];
sxdot_tmp[15] = 2.0*p[0]*sx_tmp[6] + 4.0*p[1]*sx_tmp[16] + (sx_tmp[13]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] + (sx_tmp[10]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[11]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 4.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[16] = 2.0*p[1]*sx_tmp[2] + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[8])/(pow(k[0],2)) - (1.0*sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[9])/(pow(k[0],2)) - (1.0*(3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13] + k[0]*p[0]*sx_tmp[15]*x_tmp[13];
sxdot_tmp[17] = p[2]*sx_tmp[15] + 2.0*p[2]*sx_tmp[16] + 2.0*p[1]*sx_tmp[18] - (1.0*sx_tmp[1]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[19])/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];
sxdot_tmp[18] = 2.0*p[2]*sx_tmp[2] + p[2]*sx_tmp[16] - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[17]*x_tmp[13] + k[0]*p[0]*sx_tmp[19]*x_tmp[11];
sxdot_tmp[19] = p[2]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + 2.0*p[1]*sx_tmp[18] - (1.0*sx_tmp[19]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*sx_tmp[1])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[17])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];

  } break;

  case 5: {
sxdot_tmp[0] = 2.0*p[2]*sx_tmp[17] - 2.0*p[3]*sx_tmp[0] + 4.0*p[2]*sx_tmp[18] + (p[3]*sx_tmp[1])/k[0] + (p[2]*sx_tmp[11])/k[0] + (2.0*p[2]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[1] + 2.0*p[2]*sx_tmp[12];
sxdot_tmp[2] = p[0]*sx_tmp[8] - 4.0*p[1]*sx_tmp[2] + (sx_tmp[11]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[13]*((pow(k[0],2))*p[0]*x_tmp[11] - 2.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[11] + 2.0*k[0]*p[0]*sx_tmp[16]*x_tmp[13];
sxdot_tmp[3] = 4.0*p[1]*sx_tmp[9] + (k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[3])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[3]*(2.0*(pow(k[0],2))*p[5] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] + 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[13]*(k[0]*p[5] - 4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(k[0]*p[1] - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2));
sxdot_tmp[4] = 2.0*p[1]*sx_tmp[5] + (p[1]*sx_tmp[11])/k[0] - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(4.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[10])/(pow(k[0],2)) + ((2.0*(pow(k[0],2))*p[0] - 4.0*(pow(k[0],3))*p[0]*x_tmp[10])*sx_tmp[6])/(pow(k[0],2)) - 4.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[15] + (sx_tmp[13]*(3.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 3.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + 3.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) - (1.0*sx_tmp[11]*(k[0]*p[1] - 3.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[10]*(3.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[10];
sxdot_tmp[6] = 2.0*p[1]*sx_tmp[7] + p[1]*sx_tmp[8] - 1.0*x_tmp[6] - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[0]*x_tmp[10]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[5])/(pow(k[0],2)) + (sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + k[0]*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[6]))/(pow(k[0],2)) - (1.0*sx_tmp[6]*((pow(k[0],2))*p[5] - 2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[10] - 2.0*k[0]*p[0]*sx_tmp[4]*x_tmp[13];
sxdot_tmp[7] = p[1]*sx_tmp[16] - (1.0*sx_tmp[10]*(p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[6])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[7])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[7]))/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[12]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[6]*x_tmp[11] + k[0]*p[0]*sx_tmp[5]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[10];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + 2.0*p[1]*sx_tmp[16] - 1.0*x_tmp[8] + (2.0*p[1]*sx_tmp[12])/k[0] - (1.0*sx_tmp[11]*(k[0]*p[1] - 2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[3])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 2.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[6] + (pow(k[0],3))*p[0]*x_tmp[8] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[6])/(pow(k[0],2)) + (sx_tmp[10]*(2.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[8]*((pow(k[0],2))*p[1] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[9] = 2.0*p[1]*sx_tmp[2] - 1.0*x_tmp[9] + (sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) - 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[8])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[3]*x_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13];
sxdot_tmp[10] = p[1]*sx_tmp[11] - 2.0*k[0]*p[0]*sx_tmp[6] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13] - 2.0*k[0]*p[0]*sx_tmp[13]*x_tmp[10];
sxdot_tmp[11] = 2.0*p[1]*sx_tmp[12] + ((2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11])*sx_tmp[13])/k[0] + 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - (1.0*(k[0]*p[1] + (pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] + 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[12] = k[0]*p[0]*sx_tmp[8] - 2.0*p[1]*sx_tmp[12] + k[0]*p[0]*sx_tmp[11]*x_tmp[13] + k[0]*p[0]*sx_tmp[13]*x_tmp[11];
sxdot_tmp[13] = 2.0*p[1]*sx_tmp[12] - 1.0*x_tmp[13] + ((k[0]*p[1] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13])*sx_tmp[11])/k[0] - (1.0*sx_tmp[13]*(k[0]*p[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11]))/k[0] - 2.0*k[0]*p[0]*sx_tmp[6] - 1.0*k[0]*p[0]*sx_tmp[8] - 2.0*k[0]*p[0]*sx_tmp[10]*x_tmp[13];
sxdot_tmp[14] = p[2]*sx_tmp[5] + 2.0*p[2]*sx_tmp[7] + p[1]*sx_tmp[17] - (1.0*((pow(k[0],2))*p[3] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[13])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]))/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[19]*x_tmp[10];
sxdot_tmp[15] = 2.0*p[0]*sx_tmp[6] + 4.0*p[1]*sx_tmp[16] + (sx_tmp[13]*(4.0*(pow(k[0],3))*p[0]*x_tmp[5] - 4.0*p[0]*((pow(k[0],3))*x_tmp[5] - 1.0*k[0]*((pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[10]*x_tmp[11]) + (pow(k[0],3))*x_tmp[10]*x_tmp[11]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[11] - 2.0*(pow(k[0],3))*p[0]*x_tmp[15]))/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[12])/k[0] + (sx_tmp[10]*(4.0*(pow(k[0],3))*p[0]*x_tmp[8] - 4.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*sx_tmp[11]*(4.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*k[0]*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + 2.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[15])/(pow(k[0],2)) + (sx_tmp[8]*((pow(k[0],2))*p[0] + 4.0*(pow(k[0],3))*p[0]*x_tmp[10] - 2.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + 4.0*k[0]*p[0]*sx_tmp[5]*x_tmp[13];
sxdot_tmp[16] = 2.0*p[1]*sx_tmp[2] + (sx_tmp[11]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[13]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[8])/(pow(k[0],2)) - (1.0*sx_tmp[12]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]) + 2.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[7] - 1.0*k[0]*((pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[10]*x_tmp[12]) + (pow(k[0],3))*x_tmp[10]*x_tmp[12]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[16] - 1.0*k[0]*((pow(k[0],2))*x_tmp[16] + (pow(k[0],2))*x_tmp[11]*x_tmp[12]) + (pow(k[0],3))*x_tmp[11]*x_tmp[12]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[7] + (pow(k[0],2))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[9])/(pow(k[0],2)) - (1.0*(3.0*(pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[9])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[7]*x_tmp[13] + k[0]*p[0]*sx_tmp[15]*x_tmp[13];
sxdot_tmp[17] = p[2]*sx_tmp[15] + 2.0*p[2]*sx_tmp[16] + 2.0*p[1]*sx_tmp[18] - (1.0*sx_tmp[1]*(2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13])))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[10] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[19])/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[13]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];
sxdot_tmp[18] = 2.0*p[2]*sx_tmp[2] + p[2]*sx_tmp[16] - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[17]*x_tmp[13] + k[0]*p[0]*sx_tmp[19]*x_tmp[11];
sxdot_tmp[19] = p[2]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + 2.0*p[1]*sx_tmp[18] - 1.0*x_tmp[19] - (1.0*sx_tmp[19]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[10] + (pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[6] - 1.0*k[0]*((pow(k[0],2))*x_tmp[6] + (pow(k[0],2))*x_tmp[10]*x_tmp[13]) + (pow(k[0],3))*x_tmp[10]*x_tmp[13]) + p[0]*((pow(k[0],3))*x_tmp[8] - 1.0*k[0]*((pow(k[0],2))*x_tmp[8] + (pow(k[0],2))*x_tmp[11]*x_tmp[13]) + (pow(k[0],3))*x_tmp[11]*x_tmp[13]))*sx_tmp[1])/(pow(k[0],2)) + (((pow(k[0],2))*p[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[13])*sx_tmp[17])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[11])/(pow(k[0],2)) + (sx_tmp[13]*(2.0*p[0]*((pow(k[0],3))*x_tmp[14] - 1.0*k[0]*((pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[1]*x_tmp[10]) + (pow(k[0],3))*x_tmp[1]*x_tmp[10]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[1]*x_tmp[11]) + (pow(k[0],3))*x_tmp[1]*x_tmp[11]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[19] - 1.0*k[0]*((pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[1]*x_tmp[13]) + (pow(k[0],3))*x_tmp[1]*x_tmp[13]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[19])*sx_tmp[10])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[14]*x_tmp[13];

  } break;

  }
 for (ix=0; ix<20; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_enhancer_321(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  memset(sx0_tmp,0,sizeof(double)*20);
  switch (ip) {
  }

  return;
}


void y_enhancer_321(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*1];
y[it+nt*1] = x[it+nt*0];
    
    return;
}


void dydp_enhancer_321(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }
  
  return;
}


void dydx_enhancer_321(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*40);
dydx[1] = 1.0;
dydx[2] = 1.0;
  
  return;
}


void sy_enhancer_321(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
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
int root_enhancer_321(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_enhancer_321(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root_enhancer_321(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval_enhancer_321(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval_enhancer_321(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
void deltadisc_enhancer_321(double t, int idisc, N_Vector x, void *user_data){
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double deltadisc[20];
  memset(deltadisc,0,sizeof(double)*20);
  for(ix = 0; ix<20;ix++){;
  x_tmp[ix] += deltadisc[ix];
  };
}
void sdeltadisc_enhancer_321(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data){
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
  double deltadisc[20];
  double *sdeltadisc;
  memset(deltadisc,0,sizeof(double)*20);
  sdeltadisc = mxMalloc(sizeof(double)*20*np);
  memset(sdeltadisc,0,sizeof(double)*20*np);
  for (ip=0; ip<np; ip++) {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
     switch (plist[ip]) {
     }
  }
  for(ip = 0; ip<np;ip++){
      sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
      for(ix = 0; ix<20;ix++){
      sx_tmp[ix] += sdeltadisc[plist[ip]+np*ix];
     }
  }
  for(ix = 0; ix<20;ix++){
  x_tmp[ix] += deltadisc[ix];
  };
 mxFree(sdeltadisc);
}


void dxdotdp_enhancer_321(double t, int nt, int it, double *dxdotdp, double *p, double *k, double *u, double *x, int *plist, int np, int nx){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
dxdotdp[(2+ip*nx)] = ((pow(k[0],2))*x[it+nt*8] - 2.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*12] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13] + 2.0*k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) + 2.0*k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*13]) + 2.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*16] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*12]) - 6.0*(pow(k[0],3))*x[it+nt*11]*x[it+nt*12]*x[it+nt*13])/(pow(k[0],2));
dxdotdp[(3+ip*nx)] = (2.0*(pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*8] - 4.0*(pow(k[0],3))*x[it+nt*3]*x[it+nt*10] - 2.0*(pow(k[0],3))*x[it+nt*3]*x[it+nt*11] - 4.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*13] - 2.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*13] + 2.0*(pow(k[0],2))*x[it+nt*10]*x[it+nt*13] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13])/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = (2.0*(pow(k[0],2))*x[it+nt*6] - 4.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*10] - 4.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*13] + 2.0*(pow(k[0],2))*x[it+nt*10]*x[it+nt*13])/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = -(1.0*(2.0*(pow(k[0],2))*x[it+nt*6] - 2.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*10] - 2.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*13] - 2.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*11] - 1.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*10] + 2.0*(pow(k[0],2))*x[it+nt*10]*x[it+nt*13] + 3.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*5] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*11]) + 3.0*k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) + 3.0*k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) - 9.0*(pow(k[0],3))*x[it+nt*10]*x[it+nt*11]*x[it+nt*13]))/(pow(k[0],2));
dxdotdp[(6+ip*nx)] = -(1.0*(2.0*(pow(k[0],3))*x[it+nt*3]*x[it+nt*10] - 2.0*(pow(k[0],2))*x[it+nt*6] + 2.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*10] + 2.0*(pow(k[0],3))*x[it+nt*4]*x[it+nt*13] - 1.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*10] + 2.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*13] - 2.0*(pow(k[0],2))*x[it+nt*10]*x[it+nt*13] + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*5] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*11]) + k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) + k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) - 3.0*(pow(k[0],3))*x[it+nt*10]*x[it+nt*11]*x[it+nt*13]))/(pow(k[0],2));
dxdotdp[(7+ip*nx)] = (2.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*12] - 1.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*10] + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*5] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*11]) + k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) - 2.0*k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) + k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) - 2.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*7] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*12]) - 2.0*k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*13]) - 3.0*(pow(k[0],3))*x[it+nt*10]*x[it+nt*11]*x[it+nt*13] + 6.0*(pow(k[0],3))*x[it+nt*10]*x[it+nt*12]*x[it+nt*13])/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = -(1.0*(2.0*(pow(k[0],2))*x[it+nt*6] - 1.0*(pow(k[0],2))*x[it+nt*8] - 2.0*(pow(k[0],3))*x[it+nt*3]*x[it+nt*10] + (pow(k[0],3))*x[it+nt*3]*x[it+nt*11] - 2.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*11] - 2.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*13] + (pow(k[0],3))*x[it+nt*8]*x[it+nt*11] + (pow(k[0],3))*x[it+nt*8]*x[it+nt*13] + 2.0*(pow(k[0],2))*x[it+nt*10]*x[it+nt*13] - 1.0*(pow(k[0],2))*x[it+nt*11]*x[it+nt*13] + (pow(k[0],3))*x[it+nt*13]*x[it+nt*15] + 2.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*5] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*11]) + 2.0*k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) + 2.0*k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) - 6.0*(pow(k[0],3))*x[it+nt*10]*x[it+nt*11]*x[it+nt*13]))/(pow(k[0],2));
dxdotdp[(9+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*8] - 1.0*(pow(k[0],3))*x[it+nt*3]*x[it+nt*11] - 2.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*12] - 1.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*12] - 1.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*13] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13] + 2.0*k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) + 2.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*7] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*12]) + 2.0*k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*13]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) + k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*13]) + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*16] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*12]) - 6.0*(pow(k[0],3))*x[it+nt*10]*x[it+nt*12]*x[it+nt*13] - 3.0*(pow(k[0],3))*x[it+nt*11]*x[it+nt*12]*x[it+nt*13]))/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = -(1.0*(2.0*(pow(k[0],2))*x[it+nt*6] + 2.0*(pow(k[0],2))*x[it+nt*10]*x[it+nt*13]))/k[0];
dxdotdp[(11+ip*nx)] = (2.0*(pow(k[0],2))*x[it+nt*6] - 1.0*(pow(k[0],2))*x[it+nt*8] + 2.0*(pow(k[0],2))*x[it+nt*10]*x[it+nt*13] - 1.0*(pow(k[0],2))*x[it+nt*11]*x[it+nt*13])/k[0];
dxdotdp[(12+ip*nx)] = ((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13])/k[0];
dxdotdp[(13+ip*nx)] = -(1.0*(2.0*(pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*8] + 2.0*(pow(k[0],2))*x[it+nt*10]*x[it+nt*13] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]))/k[0];
dxdotdp[(14+ip*nx)] = -(1.0*(2.0*k[0]*x[it+nt*1]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) - 2.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*6] + 2.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*14] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*10]) + 2.0*k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*19] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*13]) - 6.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*10]*x[it+nt*13]))/(pow(k[0],2));
dxdotdp[(15+ip*nx)] = (2.0*(pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*8] - 4.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*11] - 2.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*11] + 2.0*(pow(k[0],2))*x[it+nt*10]*x[it+nt*13] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13] - 2.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*15] + 4.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*5] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*11]) + 4.0*k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) + 4.0*k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) - 12.0*(pow(k[0],3))*x[it+nt*10]*x[it+nt*11]*x[it+nt*13])/(pow(k[0],2));
dxdotdp[(16+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*8] + 2.0*(pow(k[0],3))*x[it+nt*6]*x[it+nt*12] - 1.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*11] - 1.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*12] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13] - 1.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*15] - 2.0*k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) - 2.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*7] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*12]) - 2.0*k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*13]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) + k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*13]) + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*16] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*12]) + 6.0*(pow(k[0],3))*x[it+nt*10]*x[it+nt*12]*x[it+nt*13] - 3.0*(pow(k[0],3))*x[it+nt*11]*x[it+nt*12]*x[it+nt*13]))/(pow(k[0],2));
dxdotdp[(17+ip*nx)] = -(1.0*(2.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*6] - 1.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*8] - 2.0*k[0]*x[it+nt*1]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) + k[0]*x[it+nt*1]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) - 2.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*14] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*10]) + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*17] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*11]) - 2.0*k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*19] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*13]) + k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*19] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*13]) + 6.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*10]*x[it+nt*13] - 3.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*11]*x[it+nt*13]))/(pow(k[0],2));
dxdotdp[(18+ip*nx)] = (k[0]*x[it+nt*1]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) - 1.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*8] + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*17] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*11]) + k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*19] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*13]) - 3.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*11]*x[it+nt*13])/(pow(k[0],2));
dxdotdp[(19+ip*nx)] = -(1.0*(2.0*k[0]*x[it+nt*1]*((pow(k[0],2))*x[it+nt*6] + (pow(k[0],2))*x[it+nt*10]*x[it+nt*13]) - 1.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*8] - 2.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*6] + k[0]*x[it+nt*1]*((pow(k[0],2))*x[it+nt*8] + (pow(k[0],2))*x[it+nt*11]*x[it+nt*13]) + 2.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*14] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*10]) + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*17] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*11]) + 2.0*k[0]*x[it+nt*10]*((pow(k[0],2))*x[it+nt*19] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*13]) + k[0]*x[it+nt*11]*((pow(k[0],2))*x[it+nt*19] + (pow(k[0],2))*x[it+nt*1]*x[it+nt*13]) - 6.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*10]*x[it+nt*13] - 3.0*(pow(k[0],3))*x[it+nt*1]*x[it+nt*11]*x[it+nt*13]))/(pow(k[0],2));

  } break;

  case 1: {
dxdotdp[(2+ip*nx)] = (2.0*k[0]*x[it+nt*12] - 4.0*(pow(k[0],2))*x[it+nt*2])/(pow(k[0],2));
dxdotdp[(3+ip*nx)] = (k[0]*x[it+nt*11] + 2.0*k[0]*x[it+nt*12] + 2.0*(pow(k[0],2))*x[it+nt*8] + 4.0*(pow(k[0],2))*x[it+nt*9])/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = (k[0]*x[it+nt*11] + 2.0*(pow(k[0],2))*x[it+nt*5])/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = -(1.0*(k[0]*x[it+nt*11] + (pow(k[0],2))*x[it+nt*5] - 2.0*(pow(k[0],2))*x[it+nt*7] - 1.0*(pow(k[0],2))*x[it+nt*15]))/(pow(k[0],2));
dxdotdp[(6+ip*nx)] = (k[0]*x[it+nt*11] + (pow(k[0],2))*x[it+nt*5] + 2.0*(pow(k[0],2))*x[it+nt*7] + (pow(k[0],2))*x[it+nt*8])/(pow(k[0],2));
dxdotdp[(7+ip*nx)] = -(1.0*(2.0*(pow(k[0],2))*x[it+nt*7] - 1.0*(pow(k[0],2))*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = (2.0*k[0]*x[it+nt*12] - 1.0*k[0]*x[it+nt*11] - 1.0*(pow(k[0],2))*x[it+nt*8] + 2.0*(pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*15] + 2.0*(pow(k[0],2))*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(9+ip*nx)] = -(1.0*(2.0*k[0]*x[it+nt*12] - 2.0*(pow(k[0],2))*x[it+nt*2] + 2.0*(pow(k[0],2))*x[it+nt*9] - 1.0*(pow(k[0],2))*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = x[it+nt*11];
dxdotdp[(11+ip*nx)] = -(1.0*(k[0]*x[it+nt*11] - 2.0*k[0]*x[it+nt*12]))/k[0];
dxdotdp[(12+ip*nx)] = -2.0*x[it+nt*12];
dxdotdp[(13+ip*nx)] = (k[0]*x[it+nt*11] + 2.0*k[0]*x[it+nt*12])/k[0];
dxdotdp[(14+ip*nx)] = x[it+nt*17];
dxdotdp[(15+ip*nx)] = (k[0]*x[it+nt*11] + 2.0*k[0]*x[it+nt*12] - 2.0*(pow(k[0],2))*x[it+nt*15] + 4.0*(pow(k[0],2))*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(16+ip*nx)] = -(1.0*(2.0*k[0]*x[it+nt*12] - 2.0*(pow(k[0],2))*x[it+nt*2] + 3.0*(pow(k[0],2))*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(17+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*17] - 2.0*(pow(k[0],2))*x[it+nt*18]))/(pow(k[0],2));
dxdotdp[(18+ip*nx)] = -2.0*x[it+nt*18];
dxdotdp[(19+ip*nx)] = ((pow(k[0],2))*x[it+nt*17] + 2.0*(pow(k[0],2))*x[it+nt*18])/(pow(k[0],2));

  } break;

  case 2: {
dxdotdp[(0+ip*nx)] = (k[0]*x[it+nt*11] + 2.0*k[0]*x[it+nt*12] + 2.0*(pow(k[0],2))*x[it+nt*17] + 4.0*(pow(k[0],2))*x[it+nt*18])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = (k[0]*x[it+nt*11] + 2.0*k[0]*x[it+nt*12])/k[0];
dxdotdp[(14+ip*nx)] = ((pow(k[0],2))*x[it+nt*5] + 2.0*(pow(k[0],2))*x[it+nt*7])/(pow(k[0],2));
dxdotdp[(17+ip*nx)] = ((pow(k[0],2))*x[it+nt*15] + 2.0*(pow(k[0],2))*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(18+ip*nx)] = (2.0*(pow(k[0],2))*x[it+nt*2] + (pow(k[0],2))*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(19+ip*nx)] = ((pow(k[0],2))*x[it+nt*8] + 2.0*(pow(k[0],2))*x[it+nt*9])/(pow(k[0],2));

  } break;

  case 3: {
dxdotdp[(0+ip*nx)] = (k[0]*x[it+nt*1] - 2.0*(pow(k[0],2))*x[it+nt*0])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = -1.0*x[it+nt*1];
dxdotdp[(14+ip*nx)] = -1.0*x[it+nt*14];
dxdotdp[(17+ip*nx)] = -1.0*x[it+nt*17];
dxdotdp[(18+ip*nx)] = -1.0*x[it+nt*18];
dxdotdp[(19+ip*nx)] = -1.0*x[it+nt*19];

  } break;

  case 4: {
dxdotdp[(3+ip*nx)] = 16.0/(pow(k[0],2));
dxdotdp[(13+ip*nx)] = 4.0/k[0];

  } break;

  case 5: {
dxdotdp[(3+ip*nx)] = (k[0]*x[it+nt*13] - 2.0*(pow(k[0],2))*x[it+nt*3])/(pow(k[0],2));
dxdotdp[(6+ip*nx)] = -1.0*x[it+nt*6];
dxdotdp[(8+ip*nx)] = -1.0*x[it+nt*8];
dxdotdp[(9+ip*nx)] = -1.0*x[it+nt*9];
dxdotdp[(13+ip*nx)] = -1.0*x[it+nt*13];
dxdotdp[(19+ip*nx)] = -1.0*x[it+nt*19];

  } break;

  }
  }
  
  return;
}
