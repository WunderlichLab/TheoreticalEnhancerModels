#include "enhancer_3.h"
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


 int xdot_enhancer_3(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int ix;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  memset(xdot_tmp,0,sizeof(double)*27);
xdot_tmp[0] = (2.0*p[0]*((pow(k[0],3))*x_tmp[0]*x_tmp[16] + (pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[21] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]) + (pow(k[0],2))*p[1]*x_tmp[1] + (pow(k[0],2))*p[1]*x_tmp[17] - 2.0*(pow(k[0],3))*p[0]*x_tmp[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[21])/(pow(k[0],2));
xdot_tmp[1] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[0]*x_tmp[16] + (pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[21] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[16] + (pow(k[0],3))*x_tmp[13]*x_tmp[21] + (pow(k[0],3))*x_tmp[14]*x_tmp[26] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]) + (pow(k[0],2))*p[1]*x_tmp[1] - 1.0*(pow(k[0],2))*p[1]*x_tmp[19] - 1.0*(pow(k[0],3))*p[0]*x_tmp[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[1]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[21] + (pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[26]))/(pow(k[0],2));
xdot_tmp[2] = (2.0*(pow(k[0],2))*p[1]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[21] + k[0]*p[1]*x_tmp[15] - 2.0*(pow(k[0],3))*p[0]*x_tmp[2]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[21])/(pow(k[0],2));
xdot_tmp[3] = -(1.0*((pow(k[0],2))*p[1]*x_tmp[3] - 1.0*p[0]*((pow(k[0],3))*x_tmp[3]*x_tmp[16] + (pow(k[0],3))*x_tmp[15]*x_tmp[21] + (pow(k[0],3))*x_tmp[14]*x_tmp[23] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) - 1.0*k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[1]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[21] + k[0]*p[1]*x_tmp[15] - 1.0*(pow(k[0],3))*p[0]*x_tmp[2]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[21] + (pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[23]))/(pow(k[0],2));
xdot_tmp[4] = ((pow(k[0],2))*p[0]*x_tmp[21] - 2.0*(pow(k[0],2))*p[1]*x_tmp[4] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3]*x_tmp[16] + (pow(k[0],3))*x_tmp[15]*x_tmp[21] + (pow(k[0],3))*x_tmp[14]*x_tmp[23] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) - 1.0*k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]) + k[0]*p[1]*x_tmp[15] + 2.0*(pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[14]*x_tmp[16] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[23])/(pow(k[0],2));
xdot_tmp[5] = (144.0*p[4] - 2.0*(pow(k[0],2))*p[5]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[18] + (pow(k[0],2))*p[0]*x_tmp[21] + 2.0*(pow(k[0],2))*p[1]*x_tmp[23] + 2.0*(pow(k[0],2))*p[1]*x_tmp[26] + k[0]*p[1]*x_tmp[13] + k[0]*p[1]*x_tmp[15] + k[0]*p[5]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[12] - 2.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[14] + (pow(k[0],2))*p[0]*x_tmp[12]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[21])/(pow(k[0],2));
xdot_tmp[6] = (2.0*(pow(k[0],2))*p[2]*x_tmp[20] - 2.0*(pow(k[0],2))*p[3]*x_tmp[6] + 2.0*(pow(k[0],2))*p[2]*x_tmp[24] + k[0]*p[3]*x_tmp[7] + k[0]*p[2]*x_tmp[13] + k[0]*p[2]*x_tmp[15])/(pow(k[0],2));
xdot_tmp[7] = (k[0]*p[2]*x_tmp[13] - 1.0*k[0]*p[3]*x_tmp[7] + k[0]*p[2]*x_tmp[15])/k[0];
xdot_tmp[8] = (2.0*(pow(k[0],2))*p[1]*x_tmp[9] + (pow(k[0],2))*p[0]*x_tmp[18] + k[0]*p[1]*x_tmp[13] - 2.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[12]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[18])/(pow(k[0],2));
xdot_tmp[9] = -(1.0*((pow(k[0],2))*p[1]*x_tmp[9] - 1.0*p[0]*((pow(k[0],3))*x_tmp[9]*x_tmp[16] + (pow(k[0],3))*x_tmp[13]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[26] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[1]*x_tmp[10] + (pow(k[0],2))*p[0]*x_tmp[18] + k[0]*p[1]*x_tmp[13] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[9]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[12]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[18] + (pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[26]))/(pow(k[0],2));
xdot_tmp[10] = ((pow(k[0],2))*p[0]*x_tmp[18] - 2.0*(pow(k[0],2))*p[1]*x_tmp[10] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9]*x_tmp[16] + (pow(k[0],3))*x_tmp[13]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[26] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16]) + k[0]*p[1]*x_tmp[13] + 2.0*(pow(k[0],3))*p[0]*x_tmp[9]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[12]*x_tmp[16] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[26])/(pow(k[0],2));
xdot_tmp[11] = (p[0]*((pow(k[0],3))*x_tmp[7]*x_tmp[18] + (pow(k[0],3))*x_tmp[11]*x_tmp[16] + (pow(k[0],3))*x_tmp[12]*x_tmp[25] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16]) + (pow(k[0],2))*p[2]*x_tmp[9] - 1.0*(pow(k[0],2))*p[3]*x_tmp[11] + (pow(k[0],2))*p[2]*x_tmp[17] + (pow(k[0],2))*p[1]*x_tmp[20] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[25])/(pow(k[0],2));
xdot_tmp[12] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[18] - 1.0*k[0]*p[1]*x_tmp[13] + (pow(k[0],2))*p[0]*x_tmp[12]*x_tmp[16]))/k[0];
xdot_tmp[13] = ((pow(k[0],2))*p[0]*x_tmp[18] - 1.0*k[0]*p[1]*x_tmp[13] + (pow(k[0],2))*p[0]*x_tmp[12]*x_tmp[16])/k[0];
xdot_tmp[14] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[21] - 1.0*k[0]*p[1]*x_tmp[15] + (pow(k[0],2))*p[0]*x_tmp[14]*x_tmp[16]))/k[0];
xdot_tmp[15] = ((pow(k[0],2))*p[0]*x_tmp[21] - 1.0*k[0]*p[1]*x_tmp[15] + (pow(k[0],2))*p[0]*x_tmp[14]*x_tmp[16])/k[0];
xdot_tmp[16] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[18] - 12.0*p[4] + (pow(k[0],2))*p[0]*x_tmp[21] - 1.0*k[0]*p[1]*x_tmp[13] - 1.0*k[0]*p[1]*x_tmp[15] + k[0]*p[5]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[12]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[14]*x_tmp[16]))/k[0];
xdot_tmp[17] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[0]*x_tmp[16] + (pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[21] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[15]*x_tmp[18] + (pow(k[0],3))*x_tmp[16]*x_tmp[17] + (pow(k[0],3))*x_tmp[12]*x_tmp[23] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) - 1.0*k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16]) + (pow(k[0],2))*p[1]*x_tmp[17] - 1.0*(pow(k[0],2))*p[1]*x_tmp[19] - 1.0*(pow(k[0],3))*p[0]*x_tmp[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[23]))/(pow(k[0],2));
xdot_tmp[18] = (p[0]*((pow(k[0],3))*x_tmp[0]*x_tmp[16] + (pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[21] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]) + (pow(k[0],2))*p[1]*x_tmp[9] + (pow(k[0],2))*p[0]*x_tmp[18] + (pow(k[0],2))*p[1]*x_tmp[17] - 1.0*(pow(k[0],2))*p[5]*x_tmp[18] + (pow(k[0],2))*p[1]*x_tmp[26] + k[0]*p[1]*x_tmp[13] - 1.0*(pow(k[0],3))*p[0]*x_tmp[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[12] - 1.0*(pow(k[0],3))*p[0]*x_tmp[8]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[12]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[18])/(pow(k[0],2));
xdot_tmp[19] = ((pow(k[0],3))*p[0]*x_tmp[1]*x_tmp[16] - 1.0*p[0]*((pow(k[0],3))*x_tmp[15]*x_tmp[18] + (pow(k[0],3))*x_tmp[16]*x_tmp[17] + (pow(k[0],3))*x_tmp[12]*x_tmp[23] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) - 1.0*k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16]) - 2.0*(pow(k[0],2))*p[1]*x_tmp[19] - 1.0*p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[16] + (pow(k[0],3))*x_tmp[13]*x_tmp[21] + (pow(k[0],3))*x_tmp[14]*x_tmp[26] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[23] + (pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[26])/(pow(k[0],2));
xdot_tmp[20] = ((pow(k[0],2))*p[2]*x_tmp[10] - 1.0*p[0]*((pow(k[0],3))*x_tmp[7]*x_tmp[18] + (pow(k[0],3))*x_tmp[11]*x_tmp[16] + (pow(k[0],3))*x_tmp[12]*x_tmp[25] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[1]*x_tmp[20] + (pow(k[0],2))*p[2]*x_tmp[19] - 1.0*(pow(k[0],2))*p[3]*x_tmp[20] + (pow(k[0],3))*p[0]*x_tmp[11]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[25])/(pow(k[0],2));
xdot_tmp[21] = (p[0]*((pow(k[0],3))*x_tmp[0]*x_tmp[16] + (pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[21] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]) + (pow(k[0],2))*p[1]*x_tmp[1] + (pow(k[0],2))*p[1]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[21] + (pow(k[0],2))*p[1]*x_tmp[23] - 1.0*(pow(k[0],2))*p[5]*x_tmp[21] + k[0]*p[1]*x_tmp[15] - 1.0*(pow(k[0],3))*p[0]*x_tmp[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[2]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[14] + (pow(k[0],2))*p[0]*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[21])/(pow(k[0],2));
xdot_tmp[22] = (p[0]*((pow(k[0],3))*x_tmp[7]*x_tmp[21] + (pow(k[0],3))*x_tmp[16]*x_tmp[22] + (pow(k[0],3))*x_tmp[14]*x_tmp[25] - 1.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]) + (pow(k[0],2))*p[2]*x_tmp[1] + (pow(k[0],2))*p[2]*x_tmp[3] + (pow(k[0],2))*p[1]*x_tmp[24] - 1.0*(pow(k[0],2))*p[3]*x_tmp[22] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[22] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[25])/(pow(k[0],2));
xdot_tmp[23] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[21] - 1.0*p[0]*((pow(k[0],3))*x_tmp[15]*x_tmp[18] + (pow(k[0],3))*x_tmp[16]*x_tmp[17] + (pow(k[0],3))*x_tmp[12]*x_tmp[23] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) - 1.0*k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[1]*x_tmp[4] - 1.0*(pow(k[0],2))*p[1]*x_tmp[19] - 1.0*p[0]*((pow(k[0],3))*x_tmp[3]*x_tmp[16] + (pow(k[0],3))*x_tmp[15]*x_tmp[21] + (pow(k[0],3))*x_tmp[14]*x_tmp[23] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) - 1.0*k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]) + (pow(k[0],2))*p[1]*x_tmp[23] + (pow(k[0],2))*p[5]*x_tmp[23] + k[0]*p[1]*x_tmp[15] + (pow(k[0],3))*p[0]*x_tmp[3]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[14] + (pow(k[0],2))*p[0]*x_tmp[14]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[23] + (pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[23] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[21]))/(pow(k[0],2));
xdot_tmp[24] = ((pow(k[0],2))*p[2]*x_tmp[4] - 1.0*p[0]*((pow(k[0],3))*x_tmp[7]*x_tmp[21] + (pow(k[0],3))*x_tmp[16]*x_tmp[22] + (pow(k[0],3))*x_tmp[14]*x_tmp[25] - 1.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]) + (pow(k[0],2))*p[2]*x_tmp[19] - 1.0*(pow(k[0],2))*p[1]*x_tmp[24] - 1.0*(pow(k[0],2))*p[3]*x_tmp[24] + (pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[22] + (pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[25])/(pow(k[0],2));
xdot_tmp[25] = (p[0]*((pow(k[0],3))*x_tmp[7]*x_tmp[18] + (pow(k[0],3))*x_tmp[11]*x_tmp[16] + (pow(k[0],3))*x_tmp[12]*x_tmp[25] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[7]*x_tmp[21] + (pow(k[0],3))*x_tmp[16]*x_tmp[22] + (pow(k[0],3))*x_tmp[14]*x_tmp[25] - 1.0*k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]) + (pow(k[0],2))*p[1]*x_tmp[20] + (pow(k[0],2))*p[1]*x_tmp[24] + (pow(k[0],2))*p[2]*x_tmp[23] + (pow(k[0],2))*p[2]*x_tmp[26] - 1.0*(pow(k[0],2))*p[3]*x_tmp[25] - 1.0*(pow(k[0],2))*p[5]*x_tmp[25] - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[25] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[22] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[25])/(pow(k[0],2));
xdot_tmp[26] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[18] - 1.0*p[0]*((pow(k[0],3))*x_tmp[9]*x_tmp[16] + (pow(k[0],3))*x_tmp[13]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[26] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[1]*x_tmp[10] - 1.0*p[0]*((pow(k[0],3))*x_tmp[1]*x_tmp[16] + (pow(k[0],3))*x_tmp[13]*x_tmp[21] + (pow(k[0],3))*x_tmp[14]*x_tmp[26] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) - 1.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[1]*x_tmp[19] + (pow(k[0],2))*p[1]*x_tmp[26] + (pow(k[0],2))*p[5]*x_tmp[26] + k[0]*p[1]*x_tmp[13] + (pow(k[0],3))*p[0]*x_tmp[1]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[5]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[9]*x_tmp[16] + (pow(k[0],2))*p[0]*x_tmp[12]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]*x_tmp[18] + (pow(k[0],3))*p[0]*x_tmp[12]*x_tmp[26] + (pow(k[0],3))*p[0]*x_tmp[14]*x_tmp[26]))/(pow(k[0],2));

  for (ix=0; ix<27; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_enhancer_3(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*54);
xBdot_tmp[0] = 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[0] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[1] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[17] + k[0]*p[0]*x_tmp[16]*xB_tmp[18] + k[0]*p[0]*x_tmp[16]*xB_tmp[21];
xBdot_tmp[1] = (xB_tmp[1]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[21] - 1.0*p[2]*xB_tmp[22] - 1.0*p[1]*xB_tmp[0] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[19] + k[0]*p[0]*x_tmp[16]*xB_tmp[26];
xBdot_tmp[2] = 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[2] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[3] + k[0]*p[0]*x_tmp[16]*xB_tmp[21];
xBdot_tmp[3] = (xB_tmp[3]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[21] - 1.0*p[2]*xB_tmp[22] - 2.0*p[1]*xB_tmp[2] - 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[4] + k[0]*p[0]*x_tmp[16]*xB_tmp[23];
xBdot_tmp[4] = 2.0*p[1]*xB_tmp[4] - 1.0*p[1]*xB_tmp[3] - 1.0*p[1]*xB_tmp[23] - 1.0*p[2]*xB_tmp[24];
xBdot_tmp[5] = (xB_tmp[5]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[12]*xB_tmp[18] + k[0]*p[0]*x_tmp[14]*xB_tmp[21] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[23] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[26];
xBdot_tmp[6] = 2.0*p[3]*xB_tmp[6];
xBdot_tmp[7] = p[3]*xB_tmp[7] - (1.0*p[3]*xB_tmp[6])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*xB_tmp[25])/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[11]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*xB_tmp[20]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[22]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*xB_tmp[24]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
xBdot_tmp[8] = 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[8] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[9] + k[0]*p[0]*x_tmp[16]*xB_tmp[18];
xBdot_tmp[9] = (xB_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[11] - 1.0*p[1]*xB_tmp[18] - 2.0*p[1]*xB_tmp[8] - 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[10] + k[0]*p[0]*x_tmp[16]*xB_tmp[26];
xBdot_tmp[10] = 2.0*p[1]*xB_tmp[10] - 1.0*p[1]*xB_tmp[9] - 1.0*p[2]*xB_tmp[20] - 1.0*p[1]*xB_tmp[26];
xBdot_tmp[11] = (xB_tmp[11]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[20] + k[0]*p[0]*x_tmp[16]*xB_tmp[25];
xBdot_tmp[12] = (xB_tmp[17]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[10]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*xB_tmp[18]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*xB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*xB_tmp[26]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*xB_tmp[8]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*xB_tmp[5])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*xB_tmp[1])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[11])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*xB_tmp[21])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[25])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*xB_tmp[19])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*xB_tmp[23])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*xB_tmp[0])/(pow(k[0],2)) + k[0]*p[0]*x_tmp[16]*xB_tmp[12] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[13] + k[0]*p[0]*x_tmp[16]*xB_tmp[16];
xBdot_tmp[13] = p[1]*xB_tmp[13] - 1.0*p[1]*xB_tmp[12] - 1.0*p[2]*xB_tmp[7] - 1.0*p[1]*xB_tmp[16] - (1.0*xB_tmp[26]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*p[1]*xB_tmp[5])/k[0] - (1.0*p[2]*xB_tmp[6])/k[0] - (1.0*p[1]*xB_tmp[8])/k[0] - (1.0*p[1]*xB_tmp[18])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*xB_tmp[9])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*xB_tmp[10])/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[1]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*xB_tmp[19]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
xBdot_tmp[14] = (xB_tmp[1]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*xB_tmp[4]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[21]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[2]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*xB_tmp[5])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*xB_tmp[17])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*xB_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[22])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[25])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*xB_tmp[19])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*xB_tmp[26])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*xB_tmp[0])/(pow(k[0],2)) + k[0]*p[0]*x_tmp[16]*xB_tmp[14] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[15] + k[0]*p[0]*x_tmp[16]*xB_tmp[16];
xBdot_tmp[15] = p[1]*xB_tmp[15] - 1.0*p[1]*xB_tmp[14] - 1.0*p[2]*xB_tmp[7] - 1.0*p[1]*xB_tmp[16] - (1.0*xB_tmp[23]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*p[1]*xB_tmp[2])/k[0] - (1.0*p[1]*xB_tmp[5])/k[0] - (1.0*p[2]*xB_tmp[6])/k[0] - (1.0*p[1]*xB_tmp[21])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*xB_tmp[3])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*xB_tmp[4])/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[17]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*xB_tmp[19]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
xBdot_tmp[16] = (xB_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - (1.0*xB_tmp[4]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) - (1.0*xB_tmp[10]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) - (1.0*xB_tmp[26]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) - (1.0*xB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (xB_tmp[18]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (xB_tmp[21]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*xB_tmp[2])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*xB_tmp[8])/(pow(k[0],2)) + (xB_tmp[1]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) + (xB_tmp[17]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + (xB_tmp[19]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[25]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[11])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*xB_tmp[22])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*xB_tmp[24])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*xB_tmp[0])/(pow(k[0],2)) - (1.0*xB_tmp[5]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[12]*xB_tmp[12] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[13] + k[0]*p[0]*x_tmp[14]*xB_tmp[14] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[15];
xBdot_tmp[17] = (xB_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[11] - 1.0*p[1]*xB_tmp[18] - 1.0*p[1]*xB_tmp[0] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[19] + k[0]*p[0]*x_tmp[16]*xB_tmp[23];
xBdot_tmp[18] = (xB_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[10] + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[9])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[8])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*xB_tmp[5])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*xB_tmp[26])/(pow(k[0],2)) + k[0]*p[0]*xB_tmp[12] - 1.0*k[0]*p[0]*xB_tmp[13] + k[0]*p[0]*xB_tmp[16] + k[0]*p[0]*x_tmp[14]*xB_tmp[0] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[17];
xBdot_tmp[19] = 2.0*p[1]*xB_tmp[19] - 1.0*p[1]*xB_tmp[17] - 1.0*p[1]*xB_tmp[1] - 1.0*p[2]*xB_tmp[20] - 1.0*p[1]*xB_tmp[23] - 1.0*p[2]*xB_tmp[24] - 1.0*p[1]*xB_tmp[26];
xBdot_tmp[20] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[20])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[11] - 1.0*p[1]*xB_tmp[25] - 2.0*p[2]*xB_tmp[6];
xBdot_tmp[21] = (xB_tmp[21]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[4] + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*xB_tmp[3])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*xB_tmp[2])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*xB_tmp[5])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*xB_tmp[23])/(pow(k[0],2)) + k[0]*p[0]*xB_tmp[14] - 1.0*k[0]*p[0]*xB_tmp[15] + k[0]*p[0]*xB_tmp[16] + k[0]*p[0]*x_tmp[12]*xB_tmp[0] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[1];
xBdot_tmp[22] = (xB_tmp[22]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[24] + k[0]*p[0]*x_tmp[16]*xB_tmp[25];
xBdot_tmp[23] = (xB_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[21] - 1.0*p[2]*xB_tmp[25] - 2.0*p[1]*xB_tmp[5] + k[0]*p[0]*x_tmp[14]*xB_tmp[3] - 2.0*k[0]*p[0]*x_tmp[14]*xB_tmp[4] + k[0]*p[0]*x_tmp[12]*xB_tmp[17] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[19];
xBdot_tmp[24] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[24])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[22] - 1.0*p[1]*xB_tmp[25] - 2.0*p[2]*xB_tmp[6];
xBdot_tmp[25] = (xB_tmp[25]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[12]*xB_tmp[11] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[20] + k[0]*p[0]*x_tmp[14]*xB_tmp[22] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[24];
xBdot_tmp[26] = (xB_tmp[26]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[18] - 1.0*p[2]*xB_tmp[25] - 2.0*p[1]*xB_tmp[5] + k[0]*p[0]*x_tmp[14]*xB_tmp[1] + k[0]*p[0]*x_tmp[12]*xB_tmp[9] - 2.0*k[0]*p[0]*x_tmp[12]*xB_tmp[10] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[19];
xBdot_tmp[27] = 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[27] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[28] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[44] + k[0]*p[0]*x_tmp[16]*xB_tmp[45] + k[0]*p[0]*x_tmp[16]*xB_tmp[48];
xBdot_tmp[28] = (xB_tmp[28]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[48] - 1.0*p[2]*xB_tmp[49] - 1.0*p[1]*xB_tmp[27] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[46] + k[0]*p[0]*x_tmp[16]*xB_tmp[53];
xBdot_tmp[29] = 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[29] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[30] + k[0]*p[0]*x_tmp[16]*xB_tmp[48];
xBdot_tmp[30] = (xB_tmp[30]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[48] - 1.0*p[2]*xB_tmp[49] - 2.0*p[1]*xB_tmp[29] - 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[31] + k[0]*p[0]*x_tmp[16]*xB_tmp[50];
xBdot_tmp[31] = 2.0*p[1]*xB_tmp[31] - 1.0*p[1]*xB_tmp[30] - 1.0*p[1]*xB_tmp[50] - 1.0*p[2]*xB_tmp[51];
xBdot_tmp[32] = (xB_tmp[32]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[12]*xB_tmp[45] + k[0]*p[0]*x_tmp[14]*xB_tmp[48] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[50] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[53];
xBdot_tmp[33] = 2.0*p[3]*xB_tmp[33];
xBdot_tmp[34] = p[3]*xB_tmp[34] - (1.0*p[3]*xB_tmp[33])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*xB_tmp[52])/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[38]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*xB_tmp[47]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[49]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*xB_tmp[51]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
xBdot_tmp[35] = 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[35] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[36] + k[0]*p[0]*x_tmp[16]*xB_tmp[45];
xBdot_tmp[36] = (xB_tmp[36]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[38] - 1.0*p[1]*xB_tmp[45] - 2.0*p[1]*xB_tmp[35] - 2.0*k[0]*p[0]*x_tmp[16]*xB_tmp[37] + k[0]*p[0]*x_tmp[16]*xB_tmp[53];
xBdot_tmp[37] = 2.0*p[1]*xB_tmp[37] - 1.0*p[1]*xB_tmp[36] - 1.0*p[2]*xB_tmp[47] - 1.0*p[1]*xB_tmp[53];
xBdot_tmp[38] = (xB_tmp[38]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[47] + k[0]*p[0]*x_tmp[16]*xB_tmp[52];
xBdot_tmp[39] = (xB_tmp[44]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[37]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*xB_tmp[45]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*xB_tmp[36]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*xB_tmp[53]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*xB_tmp[35]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*xB_tmp[32])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*xB_tmp[28])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[38])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[47])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*xB_tmp[48])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[52])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*xB_tmp[46])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*xB_tmp[50])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*xB_tmp[27])/(pow(k[0],2)) + k[0]*p[0]*x_tmp[16]*xB_tmp[39] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[40] + k[0]*p[0]*x_tmp[16]*xB_tmp[43];
xBdot_tmp[40] = p[1]*xB_tmp[40] - 1.0*p[1]*xB_tmp[39] - 1.0*p[2]*xB_tmp[34] - 1.0*p[1]*xB_tmp[43] - (1.0*xB_tmp[53]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*p[1]*xB_tmp[32])/k[0] - (1.0*p[2]*xB_tmp[33])/k[0] - (1.0*p[1]*xB_tmp[35])/k[0] - (1.0*p[1]*xB_tmp[45])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*xB_tmp[36])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*xB_tmp[37])/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[28]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*xB_tmp[46]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
xBdot_tmp[41] = (xB_tmp[28]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*xB_tmp[31]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[48]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*xB_tmp[30]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[50]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[29]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*xB_tmp[32])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*xB_tmp[44])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*xB_tmp[45])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[49])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[51])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*xB_tmp[52])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*xB_tmp[46])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*xB_tmp[53])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*xB_tmp[27])/(pow(k[0],2)) + k[0]*p[0]*x_tmp[16]*xB_tmp[41] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[42] + k[0]*p[0]*x_tmp[16]*xB_tmp[43];
xBdot_tmp[42] = p[1]*xB_tmp[42] - 1.0*p[1]*xB_tmp[41] - 1.0*p[2]*xB_tmp[34] - 1.0*p[1]*xB_tmp[43] - (1.0*xB_tmp[50]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*p[1]*xB_tmp[29])/k[0] - (1.0*p[1]*xB_tmp[32])/k[0] - (1.0*p[2]*xB_tmp[33])/k[0] - (1.0*p[1]*xB_tmp[48])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*xB_tmp[30])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*xB_tmp[31])/(pow(k[0],2)) - (1.0*p[0]*xB_tmp[44]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*xB_tmp[46]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
xBdot_tmp[43] = (xB_tmp[43]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - (1.0*xB_tmp[31]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) - (1.0*xB_tmp[37]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) - (1.0*xB_tmp[53]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*xB_tmp[50]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*xB_tmp[30]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) - (1.0*xB_tmp[36]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (xB_tmp[45]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (xB_tmp[48]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*xB_tmp[29])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*xB_tmp[35])/(pow(k[0],2)) + (xB_tmp[28]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) + (xB_tmp[44]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + (xB_tmp[46]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[52]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[38])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*xB_tmp[47])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*xB_tmp[49])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*xB_tmp[51])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*xB_tmp[27])/(pow(k[0],2)) - (1.0*xB_tmp[32]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[12]*xB_tmp[39] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[40] + k[0]*p[0]*x_tmp[14]*xB_tmp[41] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[42];
xBdot_tmp[44] = (xB_tmp[44]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[38] - 1.0*p[1]*xB_tmp[45] - 1.0*p[1]*xB_tmp[27] - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[46] + k[0]*p[0]*x_tmp[16]*xB_tmp[50];
xBdot_tmp[45] = (xB_tmp[45]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[37] + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[36])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*xB_tmp[35])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*xB_tmp[32])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*xB_tmp[53])/(pow(k[0],2)) + k[0]*p[0]*xB_tmp[39] - 1.0*k[0]*p[0]*xB_tmp[40] + k[0]*p[0]*xB_tmp[43] + k[0]*p[0]*x_tmp[14]*xB_tmp[27] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[44];
xBdot_tmp[46] = 2.0*p[1]*xB_tmp[46] - 1.0*p[1]*xB_tmp[44] - 1.0*p[1]*xB_tmp[28] - 1.0*p[2]*xB_tmp[47] - 1.0*p[1]*xB_tmp[50] - 1.0*p[2]*xB_tmp[51] - 1.0*p[1]*xB_tmp[53];
xBdot_tmp[47] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[47])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[38] - 1.0*p[1]*xB_tmp[52] - 2.0*p[2]*xB_tmp[33];
xBdot_tmp[48] = (xB_tmp[48]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[31] + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*xB_tmp[30])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*xB_tmp[29])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*xB_tmp[32])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*xB_tmp[50])/(pow(k[0],2)) + k[0]*p[0]*xB_tmp[41] - 1.0*k[0]*p[0]*xB_tmp[42] + k[0]*p[0]*xB_tmp[43] + k[0]*p[0]*x_tmp[12]*xB_tmp[27] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[28];
xBdot_tmp[49] = (xB_tmp[49]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*x_tmp[16]*xB_tmp[51] + k[0]*p[0]*x_tmp[16]*xB_tmp[52];
xBdot_tmp[50] = (xB_tmp[50]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[48] - 1.0*p[2]*xB_tmp[52] - 2.0*p[1]*xB_tmp[32] + k[0]*p[0]*x_tmp[14]*xB_tmp[30] - 2.0*k[0]*p[0]*x_tmp[14]*xB_tmp[31] + k[0]*p[0]*x_tmp[12]*xB_tmp[44] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[46];
xBdot_tmp[51] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[51])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[49] - 1.0*p[1]*xB_tmp[52] - 2.0*p[2]*xB_tmp[33];
xBdot_tmp[52] = (xB_tmp[52]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + k[0]*p[0]*x_tmp[12]*xB_tmp[38] - 1.0*k[0]*p[0]*x_tmp[12]*xB_tmp[47] + k[0]*p[0]*x_tmp[14]*xB_tmp[49] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[51];
xBdot_tmp[53] = (xB_tmp[53]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[45] - 1.0*p[2]*xB_tmp[52] - 2.0*p[1]*xB_tmp[32] + k[0]*p[0]*x_tmp[14]*xB_tmp[28] + k[0]*p[0]*x_tmp[12]*xB_tmp[36] - 2.0*k[0]*p[0]*x_tmp[12]*xB_tmp[37] - 1.0*k[0]*p[0]*x_tmp[14]*xB_tmp[46];

  for (ixB=0; ixB<54; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_enhancer_3(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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
qBdot_tmp[0+ip*ny] = (xB_tmp[9]*((pow(k[0],2))*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[16] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[12]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[3]*((pow(k[0],2))*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[2]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[14]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[2]*((pow(k[0],2))*x_tmp[21] - 2.0*(pow(k[0],3))*x_tmp[2]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[14]*x_tmp[21]))/(pow(k[0],2)) - (1.0*xB_tmp[8]*((pow(k[0],2))*x_tmp[18] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[16] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[12]*x_tmp[18]))/(pow(k[0],2)) - (1.0*xB_tmp[0]*((pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[21] - 2.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 2.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 6.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[12]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]))/k[0] - (1.0*xB_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]))/k[0] + (xB_tmp[14]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]))/k[0] - (1.0*xB_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]))/k[0] + (xB_tmp[25]*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[18] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[19]*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[18]*((pow(k[0],3))*x_tmp[5]*x_tmp[12] - 1.0*(pow(k[0],2))*x_tmp[18] + (pow(k[0],3))*x_tmp[8]*x_tmp[16] - 1.0*(pow(k[0],2))*x_tmp[12]*x_tmp[16] + (pow(k[0],3))*x_tmp[12]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[12]*x_tmp[21] + (pow(k[0],3))*x_tmp[16]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[21]*((pow(k[0],3))*x_tmp[2]*x_tmp[16] - 1.0*(pow(k[0],2))*x_tmp[21] + (pow(k[0],3))*x_tmp[5]*x_tmp[14] - 1.0*(pow(k[0],2))*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[14]*x_tmp[21] + (pow(k[0],3))*x_tmp[16]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[11]*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[18] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[20]*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[18] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[22]*(k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[24]*(k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[5]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[21] - 2.0*(pow(k[0],3))*x_tmp[5]*x_tmp[12] - 2.0*(pow(k[0],3))*x_tmp[5]*x_tmp[14] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[16]*x_tmp[18] - 2.0*(pow(k[0],3))*x_tmp[16]*x_tmp[21]))/(pow(k[0],2)) + (xB_tmp[1]*((pow(k[0],3))*x_tmp[14]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[21] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[17]*((pow(k[0],3))*x_tmp[12]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[18] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[10]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[13]*x_tmp[18] + 2.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 6.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[4]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[15]*x_tmp[21] + 2.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + 2.0*k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 2.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) - 6.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[26]*((pow(k[0],2))*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[5]*x_tmp[12] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[16]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[23]*((pow(k[0],2))*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[5]*x_tmp[14] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[16]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[16]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]))/k[0];
qBdot_tmp[1+ip*ny] = (xB_tmp[36]*((pow(k[0],2))*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[16] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[12]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[30]*((pow(k[0],2))*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[2]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[14]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[29]*((pow(k[0],2))*x_tmp[21] - 2.0*(pow(k[0],3))*x_tmp[2]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[14]*x_tmp[21]))/(pow(k[0],2)) - (1.0*xB_tmp[35]*((pow(k[0],2))*x_tmp[18] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[16] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[12]*x_tmp[18]))/(pow(k[0],2)) - (1.0*xB_tmp[27]*((pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[21] - 2.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 2.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 6.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[39]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]))/k[0] - (1.0*xB_tmp[40]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]))/k[0] + (xB_tmp[41]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]))/k[0] - (1.0*xB_tmp[42]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]))/k[0] + (xB_tmp[52]*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[18] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[46]*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[45]*((pow(k[0],3))*x_tmp[5]*x_tmp[12] - 1.0*(pow(k[0],2))*x_tmp[18] + (pow(k[0],3))*x_tmp[8]*x_tmp[16] - 1.0*(pow(k[0],2))*x_tmp[12]*x_tmp[16] + (pow(k[0],3))*x_tmp[12]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[12]*x_tmp[21] + (pow(k[0],3))*x_tmp[16]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[48]*((pow(k[0],3))*x_tmp[2]*x_tmp[16] - 1.0*(pow(k[0],2))*x_tmp[21] + (pow(k[0],3))*x_tmp[5]*x_tmp[14] - 1.0*(pow(k[0],2))*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[14]*x_tmp[21] + (pow(k[0],3))*x_tmp[16]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[38]*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[18] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[47]*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[18] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[49]*(k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[51]*(k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[32]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[21] - 2.0*(pow(k[0],3))*x_tmp[5]*x_tmp[12] - 2.0*(pow(k[0],3))*x_tmp[5]*x_tmp[14] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[16]*x_tmp[18] - 2.0*(pow(k[0],3))*x_tmp[16]*x_tmp[21]))/(pow(k[0],2)) + (xB_tmp[28]*((pow(k[0],3))*x_tmp[14]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[21] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[44]*((pow(k[0],3))*x_tmp[12]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[18] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[37]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[13]*x_tmp[18] + 2.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 6.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[31]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[15]*x_tmp[21] + 2.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + 2.0*k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 2.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) - 6.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[53]*((pow(k[0],2))*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[5]*x_tmp[12] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[16]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[50]*((pow(k[0],2))*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[5]*x_tmp[14] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[16]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) + (xB_tmp[43]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]))/k[0];

  } break;

  case 1: {
qBdot_tmp[0+ip*ny] = x_tmp[13]*xB_tmp[13] - 1.0*x_tmp[13]*xB_tmp[12] - 1.0*x_tmp[15]*xB_tmp[14] + x_tmp[15]*xB_tmp[15] - 1.0*x_tmp[20]*xB_tmp[11] + 2.0*x_tmp[19]*xB_tmp[19] + x_tmp[20]*xB_tmp[20] - 1.0*x_tmp[24]*xB_tmp[22] + x_tmp[24]*xB_tmp[24] - (1.0*xB_tmp[5]*(k[0]*x_tmp[13] + k[0]*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[23] + 2.0*(pow(k[0],2))*x_tmp[26]))/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[2])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[15] - 2.0*(pow(k[0],2))*x_tmp[4])*xB_tmp[4])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[9])*xB_tmp[8])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[10])*xB_tmp[10])/(pow(k[0],2)) + (xB_tmp[23]*(k[0]*x_tmp[15] - 1.0*(pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[23]))/(pow(k[0],2)) + (xB_tmp[26]*(k[0]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[10] - 1.0*(pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[26]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[17])*xB_tmp[0])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[20] + (pow(k[0],2))*x_tmp[24])*xB_tmp[25])/(pow(k[0],2)) + (((pow(k[0],2))*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[19])*xB_tmp[1])/(pow(k[0],2)) + (((pow(k[0],2))*x_tmp[17] - 1.0*(pow(k[0],2))*x_tmp[19])*xB_tmp[17])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[13] + k[0]*x_tmp[15])*xB_tmp[16])/k[0] + (xB_tmp[3]*(k[0]*x_tmp[15] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[4]))/(pow(k[0],2)) + (xB_tmp[9]*(k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[9] - 1.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[21]*(k[0]*x_tmp[15] + (pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[18]*(k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[26]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = x_tmp[13]*xB_tmp[40] - 1.0*x_tmp[13]*xB_tmp[39] - 1.0*x_tmp[15]*xB_tmp[41] + x_tmp[15]*xB_tmp[42] - 1.0*x_tmp[20]*xB_tmp[38] + 2.0*x_tmp[19]*xB_tmp[46] + x_tmp[20]*xB_tmp[47] - 1.0*x_tmp[24]*xB_tmp[49] + x_tmp[24]*xB_tmp[51] - (1.0*xB_tmp[32]*(k[0]*x_tmp[13] + k[0]*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[23] + 2.0*(pow(k[0],2))*x_tmp[26]))/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[29])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[15] - 2.0*(pow(k[0],2))*x_tmp[4])*xB_tmp[31])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[9])*xB_tmp[35])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[10])*xB_tmp[37])/(pow(k[0],2)) + (xB_tmp[50]*(k[0]*x_tmp[15] - 1.0*(pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[23]))/(pow(k[0],2)) + (xB_tmp[53]*(k[0]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[10] - 1.0*(pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[26]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[17])*xB_tmp[27])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[20] + (pow(k[0],2))*x_tmp[24])*xB_tmp[52])/(pow(k[0],2)) + (((pow(k[0],2))*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[19])*xB_tmp[28])/(pow(k[0],2)) + (((pow(k[0],2))*x_tmp[17] - 1.0*(pow(k[0],2))*x_tmp[19])*xB_tmp[44])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[13] + k[0]*x_tmp[15])*xB_tmp[43])/k[0] + (xB_tmp[30]*(k[0]*x_tmp[15] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[4]))/(pow(k[0],2)) + (xB_tmp[36]*(k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[9] - 1.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[48]*(k[0]*x_tmp[15] + (pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[23]))/(pow(k[0],2)) - (1.0*xB_tmp[45]*(k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[26]))/(pow(k[0],2));

  } break;

  case 2: {
qBdot_tmp[0+ip*ny] = - (1.0*xB_tmp[6]*(k[0]*x_tmp[13] + k[0]*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[20] + 2.0*(pow(k[0],2))*x_tmp[24]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[3])*xB_tmp[22])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[17])*xB_tmp[11])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[19])*xB_tmp[24])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[19])*xB_tmp[20])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[26])*xB_tmp[25])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[13] + k[0]*x_tmp[15])*xB_tmp[7])/k[0];
qBdot_tmp[1+ip*ny] = - (1.0*xB_tmp[33]*(k[0]*x_tmp[13] + k[0]*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[20] + 2.0*(pow(k[0],2))*x_tmp[24]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[3])*xB_tmp[49])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[17])*xB_tmp[38])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[19])*xB_tmp[51])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[19])*xB_tmp[47])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[26])*xB_tmp[52])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[13] + k[0]*x_tmp[15])*xB_tmp[34])/k[0];

  } break;

  case 3: {
qBdot_tmp[0+ip*ny] = x_tmp[7]*xB_tmp[7] + x_tmp[11]*xB_tmp[11] + x_tmp[20]*xB_tmp[20] + x_tmp[22]*xB_tmp[22] + x_tmp[24]*xB_tmp[24] + x_tmp[25]*xB_tmp[25] - (1.0*(k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[6])*xB_tmp[6])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = x_tmp[7]*xB_tmp[34] + x_tmp[11]*xB_tmp[38] + x_tmp[20]*xB_tmp[47] + x_tmp[22]*xB_tmp[49] + x_tmp[24]*xB_tmp[51] + x_tmp[25]*xB_tmp[52] - (1.0*(k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[6])*xB_tmp[33])/(pow(k[0],2));

  } break;

  case 4: {
qBdot_tmp[0+ip*ny] = - (144.0*xB_tmp[5])/(pow(k[0],2)) - (12.0*xB_tmp[16])/k[0];
qBdot_tmp[1+ip*ny] = - (144.0*xB_tmp[32])/(pow(k[0],2)) - (12.0*xB_tmp[43])/k[0];

  } break;

  case 5: {
qBdot_tmp[0+ip*ny] = x_tmp[16]*xB_tmp[16] + x_tmp[18]*xB_tmp[18] + x_tmp[21]*xB_tmp[21] + x_tmp[23]*xB_tmp[23] + x_tmp[25]*xB_tmp[25] + x_tmp[26]*xB_tmp[26] - (1.0*(k[0]*x_tmp[16] - 2.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[5])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = x_tmp[16]*xB_tmp[43] + x_tmp[18]*xB_tmp[45] + x_tmp[21]*xB_tmp[48] + x_tmp[23]*xB_tmp[50] + x_tmp[25]*xB_tmp[52] + x_tmp[26]*xB_tmp[53] - (1.0*(k[0]*x_tmp[16] - 2.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[32])/(pow(k[0],2));

  } break;

  }
  }

  for (iyp=0; iyp<2*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_enhancer_3(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*27);
x0_tmp[0] = (k[9]*k[36])/(pow(k[0],2));
x0_tmp[1] = (k[14]*k[41])/(pow(k[0],2));
x0_tmp[2] = (k[18]*k[45])/(pow(k[0],2));
x0_tmp[3] = (k[19]*k[46])/(pow(k[0],2));
x0_tmp[4] = (k[22]*k[49])/(pow(k[0],2));
x0_tmp[5] = (k[25]*k[52])/(pow(k[0],2));
x0_tmp[6] = (k[27]*k[54])/(pow(k[0],2));
x0_tmp[7] = (k[6]*k[33])/k[0];
x0_tmp[8] = (k[7]*k[34])/(pow(k[0],2));
x0_tmp[9] = (k[8]*k[35])/(pow(k[0],2));
x0_tmp[10] = (k[13]*k[40])/(pow(k[0],2));
x0_tmp[11] = (k[12]*k[39])/(pow(k[0],2));
x0_tmp[12] = (k[1]*k[28] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[13] = (k[2]*k[29])/k[0];
x0_tmp[14] = (k[3]*k[30] - 1.0*k[3] + 1.0)/k[0];
x0_tmp[15] = (k[4]*k[31])/k[0];
x0_tmp[16] = (k[5]*k[32])/k[0];
x0_tmp[17] = (k[10]*k[37])/(pow(k[0],2));
x0_tmp[18] = (k[11]*k[38])/(pow(k[0],2));
x0_tmp[19] = (k[15]*k[42])/(pow(k[0],2));
x0_tmp[20] = (k[17]*k[44])/(pow(k[0],2));
x0_tmp[21] = (k[20]*k[47])/(pow(k[0],2));
x0_tmp[22] = (k[21]*k[48])/(pow(k[0],2));
x0_tmp[23] = (k[23]*k[50])/(pow(k[0],2));
x0_tmp[24] = (k[24]*k[51])/(pow(k[0],2));
x0_tmp[25] = (k[26]*k[53])/(pow(k[0],2));
x0_tmp[26] = (k[16]*k[43])/(pow(k[0],2));
  
  
  return;
}


 int Jv_enhancer_3(N_Vector v, N_Vector Jv, realtype t,
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
  memset(Jv_tmp,0,sizeof(double)*27);
Jv_tmp[0] = p[1]*v_tmp[1] + p[1]*v_tmp[17] + ((2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*v_tmp[16])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*v_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*v_tmp[12])/(pow(k[0],2)) - 2.0*k[0]*p[0]*v_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*v_tmp[21]*x_tmp[12];
Jv_tmp[1] = p[1]*v_tmp[19] - (1.0*v_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*v_tmp[1]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*v_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*v_tmp[12])/(pow(k[0],2)) + (p[0]*v_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*v_tmp[0]*x_tmp[16] + k[0]*p[0]*v_tmp[21]*x_tmp[12] - 1.0*k[0]*p[0]*v_tmp[26]*x_tmp[14];
Jv_tmp[2] = 2.0*p[1]*v_tmp[3] + (p[1]*v_tmp[15])/k[0] + (v_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*v_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*v_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*v_tmp[2]*x_tmp[16];
Jv_tmp[3] = p[1]*v_tmp[4] - (1.0*v_tmp[3]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (v_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (v_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*v_tmp[15])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*v_tmp[21])/(pow(k[0],2)) + k[0]*p[0]*v_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[23]*x_tmp[14];
Jv_tmp[4] = p[0]*v_tmp[21] - 2.0*p[1]*v_tmp[4] + (v_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (v_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*v_tmp[15])/(pow(k[0],2)) + 2.0*k[0]*p[0]*v_tmp[3]*x_tmp[16] + 2.0*k[0]*p[0]*v_tmp[23]*x_tmp[14];
Jv_tmp[5] = 2.0*p[1]*v_tmp[23] + 2.0*p[1]*v_tmp[26] + (p[1]*v_tmp[13])/k[0] + (p[1]*v_tmp[15])/k[0] - (1.0*v_tmp[5]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*v_tmp[18])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*v_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*v_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*v_tmp[14])/(pow(k[0],2)) + (v_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
Jv_tmp[6] = 2.0*p[2]*v_tmp[20] - 2.0*p[3]*v_tmp[6] + 2.0*p[2]*v_tmp[24] + (p[3]*v_tmp[7])/k[0] + (p[2]*v_tmp[13])/k[0] + (p[2]*v_tmp[15])/k[0];
Jv_tmp[7] = p[2]*v_tmp[13] - 1.0*p[3]*v_tmp[7] + p[2]*v_tmp[15];
Jv_tmp[8] = 2.0*p[1]*v_tmp[9] + (p[1]*v_tmp[13])/k[0] + (v_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*v_tmp[18])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*v_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*v_tmp[8]*x_tmp[16];
Jv_tmp[9] = p[1]*v_tmp[10] - (1.0*v_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (v_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (v_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*v_tmp[13])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*v_tmp[18])/(pow(k[0],2)) + k[0]*p[0]*v_tmp[8]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[26]*x_tmp[12];
Jv_tmp[10] = p[0]*v_tmp[18] - 2.0*p[1]*v_tmp[10] + (v_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (v_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*v_tmp[13])/(pow(k[0],2)) + 2.0*k[0]*p[0]*v_tmp[9]*x_tmp[16] + 2.0*k[0]*p[0]*v_tmp[26]*x_tmp[12];
Jv_tmp[11] = p[2]*v_tmp[9] + p[2]*v_tmp[17] + p[1]*v_tmp[20] - (1.0*v_tmp[11]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*v_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*v_tmp[12])/(pow(k[0],2)) + (p[0]*v_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[25]*x_tmp[12];
Jv_tmp[12] = p[1]*v_tmp[13] - 1.0*k[0]*p[0]*v_tmp[18] - 1.0*k[0]*p[0]*v_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[16]*x_tmp[12];
Jv_tmp[13] = k[0]*p[0]*v_tmp[18] - 1.0*p[1]*v_tmp[13] + k[0]*p[0]*v_tmp[12]*x_tmp[16] + k[0]*p[0]*v_tmp[16]*x_tmp[12];
Jv_tmp[14] = p[1]*v_tmp[15] - 1.0*k[0]*p[0]*v_tmp[21] - 1.0*k[0]*p[0]*v_tmp[14]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[16]*x_tmp[14];
Jv_tmp[15] = k[0]*p[0]*v_tmp[21] - 1.0*p[1]*v_tmp[15] + k[0]*p[0]*v_tmp[14]*x_tmp[16] + k[0]*p[0]*v_tmp[16]*x_tmp[14];
Jv_tmp[16] = p[1]*v_tmp[13] + p[1]*v_tmp[15] - (1.0*v_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - 1.0*k[0]*p[0]*v_tmp[18] - 1.0*k[0]*p[0]*v_tmp[21] - 1.0*k[0]*p[0]*v_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[14]*x_tmp[16];
Jv_tmp[17] = p[1]*v_tmp[19] - (1.0*v_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*v_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*v_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*v_tmp[14])/(pow(k[0],2)) + (p[0]*v_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*v_tmp[0]*x_tmp[16] + k[0]*p[0]*v_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*v_tmp[23]*x_tmp[12];
Jv_tmp[18] = p[1]*v_tmp[9] + p[1]*v_tmp[17] + p[1]*v_tmp[26] + (p[1]*v_tmp[13])/k[0] - (1.0*v_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (v_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*v_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*v_tmp[14])/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*v_tmp[8]*x_tmp[16];
Jv_tmp[19] = k[0]*p[0]*v_tmp[1]*x_tmp[16] - (1.0*v_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*v_tmp[12])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*v_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*v_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*v_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 2.0*p[1]*v_tmp[19] + k[0]*p[0]*v_tmp[17]*x_tmp[16] + k[0]*p[0]*v_tmp[23]*x_tmp[12] + k[0]*p[0]*v_tmp[26]*x_tmp[14];
Jv_tmp[20] = p[2]*v_tmp[10] + p[2]*v_tmp[19] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*v_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*v_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*v_tmp[12])/(pow(k[0],2)) - (1.0*p[0]*v_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*v_tmp[11]*x_tmp[16] + k[0]*p[0]*v_tmp[25]*x_tmp[12];
Jv_tmp[21] = p[1]*v_tmp[1] + p[1]*v_tmp[3] + p[1]*v_tmp[23] + (p[1]*v_tmp[15])/k[0] - (1.0*v_tmp[21]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (v_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*v_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*v_tmp[12])/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[5]*x_tmp[14];
Jv_tmp[22] = p[2]*v_tmp[1] + p[2]*v_tmp[3] + p[1]*v_tmp[24] - (1.0*v_tmp[22]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*v_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*v_tmp[14])/(pow(k[0],2)) + (p[0]*v_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[25]*x_tmp[14];
Jv_tmp[23] = p[1]*v_tmp[4] + p[1]*v_tmp[19] + (v_tmp[15]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (v_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (v_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*v_tmp[21])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*v_tmp[12])/(pow(k[0],2)) - (1.0*v_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[3]*x_tmp[16] + k[0]*p[0]*v_tmp[5]*x_tmp[14] - 1.0*k[0]*p[0]*v_tmp[17]*x_tmp[16];
Jv_tmp[24] = p[2]*v_tmp[4] + p[2]*v_tmp[19] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*v_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*v_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*v_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*v_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*v_tmp[22]*x_tmp[16] + k[0]*p[0]*v_tmp[25]*x_tmp[14];
Jv_tmp[25] = p[1]*v_tmp[20] + p[1]*v_tmp[24] + p[2]*v_tmp[23] + p[2]*v_tmp[26] + (v_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*v_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*v_tmp[14])/(pow(k[0],2)) - (1.0*v_tmp[25]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*v_tmp[7])/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[11]*x_tmp[16] - 1.0*k[0]*p[0]*v_tmp[22]*x_tmp[16];
Jv_tmp[26] = p[1]*v_tmp[10] + p[1]*v_tmp[19] + (v_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (v_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (v_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*v_tmp[18])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*v_tmp[14])/(pow(k[0],2)) - (1.0*v_tmp[26]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*v_tmp[1]*x_tmp[16] + k[0]*p[0]*v_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*v_tmp[9]*x_tmp[16];

  for (ix=0; ix<27; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_enhancer_3(N_Vector vB, N_Vector JvB, realtype t,
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
  memset(JvB_tmp,0,sizeof(double)*27);
JvB_tmp[0] = 2.0*k[0]*p[0]*vB_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*vB_tmp[1]*x_tmp[16] - 1.0*k[0]*p[0]*vB_tmp[17]*x_tmp[16] + k[0]*p[0]*vB_tmp[18]*x_tmp[16] + k[0]*p[0]*vB_tmp[21]*x_tmp[16];
JvB_tmp[1] = (vB_tmp[1]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[21] - 1.0*p[2]*vB_tmp[22] - 1.0*p[1]*vB_tmp[0] - 1.0*k[0]*p[0]*vB_tmp[19]*x_tmp[16] + k[0]*p[0]*vB_tmp[26]*x_tmp[16];
JvB_tmp[2] = 2.0*k[0]*p[0]*vB_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*vB_tmp[3]*x_tmp[16] + k[0]*p[0]*vB_tmp[21]*x_tmp[16];
JvB_tmp[3] = (vB_tmp[3]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[21] - 1.0*p[2]*vB_tmp[22] - 2.0*p[1]*vB_tmp[2] - 2.0*k[0]*p[0]*vB_tmp[4]*x_tmp[16] + k[0]*p[0]*vB_tmp[23]*x_tmp[16];
JvB_tmp[4] = 2.0*p[1]*vB_tmp[4] - 1.0*p[1]*vB_tmp[3] - 1.0*p[1]*vB_tmp[23] - 1.0*p[2]*vB_tmp[24];
JvB_tmp[5] = (vB_tmp[5]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[18]*x_tmp[12] + k[0]*p[0]*vB_tmp[21]*x_tmp[14] - 1.0*k[0]*p[0]*vB_tmp[23]*x_tmp[14] - 1.0*k[0]*p[0]*vB_tmp[26]*x_tmp[12];
JvB_tmp[6] = 2.0*p[3]*vB_tmp[6];
JvB_tmp[7] = p[3]*vB_tmp[7] - (1.0*p[3]*vB_tmp[6])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*vB_tmp[25])/(pow(k[0],2)) - (1.0*p[0]*vB_tmp[11]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*vB_tmp[20]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*vB_tmp[22]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*vB_tmp[24]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
JvB_tmp[8] = 2.0*k[0]*p[0]*vB_tmp[8]*x_tmp[16] - 1.0*k[0]*p[0]*vB_tmp[9]*x_tmp[16] + k[0]*p[0]*vB_tmp[18]*x_tmp[16];
JvB_tmp[9] = (vB_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[2]*vB_tmp[11] - 1.0*p[1]*vB_tmp[18] - 2.0*p[1]*vB_tmp[8] - 2.0*k[0]*p[0]*vB_tmp[10]*x_tmp[16] + k[0]*p[0]*vB_tmp[26]*x_tmp[16];
JvB_tmp[10] = 2.0*p[1]*vB_tmp[10] - 1.0*p[1]*vB_tmp[9] - 1.0*p[2]*vB_tmp[20] - 1.0*p[1]*vB_tmp[26];
JvB_tmp[11] = (vB_tmp[11]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*vB_tmp[20]*x_tmp[16] + k[0]*p[0]*vB_tmp[25]*x_tmp[16];
JvB_tmp[12] = (vB_tmp[17]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*vB_tmp[10]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*vB_tmp[18]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*vB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*vB_tmp[26]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*vB_tmp[8]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*vB_tmp[5])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*vB_tmp[1])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*vB_tmp[11])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*vB_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*vB_tmp[21])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*vB_tmp[25])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*vB_tmp[19])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*vB_tmp[23])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*vB_tmp[0])/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*vB_tmp[13]*x_tmp[16] + k[0]*p[0]*vB_tmp[16]*x_tmp[16];
JvB_tmp[13] = p[1]*vB_tmp[13] - 1.0*p[1]*vB_tmp[12] - 1.0*p[2]*vB_tmp[7] - 1.0*p[1]*vB_tmp[16] - (1.0*vB_tmp[26]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*p[1]*vB_tmp[5])/k[0] - (1.0*p[2]*vB_tmp[6])/k[0] - (1.0*p[1]*vB_tmp[8])/k[0] - (1.0*p[1]*vB_tmp[18])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*vB_tmp[9])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*vB_tmp[10])/(pow(k[0],2)) - (1.0*p[0]*vB_tmp[1]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*vB_tmp[19]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
JvB_tmp[14] = (vB_tmp[1]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*vB_tmp[4]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*vB_tmp[21]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*vB_tmp[3]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*vB_tmp[23]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*vB_tmp[2]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*vB_tmp[5])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*vB_tmp[17])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*vB_tmp[18])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*vB_tmp[22])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*vB_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*vB_tmp[25])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*vB_tmp[19])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*vB_tmp[26])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*vB_tmp[0])/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[14]*x_tmp[16] - 1.0*k[0]*p[0]*vB_tmp[15]*x_tmp[16] + k[0]*p[0]*vB_tmp[16]*x_tmp[16];
JvB_tmp[15] = p[1]*vB_tmp[15] - 1.0*p[1]*vB_tmp[14] - 1.0*p[2]*vB_tmp[7] - 1.0*p[1]*vB_tmp[16] - (1.0*vB_tmp[23]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) - (1.0*p[1]*vB_tmp[2])/k[0] - (1.0*p[1]*vB_tmp[5])/k[0] - (1.0*p[2]*vB_tmp[6])/k[0] - (1.0*p[1]*vB_tmp[21])/k[0] - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*vB_tmp[3])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*vB_tmp[4])/(pow(k[0],2)) - (1.0*p[0]*vB_tmp[17]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + (p[0]*vB_tmp[19]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
JvB_tmp[16] = (vB_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - (1.0*vB_tmp[4]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) - (1.0*vB_tmp[10]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) - (1.0*vB_tmp[26]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*vB_tmp[23]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*vB_tmp[3]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) - (1.0*vB_tmp[9]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (vB_tmp[18]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (vB_tmp[21]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*vB_tmp[2])/(pow(k[0],2)) + ((2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*vB_tmp[8])/(pow(k[0],2)) + (vB_tmp[1]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) + (vB_tmp[17]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) + (vB_tmp[19]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*vB_tmp[25]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*vB_tmp[11])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*vB_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*vB_tmp[22])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*vB_tmp[24])/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*vB_tmp[0])/(pow(k[0],2)) - (1.0*vB_tmp[5]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[12]*x_tmp[12] - 1.0*k[0]*p[0]*vB_tmp[13]*x_tmp[12] + k[0]*p[0]*vB_tmp[14]*x_tmp[14] - 1.0*k[0]*p[0]*vB_tmp[15]*x_tmp[14];
JvB_tmp[17] = (vB_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[2]*vB_tmp[11] - 1.0*p[1]*vB_tmp[18] - 1.0*p[1]*vB_tmp[0] - 1.0*k[0]*p[0]*vB_tmp[19]*x_tmp[16] + k[0]*p[0]*vB_tmp[23]*x_tmp[16];
JvB_tmp[18] = (vB_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[10] + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*vB_tmp[9])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*vB_tmp[8])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*vB_tmp[5])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*vB_tmp[26])/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[12] - 1.0*k[0]*p[0]*vB_tmp[13] + k[0]*p[0]*vB_tmp[16] + k[0]*p[0]*vB_tmp[0]*x_tmp[14] - 1.0*k[0]*p[0]*vB_tmp[17]*x_tmp[14];
JvB_tmp[19] = 2.0*p[1]*vB_tmp[19] - 1.0*p[1]*vB_tmp[17] - 1.0*p[1]*vB_tmp[1] - 1.0*p[2]*vB_tmp[20] - 1.0*p[1]*vB_tmp[23] - 1.0*p[2]*vB_tmp[24] - 1.0*p[1]*vB_tmp[26];
JvB_tmp[20] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*vB_tmp[20])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[11] - 1.0*p[1]*vB_tmp[25] - 2.0*p[2]*vB_tmp[6];
JvB_tmp[21] = (vB_tmp[21]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[4] + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*vB_tmp[3])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*vB_tmp[2])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*vB_tmp[5])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*vB_tmp[23])/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[14] - 1.0*k[0]*p[0]*vB_tmp[15] + k[0]*p[0]*vB_tmp[16] + k[0]*p[0]*vB_tmp[0]*x_tmp[12] - 1.0*k[0]*p[0]*vB_tmp[1]*x_tmp[12];
JvB_tmp[22] = (vB_tmp[22]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*vB_tmp[24]*x_tmp[16] + k[0]*p[0]*vB_tmp[25]*x_tmp[16];
JvB_tmp[23] = (vB_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[21] - 1.0*p[2]*vB_tmp[25] - 2.0*p[1]*vB_tmp[5] + k[0]*p[0]*vB_tmp[3]*x_tmp[14] - 2.0*k[0]*p[0]*vB_tmp[4]*x_tmp[14] + k[0]*p[0]*vB_tmp[17]*x_tmp[12] - 1.0*k[0]*p[0]*vB_tmp[19]*x_tmp[12];
JvB_tmp[24] = (((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*vB_tmp[24])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[22] - 1.0*p[1]*vB_tmp[25] - 2.0*p[2]*vB_tmp[6];
JvB_tmp[25] = (vB_tmp[25]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + k[0]*p[0]*vB_tmp[11]*x_tmp[12] - 1.0*k[0]*p[0]*vB_tmp[20]*x_tmp[12] + k[0]*p[0]*vB_tmp[22]*x_tmp[14] - 1.0*k[0]*p[0]*vB_tmp[24]*x_tmp[14];
JvB_tmp[26] = (vB_tmp[26]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[18] - 1.0*p[2]*vB_tmp[25] - 2.0*p[1]*vB_tmp[5] + k[0]*p[0]*vB_tmp[1]*x_tmp[14] + k[0]*p[0]*vB_tmp[9]*x_tmp[12] - 2.0*k[0]*p[0]*vB_tmp[10]*x_tmp[12] - 1.0*k[0]*p[0]*vB_tmp[19]*x_tmp[14];

  for (ix=0; ix<27; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_enhancer_3(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_enhancer_3(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_enhancer_3(long int N, realtype t, N_Vector x,
  	N_Vector xdot, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  memset(J->data,0,sizeof(double)*729);
J->data[0] = -2.0*k[0]*p[0]*x_tmp[16];
J->data[1] = k[0]*p[0]*x_tmp[16];
J->data[17] = k[0]*p[0]*x_tmp[16];
J->data[18] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[21] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[27] = p[1];
J->data[28] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[46] = k[0]*p[0]*x_tmp[16];
J->data[48] = p[1];
J->data[49] = p[2];
J->data[53] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[56] = -2.0*k[0]*p[0]*x_tmp[16];
J->data[57] = k[0]*p[0]*x_tmp[16];
J->data[75] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[83] = 2.0*p[1];
J->data[84] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[85] = 2.0*k[0]*p[0]*x_tmp[16];
J->data[102] = p[1];
J->data[103] = p[2];
J->data[104] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[111] = p[1];
J->data[112] = -2.0*p[1];
J->data[131] = p[1];
J->data[132] = p[2];
J->data[140] = -(1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[153] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[156] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[158] = k[0]*p[0]*x_tmp[14];
J->data[161] = k[0]*p[0]*x_tmp[12];
J->data[168] = -2.0*p[3];
J->data[195] = p[3]/k[0];
J->data[196] = -1.0*p[3];
J->data[200] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
J->data[209] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
J->data[211] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[213] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[214] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[224] = -2.0*k[0]*p[0]*x_tmp[16];
J->data[225] = k[0]*p[0]*x_tmp[16];
J->data[234] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[251] = 2.0*p[1];
J->data[252] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[253] = 2.0*k[0]*p[0]*x_tmp[16];
J->data[254] = p[2];
J->data[261] = p[1];
J->data[269] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[279] = p[1];
J->data[280] = -2.0*p[1];
J->data[290] = p[2];
J->data[296] = p[1];
J->data[308] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[317] = k[0]*p[0]*x_tmp[16];
J->data[322] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[324] = (2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[325] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
J->data[329] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[332] = ((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[333] = (p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
J->data[334] = ((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
J->data[335] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
J->data[336] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[337] = k[0]*p[0]*x_tmp[16];
J->data[340] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[341] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2));
J->data[342] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[343] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2));
J->data[344] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25]))/(pow(k[0],2));
J->data[345] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[347] = (p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
J->data[349] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
J->data[350] = (p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
J->data[352] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[356] = p[1]/k[0];
J->data[357] = p[2]/k[0];
J->data[358] = p[2];
J->data[359] = p[1]/k[0];
J->data[360] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[361] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
J->data[363] = p[1];
J->data[364] = -1.0*p[1];
J->data[367] = p[1];
J->data[369] = p[1]/k[0];
J->data[370] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[377] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[378] = (2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[379] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2));
J->data[380] = ((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[381] = (p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
J->data[382] = ((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
J->data[383] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[392] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[393] = k[0]*p[0]*x_tmp[16];
J->data[394] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[395] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2));
J->data[396] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[397] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2));
J->data[399] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[400] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
J->data[401] = (p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
J->data[402] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25]))/(pow(k[0],2));
J->data[403] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
J->data[404] = (p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
J->data[407] = p[1]/k[0];
J->data[408] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[409] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
J->data[410] = p[1]/k[0];
J->data[411] = p[2]/k[0];
J->data[412] = p[2];
J->data[419] = p[1];
J->data[420] = -1.0*p[1];
J->data[421] = p[1];
J->data[422] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
J->data[424] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
J->data[426] = p[1]/k[0];
J->data[428] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[432] = (2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])/(pow(k[0],2));
J->data[433] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2));
J->data[434] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[435] = (p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])/(pow(k[0],2));
J->data[436] = (2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14])/(pow(k[0],2));
J->data[437] = (k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[440] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2));
J->data[441] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[442] = (2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[443] = (p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[444] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[445] = k[0]*p[0]*x_tmp[12];
J->data[446] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[447] = k[0]*p[0]*x_tmp[14];
J->data[448] = -(1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0];
J->data[449] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2));
J->data[450] = -(1.0*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2));
J->data[451] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2));
J->data[452] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[453] = -(1.0*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
J->data[454] = (p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])/(pow(k[0],2));
J->data[455] = (p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[456] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2));
J->data[457] = (p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])/(pow(k[0],2));
J->data[458] = (p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[459] = p[1];
J->data[470] = p[2];
J->data[476] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[477] = p[1];
J->data[478] = k[0]*p[0]*x_tmp[16];
J->data[482] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[486] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[491] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
J->data[494] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[495] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]))/(pow(k[0],2));
J->data[496] = p[0];
J->data[498] = -1.0*k[0]*p[0];
J->data[499] = k[0]*p[0];
J->data[502] = -1.0*k[0]*p[0];
J->data[503] = k[0]*p[0]*x_tmp[14];
J->data[504] = -(1.0*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[512] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[514] = p[1];
J->data[530] = p[1];
J->data[532] = -2.0*p[1];
J->data[533] = p[2];
J->data[536] = p[1];
J->data[537] = p[2];
J->data[539] = p[1];
J->data[546] = 2.0*p[2];
J->data[551] = p[1];
J->data[560] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[565] = p[1];
J->data[567] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[568] = k[0]*p[0]*x_tmp[12];
J->data[569] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])/(pow(k[0],2));
J->data[570] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[571] = p[0];
J->data[572] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
J->data[581] = -1.0*k[0]*p[0];
J->data[582] = k[0]*p[0];
J->data[583] = -1.0*k[0]*p[0];
J->data[588] = -(1.0*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[590] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[616] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[618] = k[0]*p[0]*x_tmp[16];
J->data[619] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[624] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[625] = 2.0*k[0]*p[0]*x_tmp[14];
J->data[626] = 2.0*p[1];
J->data[638] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[640] = k[0]*p[0]*x_tmp[12];
J->data[642] = p[1];
J->data[644] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[646] = p[2];
J->data[654] = 2.0*p[2];
J->data[670] = p[1];
J->data[672] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[673] = p[1];
J->data[686] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[695] = k[0]*p[0]*x_tmp[12];
J->data[697] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[699] = k[0]*p[0]*x_tmp[14];
J->data[700] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[703] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[707] = 2.0*p[1];
J->data[711] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[712] = 2.0*k[0]*p[0]*x_tmp[12];
J->data[720] = p[1];
J->data[721] = k[0]*p[0]*x_tmp[14];
J->data[727] = p[2];
J->data[728] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));

  for (iJ=0; iJ<729; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_enhancer_3(realtype t, N_Vector x,
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
  J->rowvals[1] = 1;
  J->rowvals[2] = 17;
  J->rowvals[3] = 18;
  J->rowvals[4] = 21;
  J->rowvals[5] = 0;
  J->rowvals[6] = 1;
  J->rowvals[7] = 19;
  J->rowvals[8] = 21;
  J->rowvals[9] = 22;
  J->rowvals[10] = 26;
  J->rowvals[11] = 2;
  J->rowvals[12] = 3;
  J->rowvals[13] = 21;
  J->rowvals[14] = 2;
  J->rowvals[15] = 3;
  J->rowvals[16] = 4;
  J->rowvals[17] = 21;
  J->rowvals[18] = 22;
  J->rowvals[19] = 23;
  J->rowvals[20] = 3;
  J->rowvals[21] = 4;
  J->rowvals[22] = 23;
  J->rowvals[23] = 24;
  J->rowvals[24] = 5;
  J->rowvals[25] = 18;
  J->rowvals[26] = 21;
  J->rowvals[27] = 23;
  J->rowvals[28] = 26;
  J->rowvals[29] = 6;
  J->rowvals[30] = 6;
  J->rowvals[31] = 7;
  J->rowvals[32] = 11;
  J->rowvals[33] = 20;
  J->rowvals[34] = 22;
  J->rowvals[35] = 24;
  J->rowvals[36] = 25;
  J->rowvals[37] = 8;
  J->rowvals[38] = 9;
  J->rowvals[39] = 18;
  J->rowvals[40] = 8;
  J->rowvals[41] = 9;
  J->rowvals[42] = 10;
  J->rowvals[43] = 11;
  J->rowvals[44] = 18;
  J->rowvals[45] = 26;
  J->rowvals[46] = 9;
  J->rowvals[47] = 10;
  J->rowvals[48] = 20;
  J->rowvals[49] = 26;
  J->rowvals[50] = 11;
  J->rowvals[51] = 20;
  J->rowvals[52] = 25;
  J->rowvals[53] = 0;
  J->rowvals[54] = 1;
  J->rowvals[55] = 5;
  J->rowvals[56] = 8;
  J->rowvals[57] = 9;
  J->rowvals[58] = 10;
  J->rowvals[59] = 11;
  J->rowvals[60] = 12;
  J->rowvals[61] = 13;
  J->rowvals[62] = 16;
  J->rowvals[63] = 17;
  J->rowvals[64] = 18;
  J->rowvals[65] = 19;
  J->rowvals[66] = 20;
  J->rowvals[67] = 21;
  J->rowvals[68] = 23;
  J->rowvals[69] = 25;
  J->rowvals[70] = 26;
  J->rowvals[71] = 1;
  J->rowvals[72] = 5;
  J->rowvals[73] = 6;
  J->rowvals[74] = 7;
  J->rowvals[75] = 8;
  J->rowvals[76] = 9;
  J->rowvals[77] = 10;
  J->rowvals[78] = 12;
  J->rowvals[79] = 13;
  J->rowvals[80] = 16;
  J->rowvals[81] = 18;
  J->rowvals[82] = 19;
  J->rowvals[83] = 26;
  J->rowvals[84] = 0;
  J->rowvals[85] = 1;
  J->rowvals[86] = 2;
  J->rowvals[87] = 3;
  J->rowvals[88] = 4;
  J->rowvals[89] = 5;
  J->rowvals[90] = 14;
  J->rowvals[91] = 15;
  J->rowvals[92] = 16;
  J->rowvals[93] = 17;
  J->rowvals[94] = 18;
  J->rowvals[95] = 19;
  J->rowvals[96] = 21;
  J->rowvals[97] = 22;
  J->rowvals[98] = 23;
  J->rowvals[99] = 24;
  J->rowvals[100] = 25;
  J->rowvals[101] = 26;
  J->rowvals[102] = 2;
  J->rowvals[103] = 3;
  J->rowvals[104] = 4;
  J->rowvals[105] = 5;
  J->rowvals[106] = 6;
  J->rowvals[107] = 7;
  J->rowvals[108] = 14;
  J->rowvals[109] = 15;
  J->rowvals[110] = 16;
  J->rowvals[111] = 17;
  J->rowvals[112] = 19;
  J->rowvals[113] = 21;
  J->rowvals[114] = 23;
  J->rowvals[115] = 0;
  J->rowvals[116] = 1;
  J->rowvals[117] = 2;
  J->rowvals[118] = 3;
  J->rowvals[119] = 4;
  J->rowvals[120] = 5;
  J->rowvals[121] = 8;
  J->rowvals[122] = 9;
  J->rowvals[123] = 10;
  J->rowvals[124] = 11;
  J->rowvals[125] = 12;
  J->rowvals[126] = 13;
  J->rowvals[127] = 14;
  J->rowvals[128] = 15;
  J->rowvals[129] = 16;
  J->rowvals[130] = 17;
  J->rowvals[131] = 18;
  J->rowvals[132] = 19;
  J->rowvals[133] = 20;
  J->rowvals[134] = 21;
  J->rowvals[135] = 22;
  J->rowvals[136] = 23;
  J->rowvals[137] = 24;
  J->rowvals[138] = 25;
  J->rowvals[139] = 26;
  J->rowvals[140] = 0;
  J->rowvals[141] = 11;
  J->rowvals[142] = 17;
  J->rowvals[143] = 18;
  J->rowvals[144] = 19;
  J->rowvals[145] = 23;
  J->rowvals[146] = 0;
  J->rowvals[147] = 5;
  J->rowvals[148] = 8;
  J->rowvals[149] = 9;
  J->rowvals[150] = 10;
  J->rowvals[151] = 12;
  J->rowvals[152] = 13;
  J->rowvals[153] = 16;
  J->rowvals[154] = 17;
  J->rowvals[155] = 18;
  J->rowvals[156] = 26;
  J->rowvals[157] = 1;
  J->rowvals[158] = 17;
  J->rowvals[159] = 19;
  J->rowvals[160] = 20;
  J->rowvals[161] = 23;
  J->rowvals[162] = 24;
  J->rowvals[163] = 26;
  J->rowvals[164] = 6;
  J->rowvals[165] = 11;
  J->rowvals[166] = 20;
  J->rowvals[167] = 25;
  J->rowvals[168] = 0;
  J->rowvals[169] = 1;
  J->rowvals[170] = 2;
  J->rowvals[171] = 3;
  J->rowvals[172] = 4;
  J->rowvals[173] = 5;
  J->rowvals[174] = 14;
  J->rowvals[175] = 15;
  J->rowvals[176] = 16;
  J->rowvals[177] = 21;
  J->rowvals[178] = 23;
  J->rowvals[179] = 22;
  J->rowvals[180] = 24;
  J->rowvals[181] = 25;
  J->rowvals[182] = 3;
  J->rowvals[183] = 4;
  J->rowvals[184] = 5;
  J->rowvals[185] = 17;
  J->rowvals[186] = 19;
  J->rowvals[187] = 21;
  J->rowvals[188] = 23;
  J->rowvals[189] = 25;
  J->rowvals[190] = 6;
  J->rowvals[191] = 22;
  J->rowvals[192] = 24;
  J->rowvals[193] = 25;
  J->rowvals[194] = 11;
  J->rowvals[195] = 20;
  J->rowvals[196] = 22;
  J->rowvals[197] = 24;
  J->rowvals[198] = 25;
  J->rowvals[199] = 1;
  J->rowvals[200] = 5;
  J->rowvals[201] = 9;
  J->rowvals[202] = 10;
  J->rowvals[203] = 18;
  J->rowvals[204] = 19;
  J->rowvals[205] = 25;
  J->rowvals[206] = 26;
  J->colptrs[0] = 0;
  J->colptrs[1] = 5;
  J->colptrs[2] = 11;
  J->colptrs[3] = 14;
  J->colptrs[4] = 20;
  J->colptrs[5] = 24;
  J->colptrs[6] = 29;
  J->colptrs[7] = 30;
  J->colptrs[8] = 37;
  J->colptrs[9] = 40;
  J->colptrs[10] = 46;
  J->colptrs[11] = 50;
  J->colptrs[12] = 53;
  J->colptrs[13] = 71;
  J->colptrs[14] = 84;
  J->colptrs[15] = 102;
  J->colptrs[16] = 115;
  J->colptrs[17] = 140;
  J->colptrs[18] = 146;
  J->colptrs[19] = 157;
  J->colptrs[20] = 164;
  J->colptrs[21] = 168;
  J->colptrs[22] = 179;
  J->colptrs[23] = 182;
  J->colptrs[24] = 190;
  J->colptrs[25] = 194;
  J->colptrs[26] = 199;
  J->colptrs[27] = 207;
J->data[0] = -2.0*k[0]*p[0]*x_tmp[16];
J->data[1] = k[0]*p[0]*x_tmp[16];
J->data[2] = k[0]*p[0]*x_tmp[16];
J->data[3] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[4] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[5] = p[1];
J->data[6] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[7] = k[0]*p[0]*x_tmp[16];
J->data[8] = p[1];
J->data[9] = p[2];
J->data[10] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[11] = -2.0*k[0]*p[0]*x_tmp[16];
J->data[12] = k[0]*p[0]*x_tmp[16];
J->data[13] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[14] = 2.0*p[1];
J->data[15] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[16] = 2.0*k[0]*p[0]*x_tmp[16];
J->data[17] = p[1];
J->data[18] = p[2];
J->data[19] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[20] = p[1];
J->data[21] = -2.0*p[1];
J->data[22] = p[1];
J->data[23] = p[2];
J->data[24] = -(1.0*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[25] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[26] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[27] = k[0]*p[0]*x_tmp[14];
J->data[28] = k[0]*p[0]*x_tmp[12];
J->data[29] = -2.0*p[3];
J->data[30] = p[3]/k[0];
J->data[31] = -1.0*p[3];
J->data[32] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
J->data[33] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
J->data[34] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[35] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[36] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[37] = -2.0*k[0]*p[0]*x_tmp[16];
J->data[38] = k[0]*p[0]*x_tmp[16];
J->data[39] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[40] = 2.0*p[1];
J->data[41] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[42] = 2.0*k[0]*p[0]*x_tmp[16];
J->data[43] = p[2];
J->data[44] = p[1];
J->data[45] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[46] = p[1];
J->data[47] = -2.0*p[1];
J->data[48] = p[2];
J->data[49] = p[1];
J->data[50] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[51] = k[0]*p[0]*x_tmp[16];
J->data[52] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[53] = (2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[54] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
J->data[55] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[56] = ((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[57] = (p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
J->data[58] = ((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
J->data[59] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
J->data[60] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[61] = k[0]*p[0]*x_tmp[16];
J->data[62] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[63] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2));
J->data[64] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[65] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2));
J->data[66] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25]))/(pow(k[0],2));
J->data[67] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[68] = (p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
J->data[69] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
J->data[70] = (p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
J->data[71] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[72] = p[1]/k[0];
J->data[73] = p[2]/k[0];
J->data[74] = p[2];
J->data[75] = p[1]/k[0];
J->data[76] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[77] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
J->data[78] = p[1];
J->data[79] = -1.0*p[1];
J->data[80] = p[1];
J->data[81] = p[1]/k[0];
J->data[82] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
J->data[83] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[84] = (2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[85] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2));
J->data[86] = ((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[87] = (p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
J->data[88] = ((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
J->data[89] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[90] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[91] = k[0]*p[0]*x_tmp[16];
J->data[92] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[93] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2));
J->data[94] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[95] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2));
J->data[96] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[97] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
J->data[98] = (p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
J->data[99] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25]))/(pow(k[0],2));
J->data[100] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
J->data[101] = (p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
J->data[102] = p[1]/k[0];
J->data[103] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[104] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
J->data[105] = p[1]/k[0];
J->data[106] = p[2]/k[0];
J->data[107] = p[2];
J->data[108] = p[1];
J->data[109] = -1.0*p[1];
J->data[110] = p[1];
J->data[111] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
J->data[112] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
J->data[113] = p[1]/k[0];
J->data[114] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
J->data[115] = (2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])/(pow(k[0],2));
J->data[116] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2));
J->data[117] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[118] = (p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])/(pow(k[0],2));
J->data[119] = (2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14])/(pow(k[0],2));
J->data[120] = (k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[121] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2));
J->data[122] = (p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[123] = (2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[124] = (p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
J->data[125] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[126] = k[0]*p[0]*x_tmp[12];
J->data[127] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[128] = k[0]*p[0]*x_tmp[14];
J->data[129] = -(1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0];
J->data[130] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2));
J->data[131] = -(1.0*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2));
J->data[132] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2));
J->data[133] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
J->data[134] = -(1.0*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
J->data[135] = (p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])/(pow(k[0],2));
J->data[136] = (p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
J->data[137] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2));
J->data[138] = (p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])/(pow(k[0],2));
J->data[139] = (p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
J->data[140] = p[1];
J->data[141] = p[2];
J->data[142] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[143] = p[1];
J->data[144] = k[0]*p[0]*x_tmp[16];
J->data[145] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[146] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[147] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
J->data[148] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])/(pow(k[0],2));
J->data[149] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12]))/(pow(k[0],2));
J->data[150] = p[0];
J->data[151] = -1.0*k[0]*p[0];
J->data[152] = k[0]*p[0];
J->data[153] = -1.0*k[0]*p[0];
J->data[154] = k[0]*p[0]*x_tmp[14];
J->data[155] = -(1.0*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[156] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[157] = p[1];
J->data[158] = p[1];
J->data[159] = -2.0*p[1];
J->data[160] = p[2];
J->data[161] = p[1];
J->data[162] = p[2];
J->data[163] = p[1];
J->data[164] = 2.0*p[2];
J->data[165] = p[1];
J->data[166] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[167] = p[1];
J->data[168] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[169] = k[0]*p[0]*x_tmp[12];
J->data[170] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])/(pow(k[0],2));
J->data[171] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[172] = p[0];
J->data[173] = ((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
J->data[174] = -1.0*k[0]*p[0];
J->data[175] = k[0]*p[0];
J->data[176] = -1.0*k[0]*p[0];
J->data[177] = -(1.0*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[178] = -(1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[179] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
J->data[180] = k[0]*p[0]*x_tmp[16];
J->data[181] = -1.0*k[0]*p[0]*x_tmp[16];
J->data[182] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[183] = 2.0*k[0]*p[0]*x_tmp[14];
J->data[184] = 2.0*p[1];
J->data[185] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[186] = k[0]*p[0]*x_tmp[12];
J->data[187] = p[1];
J->data[188] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[189] = p[2];
J->data[190] = 2.0*p[2];
J->data[191] = p[1];
J->data[192] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[193] = p[1];
J->data[194] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[195] = k[0]*p[0]*x_tmp[12];
J->data[196] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[197] = k[0]*p[0]*x_tmp[14];
J->data[198] = -(1.0*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
J->data[199] = -1.0*k[0]*p[0]*x_tmp[14];
J->data[200] = 2.0*p[1];
J->data[201] = -1.0*k[0]*p[0]*x_tmp[12];
J->data[202] = 2.0*k[0]*p[0]*x_tmp[12];
J->data[203] = p[1];
J->data[204] = k[0]*p[0]*x_tmp[14];
J->data[205] = p[2];
J->data[206] = -(1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
  return(0);
}


 int JBBand_enhancer_3(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_enhancer_3(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_enhancer_3(long int N, realtype t, N_Vector x,
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
JB->data[0] = 2.0*k[0]*p[0]*x_tmp[16];
JB->data[1] = -1.0*p[1];
JB->data[12] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
JB->data[14] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2));
JB->data[16] = -(1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0]))/(pow(k[0],2));
JB->data[17] = -1.0*p[1];
JB->data[18] = k[0]*p[0]*x_tmp[14];
JB->data[21] = k[0]*p[0]*x_tmp[12];
JB->data[27] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[28] = ((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[39] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
JB->data[40] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
JB->data[41] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
JB->data[43] = (p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1])/(pow(k[0],2));
JB->data[46] = -1.0*p[1];
JB->data[48] = -1.0*k[0]*p[0]*x_tmp[12];
JB->data[53] = k[0]*p[0]*x_tmp[14];
JB->data[56] = 2.0*k[0]*p[0]*x_tmp[16];
JB->data[57] = -2.0*p[1];
JB->data[68] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
JB->data[69] = -(1.0*p[1])/k[0];
JB->data[70] = (2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])/(pow(k[0],2));
JB->data[75] = -(1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2));
JB->data[83] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[84] = ((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[85] = -1.0*p[1];
JB->data[95] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2));
JB->data[96] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
JB->data[97] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2));
JB->data[102] = ((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])/(pow(k[0],2));
JB->data[104] = k[0]*p[0]*x_tmp[14];
JB->data[111] = -2.0*k[0]*p[0]*x_tmp[16];
JB->data[112] = 2.0*p[1];
JB->data[122] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2));
JB->data[123] = (2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
JB->data[124] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2));
JB->data[129] = -1.0*p[0];
JB->data[131] = -2.0*k[0]*p[0]*x_tmp[14];
JB->data[140] = (2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14])/(pow(k[0],2));
JB->data[147] = (2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[148] = -(1.0*p[1])/k[0];
JB->data[149] = (2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[150] = -(1.0*p[1])/k[0];
JB->data[151] = -(1.0*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
JB->data[153] = -(1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
JB->data[156] = -(1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2));
JB->data[158] = -2.0*p[1];
JB->data[161] = -2.0*p[1];
JB->data[168] = 2.0*p[3];
JB->data[169] = -(1.0*p[3])/k[0];
JB->data[175] = -(1.0*p[2])/k[0];
JB->data[177] = -(1.0*p[2])/k[0];
JB->data[182] = -2.0*p[2];
JB->data[186] = -2.0*p[2];
JB->data[196] = p[3];
JB->data[202] = -1.0*p[2];
JB->data[204] = -1.0*p[2];
JB->data[224] = 2.0*k[0]*p[0]*x_tmp[16];
JB->data[225] = -2.0*p[1];
JB->data[228] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2));
JB->data[229] = -(1.0*p[1])/k[0];
JB->data[232] = (2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])/(pow(k[0],2));
JB->data[234] = -(1.0*((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12]))/(pow(k[0],2));
JB->data[251] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[252] = ((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[253] = -1.0*p[1];
JB->data[255] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2));
JB->data[256] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
JB->data[259] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2));
JB->data[261] = ((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])/(pow(k[0],2));
JB->data[269] = k[0]*p[0]*x_tmp[12];
JB->data[279] = -2.0*k[0]*p[0]*x_tmp[16];
JB->data[280] = 2.0*p[1];
JB->data[282] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2));
JB->data[283] = (2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])/(pow(k[0],2));
JB->data[286] = -(1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2));
JB->data[288] = -1.0*p[0];
JB->data[296] = -2.0*k[0]*p[0]*x_tmp[12];
JB->data[304] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
JB->data[306] = -1.0*p[2];
JB->data[308] = ((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[309] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25]))/(pow(k[0],2));
JB->data[313] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11]))/(pow(k[0],2));
JB->data[314] = -1.0*p[2];
JB->data[317] = -1.0*p[1];
JB->data[322] = k[0]*p[0]*x_tmp[12];
JB->data[336] = k[0]*p[0]*x_tmp[16];
JB->data[337] = -1.0*p[1];
JB->data[340] = k[0]*p[0]*x_tmp[12];
JB->data[342] = k[0]*p[0];
JB->data[363] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[364] = p[1];
JB->data[367] = -1.0*k[0]*p[0]*x_tmp[12];
JB->data[369] = -1.0*k[0]*p[0];
JB->data[392] = k[0]*p[0]*x_tmp[16];
JB->data[393] = -1.0*p[1];
JB->data[394] = k[0]*p[0]*x_tmp[14];
JB->data[399] = k[0]*p[0];
JB->data[419] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[420] = p[1];
JB->data[421] = -1.0*k[0]*p[0]*x_tmp[14];
JB->data[426] = -1.0*k[0]*p[0];
JB->data[444] = k[0]*p[0]*x_tmp[16];
JB->data[445] = -1.0*p[1];
JB->data[446] = k[0]*p[0]*x_tmp[16];
JB->data[447] = -1.0*p[1];
JB->data[448] = (k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14])/k[0];
JB->data[450] = k[0]*p[0];
JB->data[453] = k[0]*p[0];
JB->data[459] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[471] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
JB->data[473] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
JB->data[474] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
JB->data[475] = (p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17])/(pow(k[0],2));
JB->data[476] = ((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[477] = -1.0*k[0]*p[0]*x_tmp[14];
JB->data[478] = -1.0*p[1];
JB->data[482] = k[0]*p[0]*x_tmp[12];
JB->data[486] = k[0]*p[0]*x_tmp[16];
JB->data[491] = k[0]*p[0]*x_tmp[12];
JB->data[494] = k[0]*p[0]*x_tmp[16];
JB->data[495] = -1.0*p[1];
JB->data[498] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2));
JB->data[499] = -(1.0*p[1])/k[0];
JB->data[500] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2));
JB->data[502] = ((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18])/(pow(k[0],2));
JB->data[503] = -1.0*p[1];
JB->data[504] = ((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[512] = -1.0*p[1];
JB->data[514] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[525] = (p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])/(pow(k[0],2));
JB->data[526] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
JB->data[527] = (p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])/(pow(k[0],2));
JB->data[528] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
JB->data[529] = (p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17])/(pow(k[0],2));
JB->data[530] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[532] = 2.0*p[1];
JB->data[536] = -1.0*k[0]*p[0]*x_tmp[12];
JB->data[539] = -1.0*k[0]*p[0]*x_tmp[14];
JB->data[547] = (p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2));
JB->data[550] = -1.0*p[2];
JB->data[551] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[552] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
JB->data[556] = (p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])/(pow(k[0],2));
JB->data[559] = -1.0*p[2];
JB->data[560] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])/(pow(k[0],2));
JB->data[565] = -1.0*k[0]*p[0]*x_tmp[12];
JB->data[567] = k[0]*p[0]*x_tmp[16];
JB->data[568] = -1.0*p[1];
JB->data[569] = k[0]*p[0]*x_tmp[16];
JB->data[570] = -1.0*p[1];
JB->data[572] = k[0]*p[0]*x_tmp[14];
JB->data[579] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
JB->data[581] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
JB->data[582] = -(1.0*p[1])/k[0];
JB->data[583] = ((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21])/(pow(k[0],2));
JB->data[588] = ((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[590] = -1.0*p[1];
JB->data[595] = -1.0*p[2];
JB->data[597] = -1.0*p[2];
JB->data[601] = -(1.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
JB->data[608] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25]))/(pow(k[0],2));
JB->data[610] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2));
JB->data[616] = ((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[618] = -1.0*p[1];
JB->data[619] = k[0]*p[0]*x_tmp[14];
JB->data[624] = k[0]*p[0]*x_tmp[16];
JB->data[625] = -1.0*p[1];
JB->data[626] = -1.0*k[0]*p[0]*x_tmp[14];
JB->data[633] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2));
JB->data[635] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2));
JB->data[636] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
JB->data[637] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
JB->data[638] = k[0]*p[0]*x_tmp[16];
JB->data[640] = -1.0*p[1];
JB->data[642] = ((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[644] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14])/(pow(k[0],2));
JB->data[652] = -1.0*p[2];
JB->data[655] = (p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2));
JB->data[662] = (p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])/(pow(k[0],2));
JB->data[664] = (p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])/(pow(k[0],2));
JB->data[667] = -1.0*p[2];
JB->data[670] = -1.0*k[0]*p[0]*x_tmp[16];
JB->data[672] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])/(pow(k[0],2));
JB->data[673] = -1.0*k[0]*p[0]*x_tmp[14];
JB->data[682] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16])))/(pow(k[0],2));
JB->data[686] = k[0]*p[0]*x_tmp[16];
JB->data[687] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25]))/(pow(k[0],2));
JB->data[689] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25]))/(pow(k[0],2));
JB->data[691] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2));
JB->data[695] = -1.0*p[1];
JB->data[697] = k[0]*p[0]*x_tmp[16];
JB->data[698] = -1.0*p[2];
JB->data[699] = -1.0*p[1];
JB->data[700] = ((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14])/(pow(k[0],2));
JB->data[701] = -1.0*p[2];
JB->data[703] = k[0]*p[0]*x_tmp[16];
JB->data[707] = -1.0*k[0]*p[0]*x_tmp[12];
JB->data[711] = k[0]*p[0]*x_tmp[16];
JB->data[712] = -1.0*p[1];
JB->data[714] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2));
JB->data[715] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2));
JB->data[716] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2));
JB->data[718] = -(1.0*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2));
JB->data[720] = ((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])/(pow(k[0],2));
JB->data[721] = -1.0*p[1];
JB->data[728] = ((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14])/(pow(k[0],2));

  for (iJ=0; iJ<729; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_enhancer_3(realtype t, N_Vector x,
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
  JB->rowvals[2] = 12;
  JB->rowvals[3] = 14;
  JB->rowvals[4] = 16;
  JB->rowvals[5] = 17;
  JB->rowvals[6] = 18;
  JB->rowvals[7] = 21;
  JB->rowvals[8] = 0;
  JB->rowvals[9] = 1;
  JB->rowvals[10] = 12;
  JB->rowvals[11] = 13;
  JB->rowvals[12] = 14;
  JB->rowvals[13] = 16;
  JB->rowvals[14] = 19;
  JB->rowvals[15] = 21;
  JB->rowvals[16] = 26;
  JB->rowvals[17] = 2;
  JB->rowvals[18] = 3;
  JB->rowvals[19] = 14;
  JB->rowvals[20] = 15;
  JB->rowvals[21] = 16;
  JB->rowvals[22] = 21;
  JB->rowvals[23] = 2;
  JB->rowvals[24] = 3;
  JB->rowvals[25] = 4;
  JB->rowvals[26] = 14;
  JB->rowvals[27] = 15;
  JB->rowvals[28] = 16;
  JB->rowvals[29] = 21;
  JB->rowvals[30] = 23;
  JB->rowvals[31] = 3;
  JB->rowvals[32] = 4;
  JB->rowvals[33] = 14;
  JB->rowvals[34] = 15;
  JB->rowvals[35] = 16;
  JB->rowvals[36] = 21;
  JB->rowvals[37] = 23;
  JB->rowvals[38] = 5;
  JB->rowvals[39] = 12;
  JB->rowvals[40] = 13;
  JB->rowvals[41] = 14;
  JB->rowvals[42] = 15;
  JB->rowvals[43] = 16;
  JB->rowvals[44] = 18;
  JB->rowvals[45] = 21;
  JB->rowvals[46] = 23;
  JB->rowvals[47] = 26;
  JB->rowvals[48] = 6;
  JB->rowvals[49] = 7;
  JB->rowvals[50] = 13;
  JB->rowvals[51] = 15;
  JB->rowvals[52] = 20;
  JB->rowvals[53] = 24;
  JB->rowvals[54] = 7;
  JB->rowvals[55] = 13;
  JB->rowvals[56] = 15;
  JB->rowvals[57] = 8;
  JB->rowvals[58] = 9;
  JB->rowvals[59] = 12;
  JB->rowvals[60] = 13;
  JB->rowvals[61] = 16;
  JB->rowvals[62] = 18;
  JB->rowvals[63] = 8;
  JB->rowvals[64] = 9;
  JB->rowvals[65] = 10;
  JB->rowvals[66] = 12;
  JB->rowvals[67] = 13;
  JB->rowvals[68] = 16;
  JB->rowvals[69] = 18;
  JB->rowvals[70] = 26;
  JB->rowvals[71] = 9;
  JB->rowvals[72] = 10;
  JB->rowvals[73] = 12;
  JB->rowvals[74] = 13;
  JB->rowvals[75] = 16;
  JB->rowvals[76] = 18;
  JB->rowvals[77] = 26;
  JB->rowvals[78] = 7;
  JB->rowvals[79] = 9;
  JB->rowvals[80] = 11;
  JB->rowvals[81] = 12;
  JB->rowvals[82] = 16;
  JB->rowvals[83] = 17;
  JB->rowvals[84] = 20;
  JB->rowvals[85] = 25;
  JB->rowvals[86] = 12;
  JB->rowvals[87] = 13;
  JB->rowvals[88] = 16;
  JB->rowvals[89] = 18;
  JB->rowvals[90] = 12;
  JB->rowvals[91] = 13;
  JB->rowvals[92] = 16;
  JB->rowvals[93] = 18;
  JB->rowvals[94] = 14;
  JB->rowvals[95] = 15;
  JB->rowvals[96] = 16;
  JB->rowvals[97] = 21;
  JB->rowvals[98] = 14;
  JB->rowvals[99] = 15;
  JB->rowvals[100] = 16;
  JB->rowvals[101] = 21;
  JB->rowvals[102] = 12;
  JB->rowvals[103] = 13;
  JB->rowvals[104] = 14;
  JB->rowvals[105] = 15;
  JB->rowvals[106] = 16;
  JB->rowvals[107] = 18;
  JB->rowvals[108] = 21;
  JB->rowvals[109] = 0;
  JB->rowvals[110] = 12;
  JB->rowvals[111] = 14;
  JB->rowvals[112] = 15;
  JB->rowvals[113] = 16;
  JB->rowvals[114] = 17;
  JB->rowvals[115] = 18;
  JB->rowvals[116] = 19;
  JB->rowvals[117] = 23;
  JB->rowvals[118] = 0;
  JB->rowvals[119] = 5;
  JB->rowvals[120] = 8;
  JB->rowvals[121] = 9;
  JB->rowvals[122] = 12;
  JB->rowvals[123] = 13;
  JB->rowvals[124] = 14;
  JB->rowvals[125] = 16;
  JB->rowvals[126] = 17;
  JB->rowvals[127] = 18;
  JB->rowvals[128] = 26;
  JB->rowvals[129] = 1;
  JB->rowvals[130] = 12;
  JB->rowvals[131] = 13;
  JB->rowvals[132] = 14;
  JB->rowvals[133] = 15;
  JB->rowvals[134] = 16;
  JB->rowvals[135] = 17;
  JB->rowvals[136] = 19;
  JB->rowvals[137] = 23;
  JB->rowvals[138] = 26;
  JB->rowvals[139] = 7;
  JB->rowvals[140] = 10;
  JB->rowvals[141] = 11;
  JB->rowvals[142] = 12;
  JB->rowvals[143] = 16;
  JB->rowvals[144] = 19;
  JB->rowvals[145] = 20;
  JB->rowvals[146] = 25;
  JB->rowvals[147] = 0;
  JB->rowvals[148] = 1;
  JB->rowvals[149] = 2;
  JB->rowvals[150] = 3;
  JB->rowvals[151] = 5;
  JB->rowvals[152] = 12;
  JB->rowvals[153] = 14;
  JB->rowvals[154] = 15;
  JB->rowvals[155] = 16;
  JB->rowvals[156] = 21;
  JB->rowvals[157] = 23;
  JB->rowvals[158] = 1;
  JB->rowvals[159] = 3;
  JB->rowvals[160] = 7;
  JB->rowvals[161] = 14;
  JB->rowvals[162] = 16;
  JB->rowvals[163] = 22;
  JB->rowvals[164] = 24;
  JB->rowvals[165] = 25;
  JB->rowvals[166] = 3;
  JB->rowvals[167] = 4;
  JB->rowvals[168] = 5;
  JB->rowvals[169] = 12;
  JB->rowvals[170] = 14;
  JB->rowvals[171] = 15;
  JB->rowvals[172] = 16;
  JB->rowvals[173] = 17;
  JB->rowvals[174] = 19;
  JB->rowvals[175] = 21;
  JB->rowvals[176] = 23;
  JB->rowvals[177] = 4;
  JB->rowvals[178] = 7;
  JB->rowvals[179] = 14;
  JB->rowvals[180] = 16;
  JB->rowvals[181] = 19;
  JB->rowvals[182] = 22;
  JB->rowvals[183] = 24;
  JB->rowvals[184] = 25;
  JB->rowvals[185] = 7;
  JB->rowvals[186] = 11;
  JB->rowvals[187] = 12;
  JB->rowvals[188] = 14;
  JB->rowvals[189] = 16;
  JB->rowvals[190] = 20;
  JB->rowvals[191] = 22;
  JB->rowvals[192] = 23;
  JB->rowvals[193] = 24;
  JB->rowvals[194] = 25;
  JB->rowvals[195] = 26;
  JB->rowvals[196] = 1;
  JB->rowvals[197] = 5;
  JB->rowvals[198] = 9;
  JB->rowvals[199] = 10;
  JB->rowvals[200] = 12;
  JB->rowvals[201] = 13;
  JB->rowvals[202] = 14;
  JB->rowvals[203] = 16;
  JB->rowvals[204] = 18;
  JB->rowvals[205] = 19;
  JB->rowvals[206] = 26;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 8;
  JB->colptrs[2] = 17;
  JB->colptrs[3] = 23;
  JB->colptrs[4] = 31;
  JB->colptrs[5] = 38;
  JB->colptrs[6] = 48;
  JB->colptrs[7] = 54;
  JB->colptrs[8] = 57;
  JB->colptrs[9] = 63;
  JB->colptrs[10] = 71;
  JB->colptrs[11] = 78;
  JB->colptrs[12] = 86;
  JB->colptrs[13] = 90;
  JB->colptrs[14] = 94;
  JB->colptrs[15] = 98;
  JB->colptrs[16] = 102;
  JB->colptrs[17] = 109;
  JB->colptrs[18] = 118;
  JB->colptrs[19] = 129;
  JB->colptrs[20] = 139;
  JB->colptrs[21] = 147;
  JB->colptrs[22] = 158;
  JB->colptrs[23] = 166;
  JB->colptrs[24] = 177;
  JB->colptrs[25] = 185;
  JB->colptrs[26] = 196;
  JB->colptrs[27] = 207;
  return(0);
}


 int sx_enhancer_3(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
  memset(sxdot_tmp,0,sizeof(double)*27);
  switch (plist[ip]) {
  case 0: {
sxdot_tmp[0] = p[1]*sx_tmp[1] + p[1]*sx_tmp[17] + ((pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[12]*x_tmp[21] - 2.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 2.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 6.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*sx_tmp[16])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[21]*x_tmp[12];
sxdot_tmp[1] = p[1]*sx_tmp[19] - (1.0*((pow(k[0],3))*x_tmp[14]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[21] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*sx_tmp[1]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[21]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[3] + ((pow(k[0],2))*x_tmp[21] - 2.0*(pow(k[0],3))*x_tmp[2]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[14]*x_tmp[21])/(pow(k[0],2)) + (p[1]*sx_tmp[15])/k[0] + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16];
sxdot_tmp[3] = p[1]*sx_tmp[4] - (1.0*((pow(k[0],2))*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[2]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[14]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[3]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[4] = p[0]*sx_tmp[21] - 2.0*p[1]*sx_tmp[4] + ((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[15]*x_tmp[21] + 2.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + 2.0*k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + 2.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) - 6.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16])/(pow(k[0],2)) + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[23] + 2.0*p[1]*sx_tmp[26] + ((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[21] - 2.0*(pow(k[0],3))*x_tmp[5]*x_tmp[12] - 2.0*(pow(k[0],3))*x_tmp[5]*x_tmp[14] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[16]*x_tmp[18] - 2.0*(pow(k[0],3))*x_tmp[16]*x_tmp[21])/(pow(k[0],2)) + (p[1]*sx_tmp[13])/k[0] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[5]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[14])/(pow(k[0],2)) + (sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[20] - 2.0*p[3]*sx_tmp[6] + 2.0*p[2]*sx_tmp[24] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[13])/k[0] + (p[2]*sx_tmp[15])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[13] - 1.0*p[3]*sx_tmp[7] + p[2]*sx_tmp[15];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + ((pow(k[0],2))*x_tmp[18] - 2.0*(pow(k[0],3))*x_tmp[8]*x_tmp[16] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[12]*x_tmp[18])/(pow(k[0],2)) + (p[1]*sx_tmp[13])/k[0] + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[9] = p[1]*sx_tmp[10] - (1.0*((pow(k[0],2))*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[8]*x_tmp[16] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[12]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[8]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[10] = p[0]*sx_tmp[18] - 2.0*p[1]*sx_tmp[10] + ((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 2.0*(pow(k[0],3))*x_tmp[13]*x_tmp[18] + 2.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + 2.0*k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + 2.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 6.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16])/(pow(k[0],2)) + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[11] = p[2]*sx_tmp[9] + p[2]*sx_tmp[17] + p[1]*sx_tmp[20] - (1.0*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[18] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[11]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[12] = p[1]*sx_tmp[13] - (1.0*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[13] = ((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16])/k[0] - 1.0*p[1]*sx_tmp[13] + k[0]*p[0]*sx_tmp[18] + k[0]*p[0]*sx_tmp[12]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[14] = p[1]*sx_tmp[15] - (1.0*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[15] = ((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16])/k[0] - 1.0*p[1]*sx_tmp[15] + k[0]*p[0]*sx_tmp[21] + k[0]*p[0]*sx_tmp[14]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[16] = p[1]*sx_tmp[13] + p[1]*sx_tmp[15] - (1.0*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]))/k[0] - (1.0*sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16];
sxdot_tmp[17] = p[1]*sx_tmp[19] - (1.0*((pow(k[0],3))*x_tmp[12]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[18] - 1.0*k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) - 1.0*k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[12];
sxdot_tmp[18] = p[1]*sx_tmp[9] + p[1]*sx_tmp[17] + p[1]*sx_tmp[26] - (1.0*((pow(k[0],3))*x_tmp[5]*x_tmp[12] - 1.0*(pow(k[0],2))*x_tmp[18] + (pow(k[0],3))*x_tmp[8]*x_tmp[16] - 1.0*(pow(k[0],2))*x_tmp[12]*x_tmp[16] + (pow(k[0],3))*x_tmp[12]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[12]*x_tmp[21] + (pow(k[0],3))*x_tmp[16]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (p[1]*sx_tmp[13])/k[0] - (1.0*sx_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[19] = (k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16])/(pow(k[0],2)) - 2.0*p[1]*sx_tmp[19] - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[1]*x_tmp[16] + k[0]*p[0]*sx_tmp[17]*x_tmp[16] + k[0]*p[0]*sx_tmp[23]*x_tmp[12] + k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[20] = p[2]*sx_tmp[10] + p[2]*sx_tmp[19] + (k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[18] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[11]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[21] = p[1]*sx_tmp[1] + p[1]*sx_tmp[3] + p[1]*sx_tmp[23] - (1.0*((pow(k[0],3))*x_tmp[2]*x_tmp[16] - 1.0*(pow(k[0],2))*x_tmp[21] + (pow(k[0],3))*x_tmp[5]*x_tmp[14] - 1.0*(pow(k[0],2))*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[14]*x_tmp[18] + (pow(k[0],3))*x_tmp[14]*x_tmp[21] + (pow(k[0],3))*x_tmp[16]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[21]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[14];
sxdot_tmp[22] = p[2]*sx_tmp[1] + p[2]*sx_tmp[3] + p[1]*sx_tmp[24] - (1.0*(k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[22]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[23] = p[1]*sx_tmp[4] + p[1]*sx_tmp[19] - (1.0*((pow(k[0],2))*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[5]*x_tmp[14] + (pow(k[0],2))*x_tmp[14]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[15]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[16]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + k[0]*x_tmp[15]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[15]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[14]*x_tmp[15]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[15]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*sx_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[17]*x_tmp[16];
sxdot_tmp[24] = p[2]*sx_tmp[4] + p[2]*sx_tmp[19] + (k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[21] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[22]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[25] = p[1]*sx_tmp[20] + p[1]*sx_tmp[24] + p[2]*sx_tmp[23] + p[2]*sx_tmp[26] - (1.0*(k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[7]*x_tmp[18] + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[7]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[12]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[7]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[25]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*sx_tmp[7])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[22]*x_tmp[16];
sxdot_tmp[26] = p[1]*sx_tmp[10] + p[1]*sx_tmp[19] - (1.0*((pow(k[0],2))*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[5]*x_tmp[12] + (pow(k[0],2))*x_tmp[12]*x_tmp[16] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[18] - 1.0*(pow(k[0],3))*x_tmp[13]*x_tmp[21] - 1.0*(pow(k[0],3))*x_tmp[16]*x_tmp[18] + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + k[0]*x_tmp[16]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + k[0]*x_tmp[13]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + k[0]*x_tmp[12]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + k[0]*x_tmp[14]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) - 3.0*(pow(k[0],3))*x_tmp[12]*x_tmp[13]*x_tmp[16] - 3.0*(pow(k[0],3))*x_tmp[13]*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[26]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[1]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16];

  } break;

  case 1: {
sxdot_tmp[0] = p[1]*sx_tmp[1] + p[1]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[17])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*sx_tmp[16])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[21]*x_tmp[12];
sxdot_tmp[1] = p[1]*sx_tmp[19] - (1.0*((pow(k[0],2))*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*sx_tmp[1]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[21]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[3] + (k[0]*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[3])/(pow(k[0],2)) + (p[1]*sx_tmp[15])/k[0] + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16];
sxdot_tmp[3] = p[1]*sx_tmp[4] - (1.0*(k[0]*x_tmp[15] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[4]))/(pow(k[0],2)) - (1.0*sx_tmp[3]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[4] = p[0]*sx_tmp[21] - 2.0*p[1]*sx_tmp[4] + (k[0]*x_tmp[15] - 2.0*(pow(k[0],2))*x_tmp[4])/(pow(k[0],2)) + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[23] + 2.0*p[1]*sx_tmp[26] + (k[0]*x_tmp[13] + k[0]*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[23] + 2.0*(pow(k[0],2))*x_tmp[26])/(pow(k[0],2)) + (p[1]*sx_tmp[13])/k[0] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[5]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[14])/(pow(k[0],2)) + (sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[20] - 2.0*p[3]*sx_tmp[6] + 2.0*p[2]*sx_tmp[24] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[13])/k[0] + (p[2]*sx_tmp[15])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[13] - 1.0*p[3]*sx_tmp[7] + p[2]*sx_tmp[15];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + (k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[9])/(pow(k[0],2)) + (p[1]*sx_tmp[13])/k[0] + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[9] = p[1]*sx_tmp[10] - (1.0*(k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[9] - 1.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[8]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[10] = p[0]*sx_tmp[18] - 2.0*p[1]*sx_tmp[10] + (k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[10])/(pow(k[0],2)) + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[11] = p[2]*sx_tmp[9] + p[2]*sx_tmp[17] + p[1]*sx_tmp[20] + x_tmp[20] - (1.0*sx_tmp[11]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[12] = p[1]*sx_tmp[13] + x_tmp[13] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[13] = k[0]*p[0]*sx_tmp[18] - 1.0*x_tmp[13] - 1.0*p[1]*sx_tmp[13] + k[0]*p[0]*sx_tmp[12]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[14] = p[1]*sx_tmp[15] + x_tmp[15] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[15] = k[0]*p[0]*sx_tmp[21] - 1.0*x_tmp[15] - 1.0*p[1]*sx_tmp[15] + k[0]*p[0]*sx_tmp[14]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[16] = p[1]*sx_tmp[13] + p[1]*sx_tmp[15] + (k[0]*x_tmp[13] + k[0]*x_tmp[15])/k[0] - (1.0*sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16];
sxdot_tmp[17] = p[1]*sx_tmp[19] - (1.0*((pow(k[0],2))*x_tmp[17] - 1.0*(pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[12];
sxdot_tmp[18] = p[1]*sx_tmp[9] + p[1]*sx_tmp[17] + p[1]*sx_tmp[26] + (k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[26])/(pow(k[0],2)) + (p[1]*sx_tmp[13])/k[0] - (1.0*sx_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[19] = k[0]*p[0]*sx_tmp[1]*x_tmp[16] - 2.0*x_tmp[19] - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 2.0*p[1]*sx_tmp[19] + k[0]*p[0]*sx_tmp[17]*x_tmp[16] + k[0]*p[0]*sx_tmp[23]*x_tmp[12] + k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[20] = p[2]*sx_tmp[10] + p[2]*sx_tmp[19] - 1.0*x_tmp[20] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[11]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[21] = p[1]*sx_tmp[1] + p[1]*sx_tmp[3] + p[1]*sx_tmp[23] + (k[0]*x_tmp[15] + (pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[23])/(pow(k[0],2)) + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[21]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[14];
sxdot_tmp[22] = p[2]*sx_tmp[1] + p[2]*sx_tmp[3] + p[1]*sx_tmp[24] + x_tmp[24] - (1.0*sx_tmp[22]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[23] = p[1]*sx_tmp[4] + p[1]*sx_tmp[19] - (1.0*(k[0]*x_tmp[15] - 1.0*(pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[23]))/(pow(k[0],2)) + (sx_tmp[15]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*sx_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[17]*x_tmp[16];
sxdot_tmp[24] = p[2]*sx_tmp[4] + p[2]*sx_tmp[19] - 1.0*x_tmp[24] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[22]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[25] = p[1]*sx_tmp[20] + p[1]*sx_tmp[24] + p[2]*sx_tmp[23] + p[2]*sx_tmp[26] + ((pow(k[0],2))*x_tmp[20] + (pow(k[0],2))*x_tmp[24])/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[25]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*sx_tmp[7])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[22]*x_tmp[16];
sxdot_tmp[26] = p[1]*sx_tmp[10] + p[1]*sx_tmp[19] - (1.0*(k[0]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[10] - 1.0*(pow(k[0],2))*x_tmp[19] + (pow(k[0],2))*x_tmp[26]))/(pow(k[0],2)) + (sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[26]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[1]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16];

  } break;

  case 2: {
sxdot_tmp[0] = p[1]*sx_tmp[1] + p[1]*sx_tmp[17] + ((2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*sx_tmp[16])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[21]*x_tmp[12];
sxdot_tmp[1] = p[1]*sx_tmp[19] - (1.0*sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*sx_tmp[1]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[21]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[3] + (p[1]*sx_tmp[15])/k[0] + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16];
sxdot_tmp[3] = p[1]*sx_tmp[4] - (1.0*sx_tmp[3]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[4] = p[0]*sx_tmp[21] - 2.0*p[1]*sx_tmp[4] + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[23] + 2.0*p[1]*sx_tmp[26] + (p[1]*sx_tmp[13])/k[0] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[5]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[14])/(pow(k[0],2)) + (sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[20] - 2.0*p[3]*sx_tmp[6] + 2.0*p[2]*sx_tmp[24] + (k[0]*x_tmp[13] + k[0]*x_tmp[15] + 2.0*(pow(k[0],2))*x_tmp[20] + 2.0*(pow(k[0],2))*x_tmp[24])/(pow(k[0],2)) + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[13])/k[0] + (p[2]*sx_tmp[15])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[13] - 1.0*p[3]*sx_tmp[7] + p[2]*sx_tmp[15] + (k[0]*x_tmp[13] + k[0]*x_tmp[15])/k[0];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + (p[1]*sx_tmp[13])/k[0] + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[9] = p[1]*sx_tmp[10] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[8]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[10] = p[0]*sx_tmp[18] - 2.0*p[1]*sx_tmp[10] + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[11] = p[2]*sx_tmp[9] + p[2]*sx_tmp[17] + p[1]*sx_tmp[20] + ((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[17])/(pow(k[0],2)) - (1.0*sx_tmp[11]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[12] = p[1]*sx_tmp[13] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[13] = k[0]*p[0]*sx_tmp[18] - 1.0*p[1]*sx_tmp[13] + k[0]*p[0]*sx_tmp[12]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[14] = p[1]*sx_tmp[15] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[15] = k[0]*p[0]*sx_tmp[21] - 1.0*p[1]*sx_tmp[15] + k[0]*p[0]*sx_tmp[14]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[16] = p[1]*sx_tmp[13] + p[1]*sx_tmp[15] - (1.0*sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16];
sxdot_tmp[17] = p[1]*sx_tmp[19] - (1.0*sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[12];
sxdot_tmp[18] = p[1]*sx_tmp[9] + p[1]*sx_tmp[17] + p[1]*sx_tmp[26] + (p[1]*sx_tmp[13])/k[0] - (1.0*sx_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[19] = k[0]*p[0]*sx_tmp[1]*x_tmp[16] - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 2.0*p[1]*sx_tmp[19] + k[0]*p[0]*sx_tmp[17]*x_tmp[16] + k[0]*p[0]*sx_tmp[23]*x_tmp[12] + k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[20] = p[2]*sx_tmp[10] + p[2]*sx_tmp[19] + ((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[19])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[11]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[21] = p[1]*sx_tmp[1] + p[1]*sx_tmp[3] + p[1]*sx_tmp[23] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[21]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[14];
sxdot_tmp[22] = p[2]*sx_tmp[1] + p[2]*sx_tmp[3] + p[1]*sx_tmp[24] + ((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[3])/(pow(k[0],2)) - (1.0*sx_tmp[22]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[23] = p[1]*sx_tmp[4] + p[1]*sx_tmp[19] + (sx_tmp[15]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*sx_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[17]*x_tmp[16];
sxdot_tmp[24] = p[2]*sx_tmp[4] + p[2]*sx_tmp[19] + ((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[19])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[22]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[25] = p[1]*sx_tmp[20] + p[1]*sx_tmp[24] + p[2]*sx_tmp[23] + p[2]*sx_tmp[26] + ((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[26])/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[25]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*sx_tmp[7])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[22]*x_tmp[16];
sxdot_tmp[26] = p[1]*sx_tmp[10] + p[1]*sx_tmp[19] + (sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[26]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[1]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16];

  } break;

  case 3: {
sxdot_tmp[0] = p[1]*sx_tmp[1] + p[1]*sx_tmp[17] + ((2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*sx_tmp[16])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[21]*x_tmp[12];
sxdot_tmp[1] = p[1]*sx_tmp[19] - (1.0*sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*sx_tmp[1]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[21]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[3] + (p[1]*sx_tmp[15])/k[0] + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16];
sxdot_tmp[3] = p[1]*sx_tmp[4] - (1.0*sx_tmp[3]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[4] = p[0]*sx_tmp[21] - 2.0*p[1]*sx_tmp[4] + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[23] + 2.0*p[1]*sx_tmp[26] + (p[1]*sx_tmp[13])/k[0] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[5]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[14])/(pow(k[0],2)) + (sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[20] - 2.0*p[3]*sx_tmp[6] + 2.0*p[2]*sx_tmp[24] + (k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[6])/(pow(k[0],2)) + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[13])/k[0] + (p[2]*sx_tmp[15])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[13] - 1.0*p[3]*sx_tmp[7] + p[2]*sx_tmp[15] - 1.0*x_tmp[7];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + (p[1]*sx_tmp[13])/k[0] + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[9] = p[1]*sx_tmp[10] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[8]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[10] = p[0]*sx_tmp[18] - 2.0*p[1]*sx_tmp[10] + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[11] = p[2]*sx_tmp[9] + p[2]*sx_tmp[17] + p[1]*sx_tmp[20] - 1.0*x_tmp[11] - (1.0*sx_tmp[11]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[12] = p[1]*sx_tmp[13] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[13] = k[0]*p[0]*sx_tmp[18] - 1.0*p[1]*sx_tmp[13] + k[0]*p[0]*sx_tmp[12]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[14] = p[1]*sx_tmp[15] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[15] = k[0]*p[0]*sx_tmp[21] - 1.0*p[1]*sx_tmp[15] + k[0]*p[0]*sx_tmp[14]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[16] = p[1]*sx_tmp[13] + p[1]*sx_tmp[15] - (1.0*sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16];
sxdot_tmp[17] = p[1]*sx_tmp[19] - (1.0*sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[12];
sxdot_tmp[18] = p[1]*sx_tmp[9] + p[1]*sx_tmp[17] + p[1]*sx_tmp[26] + (p[1]*sx_tmp[13])/k[0] - (1.0*sx_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[19] = k[0]*p[0]*sx_tmp[1]*x_tmp[16] - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 2.0*p[1]*sx_tmp[19] + k[0]*p[0]*sx_tmp[17]*x_tmp[16] + k[0]*p[0]*sx_tmp[23]*x_tmp[12] + k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[20] = p[2]*sx_tmp[10] + p[2]*sx_tmp[19] - 1.0*x_tmp[20] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[11]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[21] = p[1]*sx_tmp[1] + p[1]*sx_tmp[3] + p[1]*sx_tmp[23] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[21]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[14];
sxdot_tmp[22] = p[2]*sx_tmp[1] + p[2]*sx_tmp[3] + p[1]*sx_tmp[24] - 1.0*x_tmp[22] - (1.0*sx_tmp[22]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[23] = p[1]*sx_tmp[4] + p[1]*sx_tmp[19] + (sx_tmp[15]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*sx_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[17]*x_tmp[16];
sxdot_tmp[24] = p[2]*sx_tmp[4] + p[2]*sx_tmp[19] - 1.0*x_tmp[24] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[22]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[25] = p[1]*sx_tmp[20] + p[1]*sx_tmp[24] + p[2]*sx_tmp[23] + p[2]*sx_tmp[26] - 1.0*x_tmp[25] + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[25]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*sx_tmp[7])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[22]*x_tmp[16];
sxdot_tmp[26] = p[1]*sx_tmp[10] + p[1]*sx_tmp[19] + (sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[26]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[1]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16];

  } break;

  case 4: {
sxdot_tmp[0] = p[1]*sx_tmp[1] + p[1]*sx_tmp[17] + ((2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*sx_tmp[16])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[21]*x_tmp[12];
sxdot_tmp[1] = p[1]*sx_tmp[19] - (1.0*sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*sx_tmp[1]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[21]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[3] + (p[1]*sx_tmp[15])/k[0] + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16];
sxdot_tmp[3] = p[1]*sx_tmp[4] - (1.0*sx_tmp[3]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[4] = p[0]*sx_tmp[21] - 2.0*p[1]*sx_tmp[4] + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[23] + 2.0*p[1]*sx_tmp[26] + 144.0/(pow(k[0],2)) + (p[1]*sx_tmp[13])/k[0] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[5]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[14])/(pow(k[0],2)) + (sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[20] - 2.0*p[3]*sx_tmp[6] + 2.0*p[2]*sx_tmp[24] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[13])/k[0] + (p[2]*sx_tmp[15])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[13] - 1.0*p[3]*sx_tmp[7] + p[2]*sx_tmp[15];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + (p[1]*sx_tmp[13])/k[0] + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[9] = p[1]*sx_tmp[10] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[8]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[10] = p[0]*sx_tmp[18] - 2.0*p[1]*sx_tmp[10] + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[11] = p[2]*sx_tmp[9] + p[2]*sx_tmp[17] + p[1]*sx_tmp[20] - (1.0*sx_tmp[11]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[12] = p[1]*sx_tmp[13] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[13] = k[0]*p[0]*sx_tmp[18] - 1.0*p[1]*sx_tmp[13] + k[0]*p[0]*sx_tmp[12]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[14] = p[1]*sx_tmp[15] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[15] = k[0]*p[0]*sx_tmp[21] - 1.0*p[1]*sx_tmp[15] + k[0]*p[0]*sx_tmp[14]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[16] = p[1]*sx_tmp[13] + p[1]*sx_tmp[15] + 12.0/k[0] - (1.0*sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16];
sxdot_tmp[17] = p[1]*sx_tmp[19] - (1.0*sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[12];
sxdot_tmp[18] = p[1]*sx_tmp[9] + p[1]*sx_tmp[17] + p[1]*sx_tmp[26] + (p[1]*sx_tmp[13])/k[0] - (1.0*sx_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[19] = k[0]*p[0]*sx_tmp[1]*x_tmp[16] - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 2.0*p[1]*sx_tmp[19] + k[0]*p[0]*sx_tmp[17]*x_tmp[16] + k[0]*p[0]*sx_tmp[23]*x_tmp[12] + k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[20] = p[2]*sx_tmp[10] + p[2]*sx_tmp[19] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[11]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[21] = p[1]*sx_tmp[1] + p[1]*sx_tmp[3] + p[1]*sx_tmp[23] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[21]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[14];
sxdot_tmp[22] = p[2]*sx_tmp[1] + p[2]*sx_tmp[3] + p[1]*sx_tmp[24] - (1.0*sx_tmp[22]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[23] = p[1]*sx_tmp[4] + p[1]*sx_tmp[19] + (sx_tmp[15]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*sx_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[17]*x_tmp[16];
sxdot_tmp[24] = p[2]*sx_tmp[4] + p[2]*sx_tmp[19] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[22]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[25] = p[1]*sx_tmp[20] + p[1]*sx_tmp[24] + p[2]*sx_tmp[23] + p[2]*sx_tmp[26] + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[25]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*sx_tmp[7])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[22]*x_tmp[16];
sxdot_tmp[26] = p[1]*sx_tmp[10] + p[1]*sx_tmp[19] + (sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[26]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[1]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16];

  } break;

  case 5: {
sxdot_tmp[0] = p[1]*sx_tmp[1] + p[1]*sx_tmp[17] + ((2.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 2.0*(pow(k[0],3))*p[0]*x_tmp[0])*sx_tmp[16])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + ((2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[21]*x_tmp[12];
sxdot_tmp[1] = p[1]*sx_tmp[19] - (1.0*sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*sx_tmp[1]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[1]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[21]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[2] = 2.0*p[1]*sx_tmp[3] + (p[1]*sx_tmp[15])/k[0] + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16];
sxdot_tmp[3] = p[1]*sx_tmp[4] - (1.0*sx_tmp[3]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[21] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[14])*sx_tmp[21])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[4] = p[0]*sx_tmp[21] - 2.0*p[1]*sx_tmp[4] + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[3] - 2.0*p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + (pow(k[0],2))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (sx_tmp[14]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[15])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[23]*x_tmp[14];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[23] + 2.0*p[1]*sx_tmp[26] + (k[0]*x_tmp[16] - 2.0*(pow(k[0],2))*x_tmp[5])/(pow(k[0],2)) + (p[1]*sx_tmp[13])/k[0] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[5]*(2.0*(pow(k[0],2))*p[5] + 2.0*(pow(k[0],3))*p[0]*x_tmp[12] + 2.0*(pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16])*sx_tmp[14])/(pow(k[0],2)) + (sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18] - 2.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2));
sxdot_tmp[6] = 2.0*p[2]*sx_tmp[20] - 2.0*p[3]*sx_tmp[6] + 2.0*p[2]*sx_tmp[24] + (p[3]*sx_tmp[7])/k[0] + (p[2]*sx_tmp[13])/k[0] + (p[2]*sx_tmp[15])/k[0];
sxdot_tmp[7] = p[2]*sx_tmp[13] - 1.0*p[3]*sx_tmp[7] + p[2]*sx_tmp[15];
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[9] + (p[1]*sx_tmp[13])/k[0] + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (((pow(k[0],2))*p[0] - 2.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12])*sx_tmp[16])/(pow(k[0],2)) - 2.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[9] = p[1]*sx_tmp[10] - (1.0*sx_tmp[9]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] + (pow(k[0],3))*p[0]*x_tmp[18] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[12])*sx_tmp[18])/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[8]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[10] = p[0]*sx_tmp[18] - 2.0*p[1]*sx_tmp[10] + (sx_tmp[16]*(2.0*(pow(k[0],3))*p[0]*x_tmp[9] - 2.0*p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) + (pow(k[0],2))*p[0]*x_tmp[12]))/(pow(k[0],2)) + (sx_tmp[12]*((pow(k[0],2))*p[0]*x_tmp[16] - 2.0*p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + 2.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*(2.0*p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) + 2.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16] + 2.0*k[0]*p[0]*sx_tmp[26]*x_tmp[12];
sxdot_tmp[11] = p[2]*sx_tmp[9] + p[2]*sx_tmp[17] + p[1]*sx_tmp[20] - (1.0*sx_tmp[11]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[12] = p[1]*sx_tmp[13] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[13] = k[0]*p[0]*sx_tmp[18] - 1.0*p[1]*sx_tmp[13] + k[0]*p[0]*sx_tmp[12]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[12];
sxdot_tmp[14] = p[1]*sx_tmp[15] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[15] = k[0]*p[0]*sx_tmp[21] - 1.0*p[1]*sx_tmp[15] + k[0]*p[0]*sx_tmp[14]*x_tmp[16] + k[0]*p[0]*sx_tmp[16]*x_tmp[14];
sxdot_tmp[16] = p[1]*sx_tmp[13] + p[1]*sx_tmp[15] - 1.0*x_tmp[16] - (1.0*sx_tmp[16]*(k[0]*p[5] + (pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],2))*p[0]*x_tmp[14]))/k[0] - 1.0*k[0]*p[0]*sx_tmp[18] - 1.0*k[0]*p[0]*sx_tmp[21] - 1.0*k[0]*p[0]*sx_tmp[12]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[14]*x_tmp[16];
sxdot_tmp[17] = p[1]*sx_tmp[19] - (1.0*sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*sx_tmp[17]*((pow(k[0],2))*p[1] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) - 1.0*p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[0] + (pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[0]*x_tmp[16] + k[0]*p[0]*sx_tmp[18]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[23]*x_tmp[12];
sxdot_tmp[18] = p[1]*sx_tmp[9] + p[1]*sx_tmp[17] + p[1]*sx_tmp[26] - 1.0*x_tmp[18] + (p[1]*sx_tmp[13])/k[0] - (1.0*sx_tmp[18]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[18])*sx_tmp[14])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[8]*x_tmp[16];
sxdot_tmp[19] = k[0]*p[0]*sx_tmp[1]*x_tmp[16] - (1.0*sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17]))/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[15]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[13]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 2.0*p[1]*sx_tmp[19] + k[0]*p[0]*sx_tmp[17]*x_tmp[16] + k[0]*p[0]*sx_tmp[23]*x_tmp[12] + k[0]*p[0]*sx_tmp[26]*x_tmp[14];
sxdot_tmp[20] = p[2]*sx_tmp[10] + p[2]*sx_tmp[19] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[20])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[11]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[12];
sxdot_tmp[21] = p[1]*sx_tmp[1] + p[1]*sx_tmp[3] + p[1]*sx_tmp[23] - 1.0*x_tmp[21] + (p[1]*sx_tmp[15])/k[0] - (1.0*sx_tmp[21]*((pow(k[0],2))*p[5] - 1.0*(pow(k[0],2))*p[0] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) - (1.0*sx_tmp[16]*((pow(k[0],3))*p[0]*x_tmp[0] - 1.0*p[0]*((pow(k[0],3))*x_tmp[0] - 1.0*k[0]*((pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*x_tmp[12]*x_tmp[14]) + (pow(k[0],3))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[21])*sx_tmp[12])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[0]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[2]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[5]*x_tmp[14];
sxdot_tmp[22] = p[2]*sx_tmp[1] + p[2]*sx_tmp[3] + p[1]*sx_tmp[24] - (1.0*sx_tmp[22]*((pow(k[0],2))*p[3] + (pow(k[0],3))*p[0]*x_tmp[16]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) + (p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[23] = p[1]*sx_tmp[4] + p[1]*sx_tmp[19] - 1.0*x_tmp[23] + (sx_tmp[15]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[3] - 1.0*k[0]*((pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14]*x_tmp[15]) + (pow(k[0],3))*x_tmp[14]*x_tmp[15]) + p[0]*((pow(k[0],3))*x_tmp[17] - 1.0*k[0]*((pow(k[0],2))*x_tmp[17] + (pow(k[0],2))*x_tmp[12]*x_tmp[15]) + (pow(k[0],3))*x_tmp[12]*x_tmp[15]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],3))*p[0]*x_tmp[17] + (pow(k[0],3))*p[0]*x_tmp[21]))/(pow(k[0],2)) + (sx_tmp[14]*(p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[23]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[21])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[23] - 1.0*k[0]*((pow(k[0],2))*x_tmp[23] + (pow(k[0],2))*x_tmp[15]*x_tmp[16]) + (pow(k[0],3))*x_tmp[15]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[23])*sx_tmp[12])/(pow(k[0],2)) - (1.0*sx_tmp[23]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[3]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[14] - 1.0*k[0]*p[0]*sx_tmp[17]*x_tmp[16];
sxdot_tmp[24] = p[2]*sx_tmp[4] + p[2]*sx_tmp[19] - (1.0*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[24])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[22])*sx_tmp[16])/(pow(k[0],2)) - (1.0*(p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))/(pow(k[0],2)) + k[0]*p[0]*sx_tmp[22]*x_tmp[16] + k[0]*p[0]*sx_tmp[25]*x_tmp[14];
sxdot_tmp[25] = p[1]*sx_tmp[20] + p[1]*sx_tmp[24] + p[2]*sx_tmp[23] + p[2]*sx_tmp[26] - 1.0*x_tmp[25] + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[11] - 1.0*k[0]*((pow(k[0],2))*x_tmp[11] + (pow(k[0],2))*x_tmp[7]*x_tmp[12]) + (pow(k[0],3))*x_tmp[7]*x_tmp[12]) + p[0]*((pow(k[0],3))*x_tmp[22] - 1.0*k[0]*((pow(k[0],2))*x_tmp[22] + (pow(k[0],2))*x_tmp[7]*x_tmp[14]) + (pow(k[0],3))*x_tmp[7]*x_tmp[14]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[11] - 1.0*(pow(k[0],3))*p[0]*x_tmp[22]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[12])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[25] - 1.0*k[0]*((pow(k[0],2))*x_tmp[25] + (pow(k[0],2))*x_tmp[7]*x_tmp[16]) + (pow(k[0],3))*x_tmp[7]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[25])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[25]*((pow(k[0],2))*p[3] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]))*sx_tmp[7])/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[11]*x_tmp[16] - 1.0*k[0]*p[0]*sx_tmp[22]*x_tmp[16];
sxdot_tmp[26] = p[1]*sx_tmp[10] + p[1]*sx_tmp[19] - 1.0*x_tmp[26] + (sx_tmp[13]*(p[0]*((pow(k[0],3))*x_tmp[18] - 1.0*k[0]*((pow(k[0],2))*x_tmp[18] + (pow(k[0],2))*x_tmp[12]*x_tmp[16]) + (pow(k[0],3))*x_tmp[12]*x_tmp[16]) + p[0]*((pow(k[0],3))*x_tmp[21] - 1.0*k[0]*((pow(k[0],2))*x_tmp[21] + (pow(k[0],2))*x_tmp[14]*x_tmp[16]) + (pow(k[0],3))*x_tmp[14]*x_tmp[16]) - 1.0*k[0]*p[1]))/(pow(k[0],2)) + (sx_tmp[16]*(p[0]*((pow(k[0],3))*x_tmp[1] - 1.0*k[0]*((pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[13]*x_tmp[14]) + (pow(k[0],3))*x_tmp[13]*x_tmp[14]) + p[0]*((pow(k[0],3))*x_tmp[9] - 1.0*k[0]*((pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[12]*x_tmp[13]) + (pow(k[0],3))*x_tmp[12]*x_tmp[13]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[1] - 1.0*(pow(k[0],3))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[18]))/(pow(k[0],2)) + (sx_tmp[12]*(p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*p[0]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],3))*p[0]*x_tmp[26]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*p[0] - 1.0*(pow(k[0],3))*p[0]*x_tmp[16])*sx_tmp[18])/(pow(k[0],2)) + ((p[0]*((pow(k[0],3))*x_tmp[26] - 1.0*k[0]*((pow(k[0],2))*x_tmp[26] + (pow(k[0],2))*x_tmp[13]*x_tmp[16]) + (pow(k[0],3))*x_tmp[13]*x_tmp[16]) - 1.0*(pow(k[0],3))*p[0]*x_tmp[26])*sx_tmp[14])/(pow(k[0],2)) - (1.0*sx_tmp[26]*((pow(k[0],2))*p[1] + (pow(k[0],2))*p[5] + (pow(k[0],3))*p[0]*x_tmp[12] + (pow(k[0],3))*p[0]*x_tmp[14]))/(pow(k[0],2)) - 1.0*k[0]*p[0]*sx_tmp[1]*x_tmp[16] + k[0]*p[0]*sx_tmp[5]*x_tmp[12] - 1.0*k[0]*p[0]*sx_tmp[9]*x_tmp[16];

  } break;

  }
 for (ix=0; ix<27; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_enhancer_3(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  memset(sx0_tmp,0,sizeof(double)*27);
  switch (ip) {
  }

  return;
}


void y_enhancer_3(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*7];
y[it+nt*1] = x[it+nt*6];
    
    return;
}


void dydp_enhancer_3(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }
  
  return;
}


void dydx_enhancer_3(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*54);
dydx[13] = 1.0;
dydx[14] = 1.0;
  
  return;
}


void sy_enhancer_3(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(6+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(6+np*nx)];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(6+np*nx)];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(6+np*nx)];

  } break;

  case 4: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(6+np*nx)];

  } break;

  case 5: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(7+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(6+np*nx)];

  } break;

  }
  
  return;
}
int root_enhancer_3(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_enhancer_3(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root_enhancer_3(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval_enhancer_3(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval_enhancer_3(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
void deltadisc_enhancer_3(double t, int idisc, N_Vector x, void *user_data){
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double deltadisc[27];
  memset(deltadisc,0,sizeof(double)*27);
  for(ix = 0; ix<27;ix++){;
  x_tmp[ix] += deltadisc[ix];
  };
}
void sdeltadisc_enhancer_3(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data){
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
  double deltadisc[27];
  double *sdeltadisc;
  memset(deltadisc,0,sizeof(double)*27);
  sdeltadisc = mxMalloc(sizeof(double)*27*np);
  memset(sdeltadisc,0,sizeof(double)*27*np);
  for (ip=0; ip<np; ip++) {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
     switch (plist[ip]) {
     }
  }
  for(ip = 0; ip<np;ip++){
      sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
      for(ix = 0; ix<27;ix++){
      sx_tmp[ix] += sdeltadisc[plist[ip]+np*ix];
     }
  }
  for(ix = 0; ix<27;ix++){
  x_tmp[ix] += deltadisc[ix];
  };
 mxFree(sdeltadisc);
}


void dxdotdp_enhancer_3(double t, int nt, int it, double *dxdotdp, double *p, double *k, double *u, double *x, int *plist, int np, int nx){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
dxdotdp[(0+ip*nx)] = ((pow(k[0],3))*x[it+nt*14]*x[it+nt*18] + (pow(k[0],3))*x[it+nt*12]*x[it+nt*21] - 2.0*k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*0] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*14]) - 2.0*k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) - 2.0*k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + 6.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*14]*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = -(1.0*((pow(k[0],3))*x[it+nt*14]*x[it+nt*18] - 1.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*21] - 1.0*k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*0] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*14]) + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*1] + (pow(k[0],2))*x[it+nt*13]*x[it+nt*14]) - 1.0*k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) - 1.0*k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*26] + (pow(k[0],2))*x[it+nt*13]*x[it+nt*16]) + 3.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*14]*x[it+nt*16] - 3.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*14]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(2+ip*nx)] = ((pow(k[0],2))*x[it+nt*21] - 2.0*(pow(k[0],3))*x[it+nt*2]*x[it+nt*16] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16] - 2.0*(pow(k[0],3))*x[it+nt*14]*x[it+nt*21])/(pow(k[0],2));
dxdotdp[(3+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*21] - 1.0*(pow(k[0],3))*x[it+nt*2]*x[it+nt*16] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16] - 1.0*(pow(k[0],3))*x[it+nt*14]*x[it+nt*21] - 1.0*(pow(k[0],3))*x[it+nt*15]*x[it+nt*21] + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*3] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*15]) + k[0]*x[it+nt*15]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*23] + (pow(k[0],2))*x[it+nt*15]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*14]*x[it+nt*15]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = ((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16] - 2.0*(pow(k[0],3))*x[it+nt*15]*x[it+nt*21] + 2.0*k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*3] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*15]) + 2.0*k[0]*x[it+nt*15]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + 2.0*k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*23] + (pow(k[0],2))*x[it+nt*15]*x[it+nt*16]) - 6.0*(pow(k[0],3))*x[it+nt*14]*x[it+nt*15]*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = ((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*21] - 2.0*(pow(k[0],3))*x[it+nt*5]*x[it+nt*12] - 2.0*(pow(k[0],3))*x[it+nt*5]*x[it+nt*14] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16] - 2.0*(pow(k[0],3))*x[it+nt*16]*x[it+nt*18] - 2.0*(pow(k[0],3))*x[it+nt*16]*x[it+nt*21])/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = ((pow(k[0],2))*x[it+nt*18] - 2.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*16] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16] - 2.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*18])/(pow(k[0],2));
dxdotdp[(9+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*18] - 1.0*(pow(k[0],3))*x[it+nt*8]*x[it+nt*16] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16] - 1.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*18] - 1.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*18] + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*13]) + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*26] + (pow(k[0],2))*x[it+nt*13]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*13]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = ((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16] - 2.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*18] + 2.0*k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*13]) + 2.0*k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + 2.0*k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*26] + (pow(k[0],2))*x[it+nt*13]*x[it+nt*16]) - 6.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*13]*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(11+ip*nx)] = -(1.0*(k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*11] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*12]) - 1.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*18] + k[0]*x[it+nt*7]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*25] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*12]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(12+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]))/k[0];
dxdotdp[(13+ip*nx)] = ((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16])/k[0];
dxdotdp[(14+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]))/k[0];
dxdotdp[(15+ip*nx)] = ((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16])/k[0];
dxdotdp[(16+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]))/k[0];
dxdotdp[(17+ip*nx)] = -(1.0*((pow(k[0],3))*x[it+nt*12]*x[it+nt*21] - 1.0*(pow(k[0],3))*x[it+nt*15]*x[it+nt*18] - 1.0*k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*0] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*14]) - 1.0*k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*17] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*15]) + k[0]*x[it+nt*15]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) - 1.0*k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*23] + (pow(k[0],2))*x[it+nt*15]*x[it+nt*16]) + 3.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*14]*x[it+nt*16] - 3.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*15]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(18+ip*nx)] = -(1.0*((pow(k[0],3))*x[it+nt*5]*x[it+nt*12] - 1.0*(pow(k[0],2))*x[it+nt*18] + (pow(k[0],3))*x[it+nt*8]*x[it+nt*16] - 1.0*(pow(k[0],2))*x[it+nt*12]*x[it+nt*16] + (pow(k[0],3))*x[it+nt*12]*x[it+nt*18] - 1.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*21] + (pow(k[0],3))*x[it+nt*16]*x[it+nt*18] + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*0] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*14]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*14]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(19+ip*nx)] = (k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*1] + (pow(k[0],2))*x[it+nt*13]*x[it+nt*14]) - 1.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*21] - 1.0*(pow(k[0],3))*x[it+nt*15]*x[it+nt*18] + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*17] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*15]) + k[0]*x[it+nt*15]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*23] + (pow(k[0],2))*x[it+nt*15]*x[it+nt*16]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*26] + (pow(k[0],2))*x[it+nt*13]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*15]*x[it+nt*16] - 3.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*14]*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(20+ip*nx)] = (k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*11] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*12]) - 1.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*18] + k[0]*x[it+nt*7]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*25] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*12]*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(21+ip*nx)] = -(1.0*((pow(k[0],3))*x[it+nt*2]*x[it+nt*16] - 1.0*(pow(k[0],2))*x[it+nt*21] + (pow(k[0],3))*x[it+nt*5]*x[it+nt*14] - 1.0*(pow(k[0],2))*x[it+nt*14]*x[it+nt*16] - 1.0*(pow(k[0],3))*x[it+nt*14]*x[it+nt*18] + (pow(k[0],3))*x[it+nt*14]*x[it+nt*21] + (pow(k[0],3))*x[it+nt*16]*x[it+nt*21] + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*0] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*14]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*14]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(22+ip*nx)] = -(1.0*(k[0]*x[it+nt*7]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) - 1.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*21] + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*22] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*14]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*25] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*14]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(23+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*21] - 1.0*(pow(k[0],3))*x[it+nt*5]*x[it+nt*14] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16] - 1.0*(pow(k[0],3))*x[it+nt*15]*x[it+nt*18] - 1.0*(pow(k[0],3))*x[it+nt*15]*x[it+nt*21] - 1.0*(pow(k[0],3))*x[it+nt*16]*x[it+nt*21] + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*3] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*15]) + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*17] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*15]) + k[0]*x[it+nt*15]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*23] + (pow(k[0],2))*x[it+nt*15]*x[it+nt*16]) + k[0]*x[it+nt*15]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*23] + (pow(k[0],2))*x[it+nt*15]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*15]*x[it+nt*16] - 3.0*(pow(k[0],3))*x[it+nt*14]*x[it+nt*15]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(24+ip*nx)] = (k[0]*x[it+nt*7]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) - 1.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*21] + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*22] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*14]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*25] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*14]*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(25+ip*nx)] = -(1.0*(k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*11] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*12]) - 1.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*21] - 1.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*18] + k[0]*x[it+nt*7]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*7]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*22] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*14]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*25] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*16]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*25] + (pow(k[0],2))*x[it+nt*7]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*12]*x[it+nt*16] - 3.0*(pow(k[0],3))*x[it+nt*7]*x[it+nt*14]*x[it+nt*16]))/(pow(k[0],2));
dxdotdp[(26+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*18] - 1.0*(pow(k[0],3))*x[it+nt*5]*x[it+nt*12] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16] - 1.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*18] - 1.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*21] - 1.0*(pow(k[0],3))*x[it+nt*16]*x[it+nt*18] + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*1] + (pow(k[0],2))*x[it+nt*13]*x[it+nt*14]) + k[0]*x[it+nt*16]*((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*13]) + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*18] + (pow(k[0],2))*x[it+nt*12]*x[it+nt*16]) + k[0]*x[it+nt*13]*((pow(k[0],2))*x[it+nt*21] + (pow(k[0],2))*x[it+nt*14]*x[it+nt*16]) + k[0]*x[it+nt*12]*((pow(k[0],2))*x[it+nt*26] + (pow(k[0],2))*x[it+nt*13]*x[it+nt*16]) + k[0]*x[it+nt*14]*((pow(k[0],2))*x[it+nt*26] + (pow(k[0],2))*x[it+nt*13]*x[it+nt*16]) - 3.0*(pow(k[0],3))*x[it+nt*12]*x[it+nt*13]*x[it+nt*16] - 3.0*(pow(k[0],3))*x[it+nt*13]*x[it+nt*14]*x[it+nt*16]))/(pow(k[0],2));

  } break;

  case 1: {
dxdotdp[(0+ip*nx)] = ((pow(k[0],2))*x[it+nt*1] + (pow(k[0],2))*x[it+nt*17])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*1] - 1.0*(pow(k[0],2))*x[it+nt*19]))/(pow(k[0],2));
dxdotdp[(2+ip*nx)] = (k[0]*x[it+nt*15] + 2.0*(pow(k[0],2))*x[it+nt*3])/(pow(k[0],2));
dxdotdp[(3+ip*nx)] = -(1.0*(k[0]*x[it+nt*15] + (pow(k[0],2))*x[it+nt*3] - 1.0*(pow(k[0],2))*x[it+nt*4]))/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = (k[0]*x[it+nt*15] - 2.0*(pow(k[0],2))*x[it+nt*4])/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = (k[0]*x[it+nt*13] + k[0]*x[it+nt*15] + 2.0*(pow(k[0],2))*x[it+nt*23] + 2.0*(pow(k[0],2))*x[it+nt*26])/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = (k[0]*x[it+nt*13] + 2.0*(pow(k[0],2))*x[it+nt*9])/(pow(k[0],2));
dxdotdp[(9+ip*nx)] = -(1.0*(k[0]*x[it+nt*13] + (pow(k[0],2))*x[it+nt*9] - 1.0*(pow(k[0],2))*x[it+nt*10]))/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = (k[0]*x[it+nt*13] - 2.0*(pow(k[0],2))*x[it+nt*10])/(pow(k[0],2));
dxdotdp[(11+ip*nx)] = x[it+nt*20];
dxdotdp[(12+ip*nx)] = x[it+nt*13];
dxdotdp[(13+ip*nx)] = -1.0*x[it+nt*13];
dxdotdp[(14+ip*nx)] = x[it+nt*15];
dxdotdp[(15+ip*nx)] = -1.0*x[it+nt*15];
dxdotdp[(16+ip*nx)] = (k[0]*x[it+nt*13] + k[0]*x[it+nt*15])/k[0];
dxdotdp[(17+ip*nx)] = -(1.0*((pow(k[0],2))*x[it+nt*17] - 1.0*(pow(k[0],2))*x[it+nt*19]))/(pow(k[0],2));
dxdotdp[(18+ip*nx)] = (k[0]*x[it+nt*13] + (pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*17] + (pow(k[0],2))*x[it+nt*26])/(pow(k[0],2));
dxdotdp[(19+ip*nx)] = -2.0*x[it+nt*19];
dxdotdp[(20+ip*nx)] = -1.0*x[it+nt*20];
dxdotdp[(21+ip*nx)] = (k[0]*x[it+nt*15] + (pow(k[0],2))*x[it+nt*1] + (pow(k[0],2))*x[it+nt*3] + (pow(k[0],2))*x[it+nt*23])/(pow(k[0],2));
dxdotdp[(22+ip*nx)] = x[it+nt*24];
dxdotdp[(23+ip*nx)] = -(1.0*(k[0]*x[it+nt*15] - 1.0*(pow(k[0],2))*x[it+nt*4] - 1.0*(pow(k[0],2))*x[it+nt*19] + (pow(k[0],2))*x[it+nt*23]))/(pow(k[0],2));
dxdotdp[(24+ip*nx)] = -1.0*x[it+nt*24];
dxdotdp[(25+ip*nx)] = ((pow(k[0],2))*x[it+nt*20] + (pow(k[0],2))*x[it+nt*24])/(pow(k[0],2));
dxdotdp[(26+ip*nx)] = -(1.0*(k[0]*x[it+nt*13] - 1.0*(pow(k[0],2))*x[it+nt*10] - 1.0*(pow(k[0],2))*x[it+nt*19] + (pow(k[0],2))*x[it+nt*26]))/(pow(k[0],2));

  } break;

  case 2: {
dxdotdp[(6+ip*nx)] = (k[0]*x[it+nt*13] + k[0]*x[it+nt*15] + 2.0*(pow(k[0],2))*x[it+nt*20] + 2.0*(pow(k[0],2))*x[it+nt*24])/(pow(k[0],2));
dxdotdp[(7+ip*nx)] = (k[0]*x[it+nt*13] + k[0]*x[it+nt*15])/k[0];
dxdotdp[(11+ip*nx)] = ((pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*17])/(pow(k[0],2));
dxdotdp[(20+ip*nx)] = ((pow(k[0],2))*x[it+nt*10] + (pow(k[0],2))*x[it+nt*19])/(pow(k[0],2));
dxdotdp[(22+ip*nx)] = ((pow(k[0],2))*x[it+nt*1] + (pow(k[0],2))*x[it+nt*3])/(pow(k[0],2));
dxdotdp[(24+ip*nx)] = ((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*19])/(pow(k[0],2));
dxdotdp[(25+ip*nx)] = ((pow(k[0],2))*x[it+nt*23] + (pow(k[0],2))*x[it+nt*26])/(pow(k[0],2));

  } break;

  case 3: {
dxdotdp[(6+ip*nx)] = (k[0]*x[it+nt*7] - 2.0*(pow(k[0],2))*x[it+nt*6])/(pow(k[0],2));
dxdotdp[(7+ip*nx)] = -1.0*x[it+nt*7];
dxdotdp[(11+ip*nx)] = -1.0*x[it+nt*11];
dxdotdp[(20+ip*nx)] = -1.0*x[it+nt*20];
dxdotdp[(22+ip*nx)] = -1.0*x[it+nt*22];
dxdotdp[(24+ip*nx)] = -1.0*x[it+nt*24];
dxdotdp[(25+ip*nx)] = -1.0*x[it+nt*25];

  } break;

  case 4: {
dxdotdp[(5+ip*nx)] = 144.0/(pow(k[0],2));
dxdotdp[(16+ip*nx)] = 12.0/k[0];

  } break;

  case 5: {
dxdotdp[(5+ip*nx)] = (k[0]*x[it+nt*16] - 2.0*(pow(k[0],2))*x[it+nt*5])/(pow(k[0],2));
dxdotdp[(16+ip*nx)] = -1.0*x[it+nt*16];
dxdotdp[(18+ip*nx)] = -1.0*x[it+nt*18];
dxdotdp[(21+ip*nx)] = -1.0*x[it+nt*21];
dxdotdp[(23+ip*nx)] = -1.0*x[it+nt*23];
dxdotdp[(25+ip*nx)] = -1.0*x[it+nt*25];
dxdotdp[(26+ip*nx)] = -1.0*x[it+nt*26];

  } break;

  }
  }
  
  return;
}
