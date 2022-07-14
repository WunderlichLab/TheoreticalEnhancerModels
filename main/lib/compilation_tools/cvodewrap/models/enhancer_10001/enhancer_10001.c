#include "enhancer_10001.h"
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


 int xdot_enhancer_10001(realtype t, N_Vector x, N_Vector xdot, void *user_data)
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
xdot_tmp[0] = (2.0*(pow(k[0],2))*p[1]*x_tmp[14] - 4.0*(pow(k[0],2))*p[0]*x_tmp[0] + 2.0*(pow(k[0],2))*p[1]*x_tmp[15] + 2.0*k[0]*p[0]*x_tmp[1] + k[0]*p[1]*x_tmp[11] + k[0]*p[1]*x_tmp[12])/(pow(k[0],2));
xdot_tmp[1] = (k[0]*p[1]*x_tmp[11] - 2.0*k[0]*p[0]*x_tmp[1] + k[0]*p[1]*x_tmp[12])/k[0];
xdot_tmp[2] = (2.0*(pow(k[0],2))*p[2]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[2] - 1.0*(pow(k[0],2))*p[3]*x_tmp[2] + (pow(k[0],2))*p[1]*x_tmp[10] + (pow(k[0],2))*p[2]*x_tmp[14] + (pow(k[0],2))*p[2]*x_tmp[15] + (pow(k[0],2))*p[1]*x_tmp[18])/(pow(k[0],2));
xdot_tmp[3] = ((pow(k[0],2))*p[1]*x_tmp[9] - 2.0*(pow(k[0],2))*p[1]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],2))*p[0]*x_tmp[15] + (pow(k[0],2))*p[1]*x_tmp[17])/(pow(k[0],2));
xdot_tmp[4] = ((pow(k[0],2))*p[1]*x_tmp[9] - 2.0*(pow(k[0],2))*p[1]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[4] + (pow(k[0],2))*p[0]*x_tmp[14] + (pow(k[0],2))*p[0]*x_tmp[15] + (pow(k[0],2))*p[1]*x_tmp[17])/(pow(k[0],2));
xdot_tmp[5] = (2.0*(pow(k[0],2))*p[0]*x_tmp[9] - 4.0*(pow(k[0],2))*p[1]*x_tmp[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[17] + k[0]*p[0]*x_tmp[11] + k[0]*p[0]*x_tmp[12] + 2.0*k[0]*p[1]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[6] = (k[0]*p[2]*x_tmp[11] - 1.0*k[0]*p[3]*x_tmp[6] + k[0]*p[2]*x_tmp[12] + 2.0*k[0]*p[2]*x_tmp[13])/k[0];
xdot_tmp[7] = (4.0*(pow(k[0],2))*p[2]*x_tmp[8] - 2.0*(pow(k[0],2))*p[3]*x_tmp[7] + 2.0*(pow(k[0],2))*p[2]*x_tmp[10] + 2.0*(pow(k[0],2))*p[2]*x_tmp[18] + k[0]*p[3]*x_tmp[6] + k[0]*p[2]*x_tmp[11] + k[0]*p[2]*x_tmp[12] + 2.0*k[0]*p[2]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[8] = (2.0*(pow(k[0],2))*p[2]*x_tmp[5] - 2.0*(pow(k[0],2))*p[1]*x_tmp[8] + (pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[2]*x_tmp[9] - 1.0*(pow(k[0],2))*p[3]*x_tmp[8] + (pow(k[0],2))*p[0]*x_tmp[18] + (pow(k[0],2))*p[2]*x_tmp[17])/(pow(k[0],2));
xdot_tmp[9] = ((pow(k[0],2))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[4] + (pow(k[0],2))*p[1]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[9] - 3.0*(pow(k[0],2))*p[1]*x_tmp[9] + (pow(k[0],2))*p[0]*x_tmp[19] - 1.0*k[0]*p[0]*x_tmp[11] - 1.0*k[0]*p[1]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[10] = ((pow(k[0],2))*p[0]*x_tmp[2] + (pow(k[0],2))*p[2]*x_tmp[4] + (pow(k[0],2))*p[1]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[1]*x_tmp[10] + 2.0*(pow(k[0],2))*p[2]*x_tmp[9] - 1.0*(pow(k[0],2))*p[3]*x_tmp[10] + (pow(k[0],2))*p[2]*x_tmp[19])/(pow(k[0],2));
xdot_tmp[11] = (k[0]*p[0]*x_tmp[1] - 1.0*k[0]*p[0]*x_tmp[11] - 1.0*k[0]*p[1]*x_tmp[11] + k[0]*p[1]*x_tmp[13])/k[0];
xdot_tmp[12] = (k[0]*p[0]*x_tmp[1] - 1.0*k[0]*p[0]*x_tmp[12] - 1.0*k[0]*p[1]*x_tmp[12] + k[0]*p[1]*x_tmp[13])/k[0];
xdot_tmp[13] = (k[0]*p[0]*x_tmp[11] + k[0]*p[0]*x_tmp[12] - 2.0*k[0]*p[1]*x_tmp[13])/k[0];
xdot_tmp[14] = ((pow(k[0],2))*p[0]*x_tmp[0] + (pow(k[0],2))*p[1]*x_tmp[3] + (pow(k[0],2))*p[1]*x_tmp[4] - 3.0*(pow(k[0],2))*p[0]*x_tmp[14] - 1.0*(pow(k[0],2))*p[1]*x_tmp[14] + (pow(k[0],2))*p[1]*x_tmp[19] - 1.0*k[0]*p[0]*x_tmp[1] - 1.0*k[0]*p[1]*x_tmp[11])/(pow(k[0],2));
xdot_tmp[15] = ((pow(k[0],2))*p[0]*x_tmp[0] + (pow(k[0],2))*p[1]*x_tmp[3] + (pow(k[0],2))*p[1]*x_tmp[4] - 3.0*(pow(k[0],2))*p[0]*x_tmp[15] - 1.0*(pow(k[0],2))*p[1]*x_tmp[15] + (pow(k[0],2))*p[1]*x_tmp[16] - 1.0*k[0]*p[0]*x_tmp[1] - 1.0*k[0]*p[1]*x_tmp[12])/(pow(k[0],2));
xdot_tmp[16] = (2.0*(pow(k[0],2))*p[0]*x_tmp[15] - 2.0*(pow(k[0],2))*p[0]*x_tmp[16] - 2.0*(pow(k[0],2))*p[1]*x_tmp[16] + 2.0*(pow(k[0],2))*p[1]*x_tmp[17] + k[0]*p[0]*x_tmp[1] + k[0]*p[0]*x_tmp[12] + k[0]*p[1]*x_tmp[12] + k[0]*p[1]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[17] = ((pow(k[0],2))*p[0]*x_tmp[3] + (pow(k[0],2))*p[0]*x_tmp[4] + (pow(k[0],2))*p[1]*x_tmp[5] + (pow(k[0],2))*p[0]*x_tmp[16] - 1.0*(pow(k[0],2))*p[0]*x_tmp[17] - 3.0*(pow(k[0],2))*p[1]*x_tmp[17] - 1.0*k[0]*p[0]*x_tmp[12] - 1.0*k[0]*p[1]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[18] = ((pow(k[0],2))*p[0]*x_tmp[2] + (pow(k[0],2))*p[2]*x_tmp[4] + (pow(k[0],2))*p[1]*x_tmp[8] - 1.0*(pow(k[0],2))*p[0]*x_tmp[18] + (pow(k[0],2))*p[2]*x_tmp[16] - 1.0*(pow(k[0],2))*p[1]*x_tmp[18] + 2.0*(pow(k[0],2))*p[2]*x_tmp[17] - 1.0*(pow(k[0],2))*p[3]*x_tmp[18])/(pow(k[0],2));
xdot_tmp[19] = (2.0*(pow(k[0],2))*p[1]*x_tmp[9] + 2.0*(pow(k[0],2))*p[0]*x_tmp[14] - 2.0*(pow(k[0],2))*p[0]*x_tmp[19] - 2.0*(pow(k[0],2))*p[1]*x_tmp[19] + k[0]*p[0]*x_tmp[1] + k[0]*p[0]*x_tmp[11] + k[0]*p[1]*x_tmp[11] + k[0]*p[1]*x_tmp[13])/(pow(k[0],2));

  for (ix=0; ix<20; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_enhancer_10001(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
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
xBdot_tmp[0] = 4.0*p[0]*xB_tmp[0] - 1.0*p[0]*xB_tmp[14] - 1.0*p[0]*xB_tmp[15];
xBdot_tmp[1] = 2.0*p[0]*xB_tmp[1] - 1.0*p[0]*xB_tmp[11] - 1.0*p[0]*xB_tmp[12] - (2.0*p[0]*xB_tmp[0])/k[0] + (p[0]*xB_tmp[14])/k[0] + (p[0]*xB_tmp[15])/k[0] - (1.0*p[0]*xB_tmp[16])/k[0] - (1.0*p[0]*xB_tmp[19])/k[0];
xBdot_tmp[2] = ((2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3])*xB_tmp[2])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[18] - 1.0*p[0]*xB_tmp[10];
xBdot_tmp[3] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[3])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[9] - 1.0*p[1]*xB_tmp[14] - 1.0*p[1]*xB_tmp[15] - 1.0*p[0]*xB_tmp[17] - 2.0*p[2]*xB_tmp[2];
xBdot_tmp[4] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[4])/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[10] - 1.0*p[1]*xB_tmp[14] - 1.0*p[1]*xB_tmp[15] - 1.0*p[0]*xB_tmp[17] - 1.0*p[2]*xB_tmp[18] - 1.0*p[0]*xB_tmp[9];
xBdot_tmp[5] = 4.0*p[1]*xB_tmp[5] - 1.0*p[1]*xB_tmp[9] - 2.0*p[2]*xB_tmp[8] - 1.0*p[1]*xB_tmp[17];
xBdot_tmp[6] = p[3]*xB_tmp[6] - (1.0*p[3]*xB_tmp[7])/k[0];
xBdot_tmp[7] = 2.0*p[3]*xB_tmp[7];
xBdot_tmp[8] = ((2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[8])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[10] - 1.0*p[1]*xB_tmp[18] - 4.0*p[2]*xB_tmp[7];
xBdot_tmp[9] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[9])/(pow(k[0],2)) - 2.0*p[0]*xB_tmp[5] - 1.0*p[1]*xB_tmp[4] - 1.0*p[2]*xB_tmp[8] - 2.0*p[2]*xB_tmp[10] - 2.0*p[1]*xB_tmp[19] - 1.0*p[1]*xB_tmp[3];
xBdot_tmp[10] = (xB_tmp[10]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[8] - 2.0*p[2]*xB_tmp[7] - 1.0*p[1]*xB_tmp[2];
xBdot_tmp[11] = (p[0]*xB_tmp[9])/k[0] - 1.0*p[2]*xB_tmp[6] - 1.0*p[0]*xB_tmp[13] - (1.0*p[1]*xB_tmp[0])/k[0] - (1.0*p[0]*xB_tmp[5])/k[0] - 1.0*p[1]*xB_tmp[1] - (1.0*p[2]*xB_tmp[7])/k[0] + (p[1]*xB_tmp[14])/k[0] + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[11])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[19])/(pow(k[0],2));
xBdot_tmp[12] = (p[1]*xB_tmp[15])/k[0] - 1.0*p[2]*xB_tmp[6] - 1.0*p[0]*xB_tmp[13] - (1.0*p[1]*xB_tmp[0])/k[0] - (1.0*p[0]*xB_tmp[5])/k[0] - (1.0*p[2]*xB_tmp[7])/k[0] - 1.0*p[1]*xB_tmp[1] + (p[0]*xB_tmp[17])/k[0] + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[12])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[16])/(pow(k[0],2));
xBdot_tmp[13] = 2.0*p[1]*xB_tmp[13] - 1.0*p[1]*xB_tmp[11] - 1.0*p[1]*xB_tmp[12] - 2.0*p[2]*xB_tmp[6] - (2.0*p[1]*xB_tmp[5])/k[0] - (2.0*p[2]*xB_tmp[7])/k[0] + (p[1]*xB_tmp[9])/k[0] - (1.0*p[1]*xB_tmp[16])/k[0] + (p[1]*xB_tmp[17])/k[0] - (1.0*p[1]*xB_tmp[19])/k[0];
xBdot_tmp[14] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[14])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[3] - 1.0*p[0]*xB_tmp[4] - 1.0*p[2]*xB_tmp[2] - 2.0*p[0]*xB_tmp[19] - 2.0*p[1]*xB_tmp[0];
xBdot_tmp[15] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[15])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[3] - 1.0*p[0]*xB_tmp[4] - 1.0*p[2]*xB_tmp[2] - 2.0*p[0]*xB_tmp[16] - 2.0*p[1]*xB_tmp[0];
xBdot_tmp[16] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[16])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[17] - 1.0*p[2]*xB_tmp[18] - 1.0*p[1]*xB_tmp[15];
xBdot_tmp[17] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[17])/(pow(k[0],2)) - 2.0*p[0]*xB_tmp[5] - 1.0*p[1]*xB_tmp[4] - 1.0*p[2]*xB_tmp[8] - 2.0*p[1]*xB_tmp[16] - 2.0*p[2]*xB_tmp[18] - 1.0*p[1]*xB_tmp[3];
xBdot_tmp[18] = (xB_tmp[18]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[8] - 2.0*p[2]*xB_tmp[7] - 1.0*p[1]*xB_tmp[2];
xBdot_tmp[19] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[19])/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[10] - 1.0*p[1]*xB_tmp[14] - 1.0*p[0]*xB_tmp[9];
xBdot_tmp[20] = 4.0*p[0]*xB_tmp[20] - 1.0*p[0]*xB_tmp[34] - 1.0*p[0]*xB_tmp[35];
xBdot_tmp[21] = 2.0*p[0]*xB_tmp[21] - 1.0*p[0]*xB_tmp[31] - 1.0*p[0]*xB_tmp[32] - (2.0*p[0]*xB_tmp[20])/k[0] + (p[0]*xB_tmp[34])/k[0] + (p[0]*xB_tmp[35])/k[0] - (1.0*p[0]*xB_tmp[36])/k[0] - (1.0*p[0]*xB_tmp[39])/k[0];
xBdot_tmp[22] = ((2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3])*xB_tmp[22])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[38] - 1.0*p[0]*xB_tmp[30];
xBdot_tmp[23] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[23])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[29] - 1.0*p[1]*xB_tmp[34] - 1.0*p[1]*xB_tmp[35] - 1.0*p[0]*xB_tmp[37] - 2.0*p[2]*xB_tmp[22];
xBdot_tmp[24] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[24])/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[30] - 1.0*p[1]*xB_tmp[34] - 1.0*p[1]*xB_tmp[35] - 1.0*p[0]*xB_tmp[37] - 1.0*p[2]*xB_tmp[38] - 1.0*p[0]*xB_tmp[29];
xBdot_tmp[25] = 4.0*p[1]*xB_tmp[25] - 1.0*p[1]*xB_tmp[29] - 2.0*p[2]*xB_tmp[28] - 1.0*p[1]*xB_tmp[37];
xBdot_tmp[26] = p[3]*xB_tmp[26] - (1.0*p[3]*xB_tmp[27])/k[0];
xBdot_tmp[27] = 2.0*p[3]*xB_tmp[27];
xBdot_tmp[28] = ((2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*xB_tmp[28])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[30] - 1.0*p[1]*xB_tmp[38] - 4.0*p[2]*xB_tmp[27];
xBdot_tmp[29] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[29])/(pow(k[0],2)) - 2.0*p[0]*xB_tmp[25] - 1.0*p[1]*xB_tmp[24] - 1.0*p[2]*xB_tmp[28] - 2.0*p[2]*xB_tmp[30] - 2.0*p[1]*xB_tmp[39] - 1.0*p[1]*xB_tmp[23];
xBdot_tmp[30] = (xB_tmp[30]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[28] - 2.0*p[2]*xB_tmp[27] - 1.0*p[1]*xB_tmp[22];
xBdot_tmp[31] = (p[0]*xB_tmp[29])/k[0] - 1.0*p[2]*xB_tmp[26] - 1.0*p[0]*xB_tmp[33] - (1.0*p[1]*xB_tmp[20])/k[0] - (1.0*p[0]*xB_tmp[25])/k[0] - 1.0*p[1]*xB_tmp[21] - (1.0*p[2]*xB_tmp[27])/k[0] + (p[1]*xB_tmp[34])/k[0] + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[31])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[39])/(pow(k[0],2));
xBdot_tmp[32] = (p[1]*xB_tmp[35])/k[0] - 1.0*p[2]*xB_tmp[26] - 1.0*p[0]*xB_tmp[33] - (1.0*p[1]*xB_tmp[20])/k[0] - (1.0*p[0]*xB_tmp[25])/k[0] - (1.0*p[2]*xB_tmp[27])/k[0] - 1.0*p[1]*xB_tmp[21] + (p[0]*xB_tmp[37])/k[0] + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[32])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[36])/(pow(k[0],2));
xBdot_tmp[33] = 2.0*p[1]*xB_tmp[33] - 1.0*p[1]*xB_tmp[31] - 1.0*p[1]*xB_tmp[32] - 2.0*p[2]*xB_tmp[26] - (2.0*p[1]*xB_tmp[25])/k[0] - (2.0*p[2]*xB_tmp[27])/k[0] + (p[1]*xB_tmp[29])/k[0] - (1.0*p[1]*xB_tmp[36])/k[0] + (p[1]*xB_tmp[37])/k[0] - (1.0*p[1]*xB_tmp[39])/k[0];
xBdot_tmp[34] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[34])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[23] - 1.0*p[0]*xB_tmp[24] - 1.0*p[2]*xB_tmp[22] - 2.0*p[0]*xB_tmp[39] - 2.0*p[1]*xB_tmp[20];
xBdot_tmp[35] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[35])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[23] - 1.0*p[0]*xB_tmp[24] - 1.0*p[2]*xB_tmp[22] - 2.0*p[0]*xB_tmp[36] - 2.0*p[1]*xB_tmp[20];
xBdot_tmp[36] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[36])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[37] - 1.0*p[2]*xB_tmp[38] - 1.0*p[1]*xB_tmp[35];
xBdot_tmp[37] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[37])/(pow(k[0],2)) - 2.0*p[0]*xB_tmp[25] - 1.0*p[1]*xB_tmp[24] - 1.0*p[2]*xB_tmp[28] - 2.0*p[1]*xB_tmp[36] - 2.0*p[2]*xB_tmp[38] - 1.0*p[1]*xB_tmp[23];
xBdot_tmp[38] = (xB_tmp[38]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[28] - 2.0*p[2]*xB_tmp[27] - 1.0*p[1]*xB_tmp[22];
xBdot_tmp[39] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[39])/(pow(k[0],2)) - 1.0*p[2]*xB_tmp[30] - 1.0*p[1]*xB_tmp[34] - 1.0*p[0]*xB_tmp[29];

  for (ixB=0; ixB<40; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_enhancer_10001(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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
qBdot_tmp[0+ip*ny] = 2.0*x_tmp[1]*xB_tmp[1] + 2.0*x_tmp[2]*xB_tmp[2] - (1.0*(k[0]*x_tmp[1] - 1.0*k[0]*x_tmp[11])*xB_tmp[11])/k[0] - (1.0*(k[0]*x_tmp[1] - 1.0*k[0]*x_tmp[12])*xB_tmp[12])/k[0] - (1.0*xB_tmp[5]*(k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[9] + 2.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[16]*(k[0]*x_tmp[1] + k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[15] - 2.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[19]*(k[0]*x_tmp[1] + k[0]*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[14] - 2.0*(pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*((pow(k[0],2))*x_tmp[14] - 2.0*(pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[4]*((pow(k[0],2))*x_tmp[14] - 2.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) + (xB_tmp[14]*(k[0]*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[14]))/(pow(k[0],2)) + (xB_tmp[15]*(k[0]*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[18])*xB_tmp[8])/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[1] - 4.0*(pow(k[0],2))*x_tmp[0])*xB_tmp[0])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[2] - 1.0*(pow(k[0],2))*x_tmp[10])*xB_tmp[10])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[2] - 1.0*(pow(k[0],2))*x_tmp[18])*xB_tmp[18])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[11] + k[0]*x_tmp[12])*xB_tmp[13])/k[0] - (1.0*xB_tmp[9]*((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[17]*((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[12] + (pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[16] - 1.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = 2.0*x_tmp[1]*xB_tmp[21] + 2.0*x_tmp[2]*xB_tmp[22] - (1.0*(k[0]*x_tmp[1] - 1.0*k[0]*x_tmp[11])*xB_tmp[31])/k[0] - (1.0*(k[0]*x_tmp[1] - 1.0*k[0]*x_tmp[12])*xB_tmp[32])/k[0] - (1.0*xB_tmp[25]*(k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[9] + 2.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[36]*(k[0]*x_tmp[1] + k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[15] - 2.0*(pow(k[0],2))*x_tmp[16]))/(pow(k[0],2)) - (1.0*xB_tmp[39]*(k[0]*x_tmp[1] + k[0]*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[14] - 2.0*(pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*((pow(k[0],2))*x_tmp[14] - 2.0*(pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[24]*((pow(k[0],2))*x_tmp[14] - 2.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) + (xB_tmp[34]*(k[0]*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[14]))/(pow(k[0],2)) + (xB_tmp[35]*(k[0]*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[18])*xB_tmp[28])/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[1] - 4.0*(pow(k[0],2))*x_tmp[0])*xB_tmp[20])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[2] - 1.0*(pow(k[0],2))*x_tmp[10])*xB_tmp[30])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[2] - 1.0*(pow(k[0],2))*x_tmp[18])*xB_tmp[38])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[11] + k[0]*x_tmp[12])*xB_tmp[33])/k[0] - (1.0*xB_tmp[29]*((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[37]*((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[12] + (pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[16] - 1.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2));

  } break;

  case 1: {
qBdot_tmp[0+ip*ny] = 2.0*x_tmp[8]*xB_tmp[8] + 2.0*x_tmp[13]*xB_tmp[13] + ((k[0]*x_tmp[11] - 1.0*k[0]*x_tmp[13])*xB_tmp[11])/k[0] + ((k[0]*x_tmp[12] - 1.0*k[0]*x_tmp[13])*xB_tmp[12])/k[0] - (1.0*xB_tmp[0]*(k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[14] + 2.0*(pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[19]*(k[0]*x_tmp[11] + k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[16]*(k[0]*x_tmp[12] + k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[16] + 2.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[3]*((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[4]*((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) + (xB_tmp[9]*(k[0]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) + (xB_tmp[17]*(k[0]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[18])*xB_tmp[2])/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[13] - 4.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[5])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[10])*xB_tmp[10])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[18])*xB_tmp[18])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[11] + k[0]*x_tmp[12])*xB_tmp[1])/k[0] - (1.0*xB_tmp[14]*((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[15]*((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[12] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[15] + (pow(k[0],2))*x_tmp[16]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = 2.0*x_tmp[8]*xB_tmp[28] + 2.0*x_tmp[13]*xB_tmp[33] + ((k[0]*x_tmp[11] - 1.0*k[0]*x_tmp[13])*xB_tmp[31])/k[0] + ((k[0]*x_tmp[12] - 1.0*k[0]*x_tmp[13])*xB_tmp[32])/k[0] - (1.0*xB_tmp[20]*(k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[14] + 2.0*(pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[39]*(k[0]*x_tmp[11] + k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[36]*(k[0]*x_tmp[12] + k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[16] + 2.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[23]*((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[24]*((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) + (xB_tmp[29]*(k[0]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) + (xB_tmp[37]*(k[0]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[18])*xB_tmp[22])/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[13] - 4.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[25])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[10])*xB_tmp[30])/(pow(k[0],2)) - (1.0*((pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[18])*xB_tmp[38])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[11] + k[0]*x_tmp[12])*xB_tmp[21])/k[0] - (1.0*xB_tmp[34]*((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[35]*((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[12] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[15] + (pow(k[0],2))*x_tmp[16]))/(pow(k[0],2));

  } break;

  case 2: {
qBdot_tmp[0+ip*ny] = - (1.0*xB_tmp[6]*(k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*k[0]*x_tmp[13]))/k[0] - (1.0*xB_tmp[2]*(2.0*(pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[8]*(2.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[10]*((pow(k[0],2))*x_tmp[4] + 2.0*(pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[18]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[16] + 2.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[7]*(k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*k[0]*x_tmp[13] + 4.0*(pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[10] + 2.0*(pow(k[0],2))*x_tmp[18]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = - (1.0*xB_tmp[26]*(k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*k[0]*x_tmp[13]))/k[0] - (1.0*xB_tmp[22]*(2.0*(pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*xB_tmp[28]*(2.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[30]*((pow(k[0],2))*x_tmp[4] + 2.0*(pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[19]))/(pow(k[0],2)) - (1.0*xB_tmp[38]*((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[16] + 2.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*xB_tmp[27]*(k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*k[0]*x_tmp[13] + 4.0*(pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[10] + 2.0*(pow(k[0],2))*x_tmp[18]))/(pow(k[0],2));

  } break;

  case 3: {
qBdot_tmp[0+ip*ny] = x_tmp[2]*xB_tmp[2] + x_tmp[6]*xB_tmp[6] + x_tmp[8]*xB_tmp[8] + x_tmp[10]*xB_tmp[10] + x_tmp[18]*xB_tmp[18] - (1.0*(k[0]*x_tmp[6] - 2.0*(pow(k[0],2))*x_tmp[7])*xB_tmp[7])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = x_tmp[2]*xB_tmp[22] + x_tmp[6]*xB_tmp[26] + x_tmp[8]*xB_tmp[28] + x_tmp[10]*xB_tmp[30] + x_tmp[18]*xB_tmp[38] - (1.0*(k[0]*x_tmp[6] - 2.0*(pow(k[0],2))*x_tmp[7])*xB_tmp[27])/(pow(k[0],2));

  } break;

  }
  }

  for (iyp=0; iyp<2*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_enhancer_10001(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*20);
x0_tmp[0] = (k[6]*k[26])/(pow(k[0],2));
x0_tmp[1] = (k[1]*k[21] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[2] = (k[10]*k[30])/(pow(k[0],2));
x0_tmp[3] = (k[9]*k[29])/(pow(k[0],2));
x0_tmp[4] = (k[12]*k[32])/(pow(k[0],2));
x0_tmp[5] = (k[18]*k[38])/(pow(k[0],2));
x0_tmp[6] = (k[5]*k[25])/k[0];
x0_tmp[7] = (k[20]*k[40])/(pow(k[0],2));
x0_tmp[8] = (k[19]*k[39])/(pow(k[0],2));
x0_tmp[9] = (k[13]*k[33])/(pow(k[0],2));
x0_tmp[10] = (k[14]*k[34])/(pow(k[0],2));
x0_tmp[11] = (k[2]*k[22])/k[0];
x0_tmp[12] = (k[3]*k[23])/k[0];
x0_tmp[13] = (k[4]*k[24])/k[0];
x0_tmp[14] = (k[7]*k[27])/(pow(k[0],2));
x0_tmp[15] = (k[8]*k[28])/(pow(k[0],2));
x0_tmp[16] = (k[15]*k[35])/(pow(k[0],2));
x0_tmp[17] = (k[16]*k[36])/(pow(k[0],2));
x0_tmp[18] = (k[17]*k[37])/(pow(k[0],2));
x0_tmp[19] = (k[11]*k[31])/(pow(k[0],2));
  
  
  return;
}


 int Jv_enhancer_10001(N_Vector v, N_Vector Jv, realtype t,
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
Jv_tmp[0] = 2.0*p[1]*v_tmp[14] - 4.0*p[0]*v_tmp[0] + 2.0*p[1]*v_tmp[15] + (2.0*p[0]*v_tmp[1])/k[0] + (p[1]*v_tmp[11])/k[0] + (p[1]*v_tmp[12])/k[0];
Jv_tmp[1] = p[1]*v_tmp[11] - 2.0*p[0]*v_tmp[1] + p[1]*v_tmp[12];
Jv_tmp[2] = 2.0*p[2]*v_tmp[3] + p[1]*v_tmp[10] + p[2]*v_tmp[14] + p[2]*v_tmp[15] + p[1]*v_tmp[18] - (1.0*(2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3])*v_tmp[2])/(pow(k[0],2));
Jv_tmp[3] = p[1]*v_tmp[9] + p[0]*v_tmp[14] + p[0]*v_tmp[15] + p[1]*v_tmp[17] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[3])/(pow(k[0],2));
Jv_tmp[4] = p[1]*v_tmp[9] + p[0]*v_tmp[14] + p[0]*v_tmp[15] + p[1]*v_tmp[17] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[4])/(pow(k[0],2));
Jv_tmp[5] = 2.0*p[0]*v_tmp[9] - 4.0*p[1]*v_tmp[5] + 2.0*p[0]*v_tmp[17] + (p[0]*v_tmp[11])/k[0] + (p[0]*v_tmp[12])/k[0] + (2.0*p[1]*v_tmp[13])/k[0];
Jv_tmp[6] = p[2]*v_tmp[11] - 1.0*p[3]*v_tmp[6] + p[2]*v_tmp[12] + 2.0*p[2]*v_tmp[13];
Jv_tmp[7] = 4.0*p[2]*v_tmp[8] - 2.0*p[3]*v_tmp[7] + 2.0*p[2]*v_tmp[10] + 2.0*p[2]*v_tmp[18] + (p[3]*v_tmp[6])/k[0] + (p[2]*v_tmp[11])/k[0] + (p[2]*v_tmp[12])/k[0] + (2.0*p[2]*v_tmp[13])/k[0];
Jv_tmp[8] = 2.0*p[2]*v_tmp[5] + p[0]*v_tmp[10] + p[2]*v_tmp[9] + p[0]*v_tmp[18] + p[2]*v_tmp[17] - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*v_tmp[8])/(pow(k[0],2));
Jv_tmp[9] = p[0]*v_tmp[3] + p[0]*v_tmp[4] + p[1]*v_tmp[5] + p[0]*v_tmp[19] - (1.0*p[0]*v_tmp[11])/k[0] - (1.0*p[1]*v_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*v_tmp[9])/(pow(k[0],2));
Jv_tmp[10] = p[0]*v_tmp[2] + p[2]*v_tmp[4] + p[1]*v_tmp[8] + 2.0*p[2]*v_tmp[9] + p[2]*v_tmp[19] - (1.0*v_tmp[10]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
Jv_tmp[11] = p[0]*v_tmp[1] + p[1]*v_tmp[13] - (1.0*(k[0]*p[0] + k[0]*p[1])*v_tmp[11])/k[0];
Jv_tmp[12] = p[0]*v_tmp[1] + p[1]*v_tmp[13] - (1.0*(k[0]*p[0] + k[0]*p[1])*v_tmp[12])/k[0];
Jv_tmp[13] = p[0]*v_tmp[11] + p[0]*v_tmp[12] - 2.0*p[1]*v_tmp[13];
Jv_tmp[14] = p[0]*v_tmp[0] + p[1]*v_tmp[3] + p[1]*v_tmp[4] + p[1]*v_tmp[19] - (1.0*p[0]*v_tmp[1])/k[0] - (1.0*p[1]*v_tmp[11])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*v_tmp[14])/(pow(k[0],2));
Jv_tmp[15] = p[0]*v_tmp[0] + p[1]*v_tmp[3] + p[1]*v_tmp[4] + p[1]*v_tmp[16] - (1.0*p[0]*v_tmp[1])/k[0] - (1.0*p[1]*v_tmp[12])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*v_tmp[15])/(pow(k[0],2));
Jv_tmp[16] = 2.0*p[0]*v_tmp[15] + 2.0*p[1]*v_tmp[17] + (p[0]*v_tmp[1])/k[0] + (p[1]*v_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*v_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[16])/(pow(k[0],2));
Jv_tmp[17] = p[0]*v_tmp[3] + p[0]*v_tmp[4] + p[1]*v_tmp[5] + p[0]*v_tmp[16] - (1.0*p[0]*v_tmp[12])/k[0] - (1.0*p[1]*v_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*v_tmp[17])/(pow(k[0],2));
Jv_tmp[18] = p[0]*v_tmp[2] + p[2]*v_tmp[4] + p[1]*v_tmp[8] + p[2]*v_tmp[16] + 2.0*p[2]*v_tmp[17] - (1.0*v_tmp[18]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
Jv_tmp[19] = 2.0*p[1]*v_tmp[9] + 2.0*p[0]*v_tmp[14] + (p[0]*v_tmp[1])/k[0] + (p[1]*v_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*v_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[19])/(pow(k[0],2));

  for (ix=0; ix<20; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_enhancer_10001(N_Vector vB, N_Vector JvB, realtype t,
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
JvB_tmp[0] = 4.0*p[0]*vB_tmp[0] - 1.0*p[0]*vB_tmp[14] - 1.0*p[0]*vB_tmp[15];
JvB_tmp[1] = 2.0*p[0]*vB_tmp[1] - 1.0*p[0]*vB_tmp[11] - 1.0*p[0]*vB_tmp[12] - (2.0*p[0]*vB_tmp[0])/k[0] + (p[0]*vB_tmp[14])/k[0] + (p[0]*vB_tmp[15])/k[0] - (1.0*p[0]*vB_tmp[16])/k[0] - (1.0*p[0]*vB_tmp[19])/k[0];
JvB_tmp[2] = ((2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3])*vB_tmp[2])/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[18] - 1.0*p[0]*vB_tmp[10];
JvB_tmp[3] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[3])/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[9] - 1.0*p[1]*vB_tmp[14] - 1.0*p[1]*vB_tmp[15] - 1.0*p[0]*vB_tmp[17] - 2.0*p[2]*vB_tmp[2];
JvB_tmp[4] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[4])/(pow(k[0],2)) - 1.0*p[2]*vB_tmp[10] - 1.0*p[1]*vB_tmp[14] - 1.0*p[1]*vB_tmp[15] - 1.0*p[0]*vB_tmp[17] - 1.0*p[2]*vB_tmp[18] - 1.0*p[0]*vB_tmp[9];
JvB_tmp[5] = 4.0*p[1]*vB_tmp[5] - 1.0*p[1]*vB_tmp[9] - 2.0*p[2]*vB_tmp[8] - 1.0*p[1]*vB_tmp[17];
JvB_tmp[6] = p[3]*vB_tmp[6] - (1.0*p[3]*vB_tmp[7])/k[0];
JvB_tmp[7] = 2.0*p[3]*vB_tmp[7];
JvB_tmp[8] = ((2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*vB_tmp[8])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[10] - 1.0*p[1]*vB_tmp[18] - 4.0*p[2]*vB_tmp[7];
JvB_tmp[9] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*vB_tmp[9])/(pow(k[0],2)) - 2.0*p[0]*vB_tmp[5] - 1.0*p[1]*vB_tmp[4] - 1.0*p[2]*vB_tmp[8] - 2.0*p[2]*vB_tmp[10] - 2.0*p[1]*vB_tmp[19] - 1.0*p[1]*vB_tmp[3];
JvB_tmp[10] = (vB_tmp[10]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[8] - 2.0*p[2]*vB_tmp[7] - 1.0*p[1]*vB_tmp[2];
JvB_tmp[11] = (p[0]*vB_tmp[9])/k[0] - 1.0*p[2]*vB_tmp[6] - 1.0*p[0]*vB_tmp[13] - (1.0*p[1]*vB_tmp[0])/k[0] - (1.0*p[0]*vB_tmp[5])/k[0] - 1.0*p[1]*vB_tmp[1] - (1.0*p[2]*vB_tmp[7])/k[0] + (p[1]*vB_tmp[14])/k[0] + ((k[0]*p[0] + k[0]*p[1])*vB_tmp[11])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*vB_tmp[19])/(pow(k[0],2));
JvB_tmp[12] = (p[1]*vB_tmp[15])/k[0] - 1.0*p[2]*vB_tmp[6] - 1.0*p[0]*vB_tmp[13] - (1.0*p[1]*vB_tmp[0])/k[0] - (1.0*p[0]*vB_tmp[5])/k[0] - (1.0*p[2]*vB_tmp[7])/k[0] - 1.0*p[1]*vB_tmp[1] + (p[0]*vB_tmp[17])/k[0] + ((k[0]*p[0] + k[0]*p[1])*vB_tmp[12])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*vB_tmp[16])/(pow(k[0],2));
JvB_tmp[13] = 2.0*p[1]*vB_tmp[13] - 1.0*p[1]*vB_tmp[11] - 1.0*p[1]*vB_tmp[12] - 2.0*p[2]*vB_tmp[6] - (2.0*p[1]*vB_tmp[5])/k[0] - (2.0*p[2]*vB_tmp[7])/k[0] + (p[1]*vB_tmp[9])/k[0] - (1.0*p[1]*vB_tmp[16])/k[0] + (p[1]*vB_tmp[17])/k[0] - (1.0*p[1]*vB_tmp[19])/k[0];
JvB_tmp[14] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*vB_tmp[14])/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[3] - 1.0*p[0]*vB_tmp[4] - 1.0*p[2]*vB_tmp[2] - 2.0*p[0]*vB_tmp[19] - 2.0*p[1]*vB_tmp[0];
JvB_tmp[15] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*vB_tmp[15])/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[3] - 1.0*p[0]*vB_tmp[4] - 1.0*p[2]*vB_tmp[2] - 2.0*p[0]*vB_tmp[16] - 2.0*p[1]*vB_tmp[0];
JvB_tmp[16] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[16])/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[17] - 1.0*p[2]*vB_tmp[18] - 1.0*p[1]*vB_tmp[15];
JvB_tmp[17] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*vB_tmp[17])/(pow(k[0],2)) - 2.0*p[0]*vB_tmp[5] - 1.0*p[1]*vB_tmp[4] - 1.0*p[2]*vB_tmp[8] - 2.0*p[1]*vB_tmp[16] - 2.0*p[2]*vB_tmp[18] - 1.0*p[1]*vB_tmp[3];
JvB_tmp[18] = (vB_tmp[18]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[8] - 2.0*p[2]*vB_tmp[7] - 1.0*p[1]*vB_tmp[2];
JvB_tmp[19] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[19])/(pow(k[0],2)) - 1.0*p[2]*vB_tmp[10] - 1.0*p[1]*vB_tmp[14] - 1.0*p[0]*vB_tmp[9];

  for (ix=0; ix<20; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_enhancer_10001(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_enhancer_10001(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_enhancer_10001(long int N, realtype t, N_Vector x,
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
J->data[0] = -4.0*p[0];
J->data[14] = p[0];
J->data[15] = p[0];
J->data[20] = (2.0*p[0])/k[0];
J->data[21] = -2.0*p[0];
J->data[31] = p[0];
J->data[32] = p[0];
J->data[34] = -(1.0*p[0])/k[0];
J->data[35] = -(1.0*p[0])/k[0];
J->data[36] = p[0]/k[0];
J->data[39] = p[0]/k[0];
J->data[42] = -(1.0*(2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[50] = p[0];
J->data[58] = p[0];
J->data[62] = 2.0*p[2];
J->data[63] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[69] = p[0];
J->data[74] = p[1];
J->data[75] = p[1];
J->data[77] = p[0];
J->data[84] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[89] = p[0];
J->data[90] = p[2];
J->data[94] = p[1];
J->data[95] = p[1];
J->data[97] = p[0];
J->data[98] = p[2];
J->data[105] = -4.0*p[1];
J->data[108] = 2.0*p[2];
J->data[109] = p[1];
J->data[117] = p[1];
J->data[126] = -1.0*p[3];
J->data[127] = p[3]/k[0];
J->data[147] = -2.0*p[3];
J->data[167] = 4.0*p[2];
J->data[168] = -(1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[170] = p[1];
J->data[178] = p[1];
J->data[183] = p[1];
J->data[184] = p[1];
J->data[185] = 2.0*p[0];
J->data[188] = p[2];
J->data[189] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[190] = 2.0*p[2];
J->data[199] = 2.0*p[1];
J->data[202] = p[1];
J->data[207] = 2.0*p[2];
J->data[208] = p[0];
J->data[210] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[220] = p[1]/k[0];
J->data[221] = p[1];
J->data[225] = p[0]/k[0];
J->data[226] = p[2];
J->data[227] = p[2]/k[0];
J->data[229] = -(1.0*p[0])/k[0];
J->data[231] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];
J->data[233] = p[0];
J->data[234] = -(1.0*p[1])/k[0];
J->data[239] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[240] = p[1]/k[0];
J->data[241] = p[1];
J->data[245] = p[0]/k[0];
J->data[246] = p[2];
J->data[247] = p[2]/k[0];
J->data[252] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];
J->data[253] = p[0];
J->data[255] = -(1.0*p[1])/k[0];
J->data[256] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[257] = -(1.0*p[0])/k[0];
J->data[265] = (2.0*p[1])/k[0];
J->data[266] = 2.0*p[2];
J->data[267] = (2.0*p[2])/k[0];
J->data[269] = -(1.0*p[1])/k[0];
J->data[271] = p[1];
J->data[272] = p[1];
J->data[273] = -2.0*p[1];
J->data[276] = p[1]/k[0];
J->data[277] = -(1.0*p[1])/k[0];
J->data[279] = p[1]/k[0];
J->data[280] = 2.0*p[1];
J->data[282] = p[2];
J->data[283] = p[0];
J->data[284] = p[0];
J->data[294] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[299] = 2.0*p[0];
J->data[300] = 2.0*p[1];
J->data[302] = p[2];
J->data[303] = p[0];
J->data[304] = p[0];
J->data[315] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[316] = 2.0*p[0];
J->data[335] = p[1];
J->data[336] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[337] = p[0];
J->data[338] = p[2];
J->data[343] = p[1];
J->data[344] = p[1];
J->data[345] = 2.0*p[0];
J->data[348] = p[2];
J->data[356] = 2.0*p[1];
J->data[357] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[358] = 2.0*p[2];
J->data[362] = p[1];
J->data[367] = 2.0*p[2];
J->data[368] = p[0];
J->data[378] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[389] = p[0];
J->data[390] = p[2];
J->data[394] = p[1];
J->data[399] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));

  for (iJ=0; iJ<400; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_enhancer_10001(realtype t, N_Vector x,
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
  J->rowvals[1] = 14;
  J->rowvals[2] = 15;
  J->rowvals[3] = 0;
  J->rowvals[4] = 1;
  J->rowvals[5] = 11;
  J->rowvals[6] = 12;
  J->rowvals[7] = 14;
  J->rowvals[8] = 15;
  J->rowvals[9] = 16;
  J->rowvals[10] = 19;
  J->rowvals[11] = 2;
  J->rowvals[12] = 10;
  J->rowvals[13] = 18;
  J->rowvals[14] = 2;
  J->rowvals[15] = 3;
  J->rowvals[16] = 9;
  J->rowvals[17] = 14;
  J->rowvals[18] = 15;
  J->rowvals[19] = 17;
  J->rowvals[20] = 4;
  J->rowvals[21] = 9;
  J->rowvals[22] = 10;
  J->rowvals[23] = 14;
  J->rowvals[24] = 15;
  J->rowvals[25] = 17;
  J->rowvals[26] = 18;
  J->rowvals[27] = 5;
  J->rowvals[28] = 8;
  J->rowvals[29] = 9;
  J->rowvals[30] = 17;
  J->rowvals[31] = 6;
  J->rowvals[32] = 7;
  J->rowvals[33] = 7;
  J->rowvals[34] = 7;
  J->rowvals[35] = 8;
  J->rowvals[36] = 10;
  J->rowvals[37] = 18;
  J->rowvals[38] = 3;
  J->rowvals[39] = 4;
  J->rowvals[40] = 5;
  J->rowvals[41] = 8;
  J->rowvals[42] = 9;
  J->rowvals[43] = 10;
  J->rowvals[44] = 19;
  J->rowvals[45] = 2;
  J->rowvals[46] = 7;
  J->rowvals[47] = 8;
  J->rowvals[48] = 10;
  J->rowvals[49] = 0;
  J->rowvals[50] = 1;
  J->rowvals[51] = 5;
  J->rowvals[52] = 6;
  J->rowvals[53] = 7;
  J->rowvals[54] = 9;
  J->rowvals[55] = 11;
  J->rowvals[56] = 13;
  J->rowvals[57] = 14;
  J->rowvals[58] = 19;
  J->rowvals[59] = 0;
  J->rowvals[60] = 1;
  J->rowvals[61] = 5;
  J->rowvals[62] = 6;
  J->rowvals[63] = 7;
  J->rowvals[64] = 12;
  J->rowvals[65] = 13;
  J->rowvals[66] = 15;
  J->rowvals[67] = 16;
  J->rowvals[68] = 17;
  J->rowvals[69] = 5;
  J->rowvals[70] = 6;
  J->rowvals[71] = 7;
  J->rowvals[72] = 9;
  J->rowvals[73] = 11;
  J->rowvals[74] = 12;
  J->rowvals[75] = 13;
  J->rowvals[76] = 16;
  J->rowvals[77] = 17;
  J->rowvals[78] = 19;
  J->rowvals[79] = 0;
  J->rowvals[80] = 2;
  J->rowvals[81] = 3;
  J->rowvals[82] = 4;
  J->rowvals[83] = 14;
  J->rowvals[84] = 19;
  J->rowvals[85] = 0;
  J->rowvals[86] = 2;
  J->rowvals[87] = 3;
  J->rowvals[88] = 4;
  J->rowvals[89] = 15;
  J->rowvals[90] = 16;
  J->rowvals[91] = 15;
  J->rowvals[92] = 16;
  J->rowvals[93] = 17;
  J->rowvals[94] = 18;
  J->rowvals[95] = 3;
  J->rowvals[96] = 4;
  J->rowvals[97] = 5;
  J->rowvals[98] = 8;
  J->rowvals[99] = 16;
  J->rowvals[100] = 17;
  J->rowvals[101] = 18;
  J->rowvals[102] = 2;
  J->rowvals[103] = 7;
  J->rowvals[104] = 8;
  J->rowvals[105] = 18;
  J->rowvals[106] = 9;
  J->rowvals[107] = 10;
  J->rowvals[108] = 14;
  J->rowvals[109] = 19;
  J->colptrs[0] = 0;
  J->colptrs[1] = 3;
  J->colptrs[2] = 11;
  J->colptrs[3] = 14;
  J->colptrs[4] = 20;
  J->colptrs[5] = 27;
  J->colptrs[6] = 31;
  J->colptrs[7] = 33;
  J->colptrs[8] = 34;
  J->colptrs[9] = 38;
  J->colptrs[10] = 45;
  J->colptrs[11] = 49;
  J->colptrs[12] = 59;
  J->colptrs[13] = 69;
  J->colptrs[14] = 79;
  J->colptrs[15] = 85;
  J->colptrs[16] = 91;
  J->colptrs[17] = 95;
  J->colptrs[18] = 102;
  J->colptrs[19] = 106;
  J->colptrs[20] = 110;
J->data[0] = -4.0*p[0];
J->data[1] = p[0];
J->data[2] = p[0];
J->data[3] = (2.0*p[0])/k[0];
J->data[4] = -2.0*p[0];
J->data[5] = p[0];
J->data[6] = p[0];
J->data[7] = -(1.0*p[0])/k[0];
J->data[8] = -(1.0*p[0])/k[0];
J->data[9] = p[0]/k[0];
J->data[10] = p[0]/k[0];
J->data[11] = -(1.0*(2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[12] = p[0];
J->data[13] = p[0];
J->data[14] = 2.0*p[2];
J->data[15] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[16] = p[0];
J->data[17] = p[1];
J->data[18] = p[1];
J->data[19] = p[0];
J->data[20] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[21] = p[0];
J->data[22] = p[2];
J->data[23] = p[1];
J->data[24] = p[1];
J->data[25] = p[0];
J->data[26] = p[2];
J->data[27] = -4.0*p[1];
J->data[28] = 2.0*p[2];
J->data[29] = p[1];
J->data[30] = p[1];
J->data[31] = -1.0*p[3];
J->data[32] = p[3]/k[0];
J->data[33] = -2.0*p[3];
J->data[34] = 4.0*p[2];
J->data[35] = -(1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[36] = p[1];
J->data[37] = p[1];
J->data[38] = p[1];
J->data[39] = p[1];
J->data[40] = 2.0*p[0];
J->data[41] = p[2];
J->data[42] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[43] = 2.0*p[2];
J->data[44] = 2.0*p[1];
J->data[45] = p[1];
J->data[46] = 2.0*p[2];
J->data[47] = p[0];
J->data[48] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[49] = p[1]/k[0];
J->data[50] = p[1];
J->data[51] = p[0]/k[0];
J->data[52] = p[2];
J->data[53] = p[2]/k[0];
J->data[54] = -(1.0*p[0])/k[0];
J->data[55] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];
J->data[56] = p[0];
J->data[57] = -(1.0*p[1])/k[0];
J->data[58] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[59] = p[1]/k[0];
J->data[60] = p[1];
J->data[61] = p[0]/k[0];
J->data[62] = p[2];
J->data[63] = p[2]/k[0];
J->data[64] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];
J->data[65] = p[0];
J->data[66] = -(1.0*p[1])/k[0];
J->data[67] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[68] = -(1.0*p[0])/k[0];
J->data[69] = (2.0*p[1])/k[0];
J->data[70] = 2.0*p[2];
J->data[71] = (2.0*p[2])/k[0];
J->data[72] = -(1.0*p[1])/k[0];
J->data[73] = p[1];
J->data[74] = p[1];
J->data[75] = -2.0*p[1];
J->data[76] = p[1]/k[0];
J->data[77] = -(1.0*p[1])/k[0];
J->data[78] = p[1]/k[0];
J->data[79] = 2.0*p[1];
J->data[80] = p[2];
J->data[81] = p[0];
J->data[82] = p[0];
J->data[83] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[84] = 2.0*p[0];
J->data[85] = 2.0*p[1];
J->data[86] = p[2];
J->data[87] = p[0];
J->data[88] = p[0];
J->data[89] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[90] = 2.0*p[0];
J->data[91] = p[1];
J->data[92] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[93] = p[0];
J->data[94] = p[2];
J->data[95] = p[1];
J->data[96] = p[1];
J->data[97] = 2.0*p[0];
J->data[98] = p[2];
J->data[99] = 2.0*p[1];
J->data[100] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[101] = 2.0*p[2];
J->data[102] = p[1];
J->data[103] = 2.0*p[2];
J->data[104] = p[0];
J->data[105] = -(1.0*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
J->data[106] = p[0];
J->data[107] = p[2];
J->data[108] = p[1];
J->data[109] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
  return(0);
}


 int JBBand_enhancer_10001(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_enhancer_10001(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_enhancer_10001(long int N, realtype t, N_Vector x,
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
JB->data[0] = 4.0*p[0];
JB->data[1] = -(2.0*p[0])/k[0];
JB->data[11] = -(1.0*p[1])/k[0];
JB->data[12] = -(1.0*p[1])/k[0];
JB->data[14] = -2.0*p[1];
JB->data[15] = -2.0*p[1];
JB->data[21] = 2.0*p[0];
JB->data[31] = -1.0*p[1];
JB->data[32] = -1.0*p[1];
JB->data[42] = (2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3])/(pow(k[0],2));
JB->data[43] = -2.0*p[2];
JB->data[50] = -1.0*p[1];
JB->data[54] = -1.0*p[2];
JB->data[55] = -1.0*p[2];
JB->data[58] = -1.0*p[1];
JB->data[63] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[69] = -1.0*p[1];
JB->data[74] = -1.0*p[0];
JB->data[75] = -1.0*p[0];
JB->data[77] = -1.0*p[1];
JB->data[84] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[89] = -1.0*p[1];
JB->data[94] = -1.0*p[0];
JB->data[95] = -1.0*p[0];
JB->data[97] = -1.0*p[1];
JB->data[105] = 4.0*p[1];
JB->data[109] = -2.0*p[0];
JB->data[111] = -(1.0*p[0])/k[0];
JB->data[112] = -(1.0*p[0])/k[0];
JB->data[113] = -(2.0*p[1])/k[0];
JB->data[117] = -2.0*p[0];
JB->data[126] = p[3];
JB->data[131] = -1.0*p[2];
JB->data[132] = -1.0*p[2];
JB->data[133] = -2.0*p[2];
JB->data[146] = -(1.0*p[3])/k[0];
JB->data[147] = 2.0*p[3];
JB->data[148] = -4.0*p[2];
JB->data[150] = -2.0*p[2];
JB->data[151] = -(1.0*p[2])/k[0];
JB->data[152] = -(1.0*p[2])/k[0];
JB->data[153] = -(2.0*p[2])/k[0];
JB->data[158] = -2.0*p[2];
JB->data[165] = -2.0*p[2];
JB->data[168] = (2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])/(pow(k[0],2));
JB->data[169] = -1.0*p[2];
JB->data[170] = -1.0*p[0];
JB->data[177] = -1.0*p[2];
JB->data[178] = -1.0*p[0];
JB->data[183] = -1.0*p[0];
JB->data[184] = -1.0*p[0];
JB->data[185] = -1.0*p[1];
JB->data[189] = ((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[191] = p[0]/k[0];
JB->data[193] = p[1]/k[0];
JB->data[199] = -1.0*p[0];
JB->data[202] = -1.0*p[0];
JB->data[204] = -1.0*p[2];
JB->data[208] = -1.0*p[1];
JB->data[209] = -2.0*p[2];
JB->data[210] = ((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])/(pow(k[0],2));
JB->data[219] = -1.0*p[2];
JB->data[221] = -1.0*p[0];
JB->data[231] = (k[0]*p[0] + k[0]*p[1])/k[0];
JB->data[233] = -1.0*p[1];
JB->data[241] = -1.0*p[0];
JB->data[252] = (k[0]*p[0] + k[0]*p[1])/k[0];
JB->data[253] = -1.0*p[1];
JB->data[271] = -1.0*p[0];
JB->data[272] = -1.0*p[0];
JB->data[273] = 2.0*p[1];
JB->data[280] = -1.0*p[0];
JB->data[281] = p[0]/k[0];
JB->data[283] = -1.0*p[1];
JB->data[284] = -1.0*p[1];
JB->data[291] = p[1]/k[0];
JB->data[294] = (3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[299] = -1.0*p[1];
JB->data[300] = -1.0*p[0];
JB->data[301] = p[0]/k[0];
JB->data[303] = -1.0*p[1];
JB->data[304] = -1.0*p[1];
JB->data[312] = p[1]/k[0];
JB->data[315] = (3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[316] = -1.0*p[1];
JB->data[321] = -(1.0*p[0])/k[0];
JB->data[332] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/(pow(k[0],2));
JB->data[333] = -(1.0*p[1])/k[0];
JB->data[335] = -2.0*p[0];
JB->data[336] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[337] = -2.0*p[1];
JB->data[343] = -1.0*p[0];
JB->data[344] = -1.0*p[0];
JB->data[345] = -1.0*p[1];
JB->data[352] = p[0]/k[0];
JB->data[353] = p[1]/k[0];
JB->data[356] = -1.0*p[0];
JB->data[357] = ((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[362] = -1.0*p[0];
JB->data[364] = -1.0*p[2];
JB->data[368] = -1.0*p[1];
JB->data[376] = -1.0*p[2];
JB->data[377] = -2.0*p[2];
JB->data[378] = ((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])/(pow(k[0],2));
JB->data[381] = -(1.0*p[0])/k[0];
JB->data[389] = -2.0*p[1];
JB->data[391] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/(pow(k[0],2));
JB->data[393] = -(1.0*p[1])/k[0];
JB->data[394] = -2.0*p[0];
JB->data[399] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));

  for (iJ=0; iJ<400; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_enhancer_10001(realtype t, N_Vector x,
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
  JB->rowvals[4] = 14;
  JB->rowvals[5] = 15;
  JB->rowvals[6] = 1;
  JB->rowvals[7] = 11;
  JB->rowvals[8] = 12;
  JB->rowvals[9] = 2;
  JB->rowvals[10] = 3;
  JB->rowvals[11] = 10;
  JB->rowvals[12] = 14;
  JB->rowvals[13] = 15;
  JB->rowvals[14] = 18;
  JB->rowvals[15] = 3;
  JB->rowvals[16] = 9;
  JB->rowvals[17] = 14;
  JB->rowvals[18] = 15;
  JB->rowvals[19] = 17;
  JB->rowvals[20] = 4;
  JB->rowvals[21] = 9;
  JB->rowvals[22] = 14;
  JB->rowvals[23] = 15;
  JB->rowvals[24] = 17;
  JB->rowvals[25] = 5;
  JB->rowvals[26] = 9;
  JB->rowvals[27] = 11;
  JB->rowvals[28] = 12;
  JB->rowvals[29] = 13;
  JB->rowvals[30] = 17;
  JB->rowvals[31] = 6;
  JB->rowvals[32] = 11;
  JB->rowvals[33] = 12;
  JB->rowvals[34] = 13;
  JB->rowvals[35] = 6;
  JB->rowvals[36] = 7;
  JB->rowvals[37] = 8;
  JB->rowvals[38] = 10;
  JB->rowvals[39] = 11;
  JB->rowvals[40] = 12;
  JB->rowvals[41] = 13;
  JB->rowvals[42] = 18;
  JB->rowvals[43] = 5;
  JB->rowvals[44] = 8;
  JB->rowvals[45] = 9;
  JB->rowvals[46] = 10;
  JB->rowvals[47] = 17;
  JB->rowvals[48] = 18;
  JB->rowvals[49] = 3;
  JB->rowvals[50] = 4;
  JB->rowvals[51] = 5;
  JB->rowvals[52] = 9;
  JB->rowvals[53] = 11;
  JB->rowvals[54] = 13;
  JB->rowvals[55] = 19;
  JB->rowvals[56] = 2;
  JB->rowvals[57] = 4;
  JB->rowvals[58] = 8;
  JB->rowvals[59] = 9;
  JB->rowvals[60] = 10;
  JB->rowvals[61] = 19;
  JB->rowvals[62] = 1;
  JB->rowvals[63] = 11;
  JB->rowvals[64] = 13;
  JB->rowvals[65] = 1;
  JB->rowvals[66] = 12;
  JB->rowvals[67] = 13;
  JB->rowvals[68] = 11;
  JB->rowvals[69] = 12;
  JB->rowvals[70] = 13;
  JB->rowvals[71] = 0;
  JB->rowvals[72] = 1;
  JB->rowvals[73] = 3;
  JB->rowvals[74] = 4;
  JB->rowvals[75] = 11;
  JB->rowvals[76] = 14;
  JB->rowvals[77] = 19;
  JB->rowvals[78] = 0;
  JB->rowvals[79] = 1;
  JB->rowvals[80] = 3;
  JB->rowvals[81] = 4;
  JB->rowvals[82] = 12;
  JB->rowvals[83] = 15;
  JB->rowvals[84] = 16;
  JB->rowvals[85] = 1;
  JB->rowvals[86] = 12;
  JB->rowvals[87] = 13;
  JB->rowvals[88] = 15;
  JB->rowvals[89] = 16;
  JB->rowvals[90] = 17;
  JB->rowvals[91] = 3;
  JB->rowvals[92] = 4;
  JB->rowvals[93] = 5;
  JB->rowvals[94] = 12;
  JB->rowvals[95] = 13;
  JB->rowvals[96] = 16;
  JB->rowvals[97] = 17;
  JB->rowvals[98] = 2;
  JB->rowvals[99] = 4;
  JB->rowvals[100] = 8;
  JB->rowvals[101] = 16;
  JB->rowvals[102] = 17;
  JB->rowvals[103] = 18;
  JB->rowvals[104] = 1;
  JB->rowvals[105] = 9;
  JB->rowvals[106] = 11;
  JB->rowvals[107] = 13;
  JB->rowvals[108] = 14;
  JB->rowvals[109] = 19;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 6;
  JB->colptrs[2] = 9;
  JB->colptrs[3] = 15;
  JB->colptrs[4] = 20;
  JB->colptrs[5] = 25;
  JB->colptrs[6] = 31;
  JB->colptrs[7] = 35;
  JB->colptrs[8] = 43;
  JB->colptrs[9] = 49;
  JB->colptrs[10] = 56;
  JB->colptrs[11] = 62;
  JB->colptrs[12] = 65;
  JB->colptrs[13] = 68;
  JB->colptrs[14] = 71;
  JB->colptrs[15] = 78;
  JB->colptrs[16] = 85;
  JB->colptrs[17] = 91;
  JB->colptrs[18] = 98;
  JB->colptrs[19] = 104;
  JB->colptrs[20] = 110;
  return(0);
}


 int sx_enhancer_10001(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[14] - 4.0*p[0]*sx_tmp[0] + 2.0*p[1]*sx_tmp[15] + (2.0*k[0]*x_tmp[1] - 4.0*(pow(k[0],2))*x_tmp[0])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[11])/k[0] + (p[1]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[1]*sx_tmp[11] - 2.0*p[0]*sx_tmp[1] + p[1]*sx_tmp[12] - 2.0*x_tmp[1];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[3] + p[1]*sx_tmp[10] + p[2]*sx_tmp[14] + p[2]*sx_tmp[15] + p[1]*sx_tmp[18] - 2.0*x_tmp[2] - (1.0*(2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3])*sx_tmp[2])/(pow(k[0],2));
sxdot_tmp[3] = p[1]*sx_tmp[9] + p[0]*sx_tmp[14] + p[0]*sx_tmp[15] + p[1]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[14] - 2.0*(pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[15])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[0]*sx_tmp[14] + p[0]*sx_tmp[15] + p[1]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[14] - 2.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[15])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = 2.0*p[0]*sx_tmp[9] - 4.0*p[1]*sx_tmp[5] + 2.0*p[0]*sx_tmp[17] + (k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[9] + 2.0*(pow(k[0],2))*x_tmp[17])/(pow(k[0],2)) + (p[0]*sx_tmp[11])/k[0] + (p[0]*sx_tmp[12])/k[0] + (2.0*p[1]*sx_tmp[13])/k[0];
sxdot_tmp[6] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[6] + p[2]*sx_tmp[12] + 2.0*p[2]*sx_tmp[13];
sxdot_tmp[7] = 4.0*p[2]*sx_tmp[8] - 2.0*p[3]*sx_tmp[7] + 2.0*p[2]*sx_tmp[10] + 2.0*p[2]*sx_tmp[18] + (p[3]*sx_tmp[6])/k[0] + (p[2]*sx_tmp[11])/k[0] + (p[2]*sx_tmp[12])/k[0] + (2.0*p[2]*sx_tmp[13])/k[0];
sxdot_tmp[8] = 2.0*p[2]*sx_tmp[5] + p[0]*sx_tmp[10] + p[2]*sx_tmp[9] + p[0]*sx_tmp[18] + p[2]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[18])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[8])/(pow(k[0],2));
sxdot_tmp[9] = p[0]*sx_tmp[3] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] + p[0]*sx_tmp[19] + ((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[19])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[11])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[9])/(pow(k[0],2));
sxdot_tmp[10] = p[0]*sx_tmp[2] + p[2]*sx_tmp[4] + p[1]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + p[2]*sx_tmp[19] + ((pow(k[0],2))*x_tmp[2] - 1.0*(pow(k[0],2))*x_tmp[10])/(pow(k[0],2)) - (1.0*sx_tmp[10]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
sxdot_tmp[11] = p[0]*sx_tmp[1] + p[1]*sx_tmp[13] + (k[0]*x_tmp[1] - 1.0*k[0]*x_tmp[11])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[11])/k[0];
sxdot_tmp[12] = p[0]*sx_tmp[1] + p[1]*sx_tmp[13] + (k[0]*x_tmp[1] - 1.0*k[0]*x_tmp[12])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[12])/k[0];
sxdot_tmp[13] = p[0]*sx_tmp[11] + p[0]*sx_tmp[12] - 2.0*p[1]*sx_tmp[13] + (k[0]*x_tmp[11] + k[0]*x_tmp[12])/k[0];
sxdot_tmp[14] = p[0]*sx_tmp[0] + p[1]*sx_tmp[3] + p[1]*sx_tmp[4] + p[1]*sx_tmp[19] - (1.0*(k[0]*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[14]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1])/k[0] - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[14])/(pow(k[0],2));
sxdot_tmp[15] = p[0]*sx_tmp[0] + p[1]*sx_tmp[3] + p[1]*sx_tmp[4] + p[1]*sx_tmp[16] - (1.0*(k[0]*x_tmp[1] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[15]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1])/k[0] - (1.0*p[1]*sx_tmp[12])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[15])/(pow(k[0],2));
sxdot_tmp[16] = 2.0*p[0]*sx_tmp[15] + 2.0*p[1]*sx_tmp[17] + (k[0]*x_tmp[1] + k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[15] - 2.0*(pow(k[0],2))*x_tmp[16])/(pow(k[0],2)) + (p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[16])/(pow(k[0],2));
sxdot_tmp[17] = p[0]*sx_tmp[3] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] + p[0]*sx_tmp[16] + ((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[12] + (pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[16] - 1.0*(pow(k[0],2))*x_tmp[17])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[12])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[17])/(pow(k[0],2));
sxdot_tmp[18] = p[0]*sx_tmp[2] + p[2]*sx_tmp[4] + p[1]*sx_tmp[8] + p[2]*sx_tmp[16] + 2.0*p[2]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[2] - 1.0*(pow(k[0],2))*x_tmp[18])/(pow(k[0],2)) - (1.0*sx_tmp[18]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
sxdot_tmp[19] = 2.0*p[1]*sx_tmp[9] + 2.0*p[0]*sx_tmp[14] + (k[0]*x_tmp[1] + k[0]*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[14] - 2.0*(pow(k[0],2))*x_tmp[19])/(pow(k[0],2)) + (p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[19])/(pow(k[0],2));

  } break;

  case 1: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[14] - 4.0*p[0]*sx_tmp[0] + 2.0*p[1]*sx_tmp[15] + (k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*(pow(k[0],2))*x_tmp[14] + 2.0*(pow(k[0],2))*x_tmp[15])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[11])/k[0] + (p[1]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[1]*sx_tmp[11] - 2.0*p[0]*sx_tmp[1] + p[1]*sx_tmp[12] + (k[0]*x_tmp[11] + k[0]*x_tmp[12])/k[0];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[3] + p[1]*sx_tmp[10] + p[2]*sx_tmp[14] + p[2]*sx_tmp[15] + p[1]*sx_tmp[18] + ((pow(k[0],2))*x_tmp[10] + (pow(k[0],2))*x_tmp[18])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3])*sx_tmp[2])/(pow(k[0],2));
sxdot_tmp[3] = p[1]*sx_tmp[9] + p[0]*sx_tmp[14] + p[0]*sx_tmp[15] + p[1]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[17])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[0]*sx_tmp[14] + p[0]*sx_tmp[15] + p[1]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[17])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = 2.0*p[0]*sx_tmp[9] - 4.0*p[1]*sx_tmp[5] + 2.0*p[0]*sx_tmp[17] + (2.0*k[0]*x_tmp[13] - 4.0*(pow(k[0],2))*x_tmp[5])/(pow(k[0],2)) + (p[0]*sx_tmp[11])/k[0] + (p[0]*sx_tmp[12])/k[0] + (2.0*p[1]*sx_tmp[13])/k[0];
sxdot_tmp[6] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[6] + p[2]*sx_tmp[12] + 2.0*p[2]*sx_tmp[13];
sxdot_tmp[7] = 4.0*p[2]*sx_tmp[8] - 2.0*p[3]*sx_tmp[7] + 2.0*p[2]*sx_tmp[10] + 2.0*p[2]*sx_tmp[18] + (p[3]*sx_tmp[6])/k[0] + (p[2]*sx_tmp[11])/k[0] + (p[2]*sx_tmp[12])/k[0] + (2.0*p[2]*sx_tmp[13])/k[0];
sxdot_tmp[8] = 2.0*p[2]*sx_tmp[5] + p[0]*sx_tmp[10] + p[2]*sx_tmp[9] + p[0]*sx_tmp[18] + p[2]*sx_tmp[17] - 2.0*x_tmp[8] - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[8])/(pow(k[0],2));
sxdot_tmp[9] = p[0]*sx_tmp[3] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] + p[0]*sx_tmp[19] - (1.0*(k[0]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[11])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[9])/(pow(k[0],2));
sxdot_tmp[10] = p[0]*sx_tmp[2] + p[2]*sx_tmp[4] + p[1]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + p[2]*sx_tmp[19] + ((pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[10])/(pow(k[0],2)) - (1.0*sx_tmp[10]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
sxdot_tmp[11] = p[0]*sx_tmp[1] + p[1]*sx_tmp[13] - (1.0*(k[0]*x_tmp[11] - 1.0*k[0]*x_tmp[13]))/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[11])/k[0];
sxdot_tmp[12] = p[0]*sx_tmp[1] + p[1]*sx_tmp[13] - (1.0*(k[0]*x_tmp[12] - 1.0*k[0]*x_tmp[13]))/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[12])/k[0];
sxdot_tmp[13] = p[0]*sx_tmp[11] + p[0]*sx_tmp[12] - 2.0*p[1]*sx_tmp[13] - 2.0*x_tmp[13];
sxdot_tmp[14] = p[0]*sx_tmp[0] + p[1]*sx_tmp[3] + p[1]*sx_tmp[4] + p[1]*sx_tmp[19] + ((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[11] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[19])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1])/k[0] - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[14])/(pow(k[0],2));
sxdot_tmp[15] = p[0]*sx_tmp[0] + p[1]*sx_tmp[3] + p[1]*sx_tmp[4] + p[1]*sx_tmp[16] + ((pow(k[0],2))*x_tmp[3] - 1.0*k[0]*x_tmp[12] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[15] + (pow(k[0],2))*x_tmp[16])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[1])/k[0] - (1.0*p[1]*sx_tmp[12])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[15])/(pow(k[0],2));
sxdot_tmp[16] = 2.0*p[0]*sx_tmp[15] + 2.0*p[1]*sx_tmp[17] + (k[0]*x_tmp[12] + k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[16] + 2.0*(pow(k[0],2))*x_tmp[17])/(pow(k[0],2)) + (p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[16])/(pow(k[0],2));
sxdot_tmp[17] = p[0]*sx_tmp[3] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] + p[0]*sx_tmp[16] - (1.0*(k[0]*x_tmp[13] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[17]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[12])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[17])/(pow(k[0],2));
sxdot_tmp[18] = p[0]*sx_tmp[2] + p[2]*sx_tmp[4] + p[1]*sx_tmp[8] + p[2]*sx_tmp[16] + 2.0*p[2]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[18])/(pow(k[0],2)) - (1.0*sx_tmp[18]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
sxdot_tmp[19] = 2.0*p[1]*sx_tmp[9] + 2.0*p[0]*sx_tmp[14] + (k[0]*x_tmp[11] + k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[19])/(pow(k[0],2)) + (p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[19])/(pow(k[0],2));

  } break;

  case 2: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[14] - 4.0*p[0]*sx_tmp[0] + 2.0*p[1]*sx_tmp[15] + (2.0*p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[11])/k[0] + (p[1]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[1]*sx_tmp[11] - 2.0*p[0]*sx_tmp[1] + p[1]*sx_tmp[12];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[3] + p[1]*sx_tmp[10] + p[2]*sx_tmp[14] + p[2]*sx_tmp[15] + p[1]*sx_tmp[18] + (2.0*(pow(k[0],2))*x_tmp[3] + (pow(k[0],2))*x_tmp[14] + (pow(k[0],2))*x_tmp[15])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3])*sx_tmp[2])/(pow(k[0],2));
sxdot_tmp[3] = p[1]*sx_tmp[9] + p[0]*sx_tmp[14] + p[0]*sx_tmp[15] + p[1]*sx_tmp[17] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[0]*sx_tmp[14] + p[0]*sx_tmp[15] + p[1]*sx_tmp[17] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = 2.0*p[0]*sx_tmp[9] - 4.0*p[1]*sx_tmp[5] + 2.0*p[0]*sx_tmp[17] + (p[0]*sx_tmp[11])/k[0] + (p[0]*sx_tmp[12])/k[0] + (2.0*p[1]*sx_tmp[13])/k[0];
sxdot_tmp[6] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[6] + p[2]*sx_tmp[12] + 2.0*p[2]*sx_tmp[13] + (k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*k[0]*x_tmp[13])/k[0];
sxdot_tmp[7] = 4.0*p[2]*sx_tmp[8] - 2.0*p[3]*sx_tmp[7] + 2.0*p[2]*sx_tmp[10] + 2.0*p[2]*sx_tmp[18] + (k[0]*x_tmp[11] + k[0]*x_tmp[12] + 2.0*k[0]*x_tmp[13] + 4.0*(pow(k[0],2))*x_tmp[8] + 2.0*(pow(k[0],2))*x_tmp[10] + 2.0*(pow(k[0],2))*x_tmp[18])/(pow(k[0],2)) + (p[3]*sx_tmp[6])/k[0] + (p[2]*sx_tmp[11])/k[0] + (p[2]*sx_tmp[12])/k[0] + (2.0*p[2]*sx_tmp[13])/k[0];
sxdot_tmp[8] = 2.0*p[2]*sx_tmp[5] + p[0]*sx_tmp[10] + p[2]*sx_tmp[9] + p[0]*sx_tmp[18] + p[2]*sx_tmp[17] + (2.0*(pow(k[0],2))*x_tmp[5] + (pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[17])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[8])/(pow(k[0],2));
sxdot_tmp[9] = p[0]*sx_tmp[3] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] + p[0]*sx_tmp[19] - (1.0*p[0]*sx_tmp[11])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[9])/(pow(k[0],2));
sxdot_tmp[10] = p[0]*sx_tmp[2] + p[2]*sx_tmp[4] + p[1]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + p[2]*sx_tmp[19] + ((pow(k[0],2))*x_tmp[4] + 2.0*(pow(k[0],2))*x_tmp[9] + (pow(k[0],2))*x_tmp[19])/(pow(k[0],2)) - (1.0*sx_tmp[10]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
sxdot_tmp[11] = p[0]*sx_tmp[1] + p[1]*sx_tmp[13] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[11])/k[0];
sxdot_tmp[12] = p[0]*sx_tmp[1] + p[1]*sx_tmp[13] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[12])/k[0];
sxdot_tmp[13] = p[0]*sx_tmp[11] + p[0]*sx_tmp[12] - 2.0*p[1]*sx_tmp[13];
sxdot_tmp[14] = p[0]*sx_tmp[0] + p[1]*sx_tmp[3] + p[1]*sx_tmp[4] + p[1]*sx_tmp[19] - (1.0*p[0]*sx_tmp[1])/k[0] - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[14])/(pow(k[0],2));
sxdot_tmp[15] = p[0]*sx_tmp[0] + p[1]*sx_tmp[3] + p[1]*sx_tmp[4] + p[1]*sx_tmp[16] - (1.0*p[0]*sx_tmp[1])/k[0] - (1.0*p[1]*sx_tmp[12])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[15])/(pow(k[0],2));
sxdot_tmp[16] = 2.0*p[0]*sx_tmp[15] + 2.0*p[1]*sx_tmp[17] + (p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[16])/(pow(k[0],2));
sxdot_tmp[17] = p[0]*sx_tmp[3] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] + p[0]*sx_tmp[16] - (1.0*p[0]*sx_tmp[12])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[17])/(pow(k[0],2));
sxdot_tmp[18] = p[0]*sx_tmp[2] + p[2]*sx_tmp[4] + p[1]*sx_tmp[8] + p[2]*sx_tmp[16] + 2.0*p[2]*sx_tmp[17] + ((pow(k[0],2))*x_tmp[4] + (pow(k[0],2))*x_tmp[16] + 2.0*(pow(k[0],2))*x_tmp[17])/(pow(k[0],2)) - (1.0*sx_tmp[18]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
sxdot_tmp[19] = 2.0*p[1]*sx_tmp[9] + 2.0*p[0]*sx_tmp[14] + (p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[19])/(pow(k[0],2));

  } break;

  case 3: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[14] - 4.0*p[0]*sx_tmp[0] + 2.0*p[1]*sx_tmp[15] + (2.0*p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[11])/k[0] + (p[1]*sx_tmp[12])/k[0];
sxdot_tmp[1] = p[1]*sx_tmp[11] - 2.0*p[0]*sx_tmp[1] + p[1]*sx_tmp[12];
sxdot_tmp[2] = 2.0*p[2]*sx_tmp[3] + p[1]*sx_tmp[10] + p[2]*sx_tmp[14] + p[2]*sx_tmp[15] + p[1]*sx_tmp[18] - 1.0*x_tmp[2] - (1.0*(2.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[3])*sx_tmp[2])/(pow(k[0],2));
sxdot_tmp[3] = p[1]*sx_tmp[9] + p[0]*sx_tmp[14] + p[0]*sx_tmp[15] + p[1]*sx_tmp[17] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[9] + p[0]*sx_tmp[14] + p[0]*sx_tmp[15] + p[1]*sx_tmp[17] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = 2.0*p[0]*sx_tmp[9] - 4.0*p[1]*sx_tmp[5] + 2.0*p[0]*sx_tmp[17] + (p[0]*sx_tmp[11])/k[0] + (p[0]*sx_tmp[12])/k[0] + (2.0*p[1]*sx_tmp[13])/k[0];
sxdot_tmp[6] = p[2]*sx_tmp[11] - 1.0*p[3]*sx_tmp[6] + p[2]*sx_tmp[12] + 2.0*p[2]*sx_tmp[13] - 1.0*x_tmp[6];
sxdot_tmp[7] = 4.0*p[2]*sx_tmp[8] - 2.0*p[3]*sx_tmp[7] + 2.0*p[2]*sx_tmp[10] + 2.0*p[2]*sx_tmp[18] + (k[0]*x_tmp[6] - 2.0*(pow(k[0],2))*x_tmp[7])/(pow(k[0],2)) + (p[3]*sx_tmp[6])/k[0] + (p[2]*sx_tmp[11])/k[0] + (p[2]*sx_tmp[12])/k[0] + (2.0*p[2]*sx_tmp[13])/k[0];
sxdot_tmp[8] = 2.0*p[2]*sx_tmp[5] + p[0]*sx_tmp[10] + p[2]*sx_tmp[9] + p[0]*sx_tmp[18] + p[2]*sx_tmp[17] - 1.0*x_tmp[8] - (1.0*(2.0*(pow(k[0],2))*p[1] + (pow(k[0],2))*p[3])*sx_tmp[8])/(pow(k[0],2));
sxdot_tmp[9] = p[0]*sx_tmp[3] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] + p[0]*sx_tmp[19] - (1.0*p[0]*sx_tmp[11])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[9])/(pow(k[0],2));
sxdot_tmp[10] = p[0]*sx_tmp[2] + p[2]*sx_tmp[4] + p[1]*sx_tmp[8] + 2.0*p[2]*sx_tmp[9] + p[2]*sx_tmp[19] - 1.0*x_tmp[10] - (1.0*sx_tmp[10]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
sxdot_tmp[11] = p[0]*sx_tmp[1] + p[1]*sx_tmp[13] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[11])/k[0];
sxdot_tmp[12] = p[0]*sx_tmp[1] + p[1]*sx_tmp[13] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[12])/k[0];
sxdot_tmp[13] = p[0]*sx_tmp[11] + p[0]*sx_tmp[12] - 2.0*p[1]*sx_tmp[13];
sxdot_tmp[14] = p[0]*sx_tmp[0] + p[1]*sx_tmp[3] + p[1]*sx_tmp[4] + p[1]*sx_tmp[19] - (1.0*p[0]*sx_tmp[1])/k[0] - (1.0*p[1]*sx_tmp[11])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[14])/(pow(k[0],2));
sxdot_tmp[15] = p[0]*sx_tmp[0] + p[1]*sx_tmp[3] + p[1]*sx_tmp[4] + p[1]*sx_tmp[16] - (1.0*p[0]*sx_tmp[1])/k[0] - (1.0*p[1]*sx_tmp[12])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[15])/(pow(k[0],2));
sxdot_tmp[16] = 2.0*p[0]*sx_tmp[15] + 2.0*p[1]*sx_tmp[17] + (p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[16])/(pow(k[0],2));
sxdot_tmp[17] = p[0]*sx_tmp[3] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] + p[0]*sx_tmp[16] - (1.0*p[0]*sx_tmp[12])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[17])/(pow(k[0],2));
sxdot_tmp[18] = p[0]*sx_tmp[2] + p[2]*sx_tmp[4] + p[1]*sx_tmp[8] + p[2]*sx_tmp[16] + 2.0*p[2]*sx_tmp[17] - 1.0*x_tmp[18] - (1.0*sx_tmp[18]*((pow(k[0],2))*p[0] + (pow(k[0],2))*p[1] + (pow(k[0],2))*p[3]))/(pow(k[0],2));
sxdot_tmp[19] = 2.0*p[1]*sx_tmp[9] + 2.0*p[0]*sx_tmp[14] + (p[0]*sx_tmp[1])/k[0] + (p[1]*sx_tmp[13])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[11])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[19])/(pow(k[0],2));

  } break;

  }
 for (ix=0; ix<20; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_enhancer_10001(int ip, N_Vector sx0, void *user_data)
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


void y_enhancer_10001(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*6];
y[it+nt*1] = x[it+nt*7];
    
    return;
}


void dydp_enhancer_10001(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }
  
  return;
}


void dydx_enhancer_10001(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*40);
dydx[12] = 1.0;
dydx[15] = 1.0;
  
  return;
}


void sy_enhancer_10001(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];

  } break;

  case 2: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];

  } break;

  case 3: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(7+np*nx)];

  } break;

  }
  
  return;
}
int root_enhancer_10001(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_enhancer_10001(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root_enhancer_10001(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval_enhancer_10001(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval_enhancer_10001(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
void deltadisc_enhancer_10001(double t, int idisc, N_Vector x, void *user_data){
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
void sdeltadisc_enhancer_10001(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data){
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


void dxdotdp_enhancer_10001(double t, int nt, int it, double *dxdotdp, double *p, double *k, double *u, double *x, int *plist, int np, int nx){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
dxdotdp[(0+ip*nx)] = (2.0*k[0]*x[it+nt*1] - 4.0*(pow(k[0],2))*x[it+nt*0])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = -2.0*x[it+nt*1];
dxdotdp[(2+ip*nx)] = -2.0*x[it+nt*2];
dxdotdp[(3+ip*nx)] = ((pow(k[0],2))*x[it+nt*14] - 2.0*(pow(k[0],2))*x[it+nt*3] + (pow(k[0],2))*x[it+nt*15])/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = ((pow(k[0],2))*x[it+nt*14] - 2.0*(pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*15])/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = (k[0]*x[it+nt*11] + k[0]*x[it+nt*12] + 2.0*(pow(k[0],2))*x[it+nt*9] + 2.0*(pow(k[0],2))*x[it+nt*17])/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = ((pow(k[0],2))*x[it+nt*10] + (pow(k[0],2))*x[it+nt*18])/(pow(k[0],2));
dxdotdp[(9+ip*nx)] = ((pow(k[0],2))*x[it+nt*3] - 1.0*k[0]*x[it+nt*11] + (pow(k[0],2))*x[it+nt*4] - 1.0*(pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*19])/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = ((pow(k[0],2))*x[it+nt*2] - 1.0*(pow(k[0],2))*x[it+nt*10])/(pow(k[0],2));
dxdotdp[(11+ip*nx)] = (k[0]*x[it+nt*1] - 1.0*k[0]*x[it+nt*11])/k[0];
dxdotdp[(12+ip*nx)] = (k[0]*x[it+nt*1] - 1.0*k[0]*x[it+nt*12])/k[0];
dxdotdp[(13+ip*nx)] = (k[0]*x[it+nt*11] + k[0]*x[it+nt*12])/k[0];
dxdotdp[(14+ip*nx)] = -(1.0*(k[0]*x[it+nt*1] - 1.0*(pow(k[0],2))*x[it+nt*0] + 3.0*(pow(k[0],2))*x[it+nt*14]))/(pow(k[0],2));
dxdotdp[(15+ip*nx)] = -(1.0*(k[0]*x[it+nt*1] - 1.0*(pow(k[0],2))*x[it+nt*0] + 3.0*(pow(k[0],2))*x[it+nt*15]))/(pow(k[0],2));
dxdotdp[(16+ip*nx)] = (k[0]*x[it+nt*1] + k[0]*x[it+nt*12] + 2.0*(pow(k[0],2))*x[it+nt*15] - 2.0*(pow(k[0],2))*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(17+ip*nx)] = ((pow(k[0],2))*x[it+nt*3] - 1.0*k[0]*x[it+nt*12] + (pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*16] - 1.0*(pow(k[0],2))*x[it+nt*17])/(pow(k[0],2));
dxdotdp[(18+ip*nx)] = ((pow(k[0],2))*x[it+nt*2] - 1.0*(pow(k[0],2))*x[it+nt*18])/(pow(k[0],2));
dxdotdp[(19+ip*nx)] = (k[0]*x[it+nt*1] + k[0]*x[it+nt*11] + 2.0*(pow(k[0],2))*x[it+nt*14] - 2.0*(pow(k[0],2))*x[it+nt*19])/(pow(k[0],2));

  } break;

  case 1: {
dxdotdp[(0+ip*nx)] = (k[0]*x[it+nt*11] + k[0]*x[it+nt*12] + 2.0*(pow(k[0],2))*x[it+nt*14] + 2.0*(pow(k[0],2))*x[it+nt*15])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = (k[0]*x[it+nt*11] + k[0]*x[it+nt*12])/k[0];
dxdotdp[(2+ip*nx)] = ((pow(k[0],2))*x[it+nt*10] + (pow(k[0],2))*x[it+nt*18])/(pow(k[0],2));
dxdotdp[(3+ip*nx)] = ((pow(k[0],2))*x[it+nt*9] - 2.0*(pow(k[0],2))*x[it+nt*3] + (pow(k[0],2))*x[it+nt*17])/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = ((pow(k[0],2))*x[it+nt*9] - 2.0*(pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*17])/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = (2.0*k[0]*x[it+nt*13] - 4.0*(pow(k[0],2))*x[it+nt*5])/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = -2.0*x[it+nt*8];
dxdotdp[(9+ip*nx)] = -(1.0*(k[0]*x[it+nt*13] - 1.0*(pow(k[0],2))*x[it+nt*5] + 3.0*(pow(k[0],2))*x[it+nt*9]))/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = ((pow(k[0],2))*x[it+nt*8] - 1.0*(pow(k[0],2))*x[it+nt*10])/(pow(k[0],2));
dxdotdp[(11+ip*nx)] = -(1.0*(k[0]*x[it+nt*11] - 1.0*k[0]*x[it+nt*13]))/k[0];
dxdotdp[(12+ip*nx)] = -(1.0*(k[0]*x[it+nt*12] - 1.0*k[0]*x[it+nt*13]))/k[0];
dxdotdp[(13+ip*nx)] = -2.0*x[it+nt*13];
dxdotdp[(14+ip*nx)] = ((pow(k[0],2))*x[it+nt*3] - 1.0*k[0]*x[it+nt*11] + (pow(k[0],2))*x[it+nt*4] - 1.0*(pow(k[0],2))*x[it+nt*14] + (pow(k[0],2))*x[it+nt*19])/(pow(k[0],2));
dxdotdp[(15+ip*nx)] = ((pow(k[0],2))*x[it+nt*3] - 1.0*k[0]*x[it+nt*12] + (pow(k[0],2))*x[it+nt*4] - 1.0*(pow(k[0],2))*x[it+nt*15] + (pow(k[0],2))*x[it+nt*16])/(pow(k[0],2));
dxdotdp[(16+ip*nx)] = (k[0]*x[it+nt*12] + k[0]*x[it+nt*13] - 2.0*(pow(k[0],2))*x[it+nt*16] + 2.0*(pow(k[0],2))*x[it+nt*17])/(pow(k[0],2));
dxdotdp[(17+ip*nx)] = -(1.0*(k[0]*x[it+nt*13] - 1.0*(pow(k[0],2))*x[it+nt*5] + 3.0*(pow(k[0],2))*x[it+nt*17]))/(pow(k[0],2));
dxdotdp[(18+ip*nx)] = ((pow(k[0],2))*x[it+nt*8] - 1.0*(pow(k[0],2))*x[it+nt*18])/(pow(k[0],2));
dxdotdp[(19+ip*nx)] = (k[0]*x[it+nt*11] + k[0]*x[it+nt*13] + 2.0*(pow(k[0],2))*x[it+nt*9] - 2.0*(pow(k[0],2))*x[it+nt*19])/(pow(k[0],2));

  } break;

  case 2: {
dxdotdp[(2+ip*nx)] = (2.0*(pow(k[0],2))*x[it+nt*3] + (pow(k[0],2))*x[it+nt*14] + (pow(k[0],2))*x[it+nt*15])/(pow(k[0],2));
dxdotdp[(6+ip*nx)] = (k[0]*x[it+nt*11] + k[0]*x[it+nt*12] + 2.0*k[0]*x[it+nt*13])/k[0];
dxdotdp[(7+ip*nx)] = (k[0]*x[it+nt*11] + k[0]*x[it+nt*12] + 2.0*k[0]*x[it+nt*13] + 4.0*(pow(k[0],2))*x[it+nt*8] + 2.0*(pow(k[0],2))*x[it+nt*10] + 2.0*(pow(k[0],2))*x[it+nt*18])/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = (2.0*(pow(k[0],2))*x[it+nt*5] + (pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*17])/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = ((pow(k[0],2))*x[it+nt*4] + 2.0*(pow(k[0],2))*x[it+nt*9] + (pow(k[0],2))*x[it+nt*19])/(pow(k[0],2));
dxdotdp[(18+ip*nx)] = ((pow(k[0],2))*x[it+nt*4] + (pow(k[0],2))*x[it+nt*16] + 2.0*(pow(k[0],2))*x[it+nt*17])/(pow(k[0],2));

  } break;

  case 3: {
dxdotdp[(2+ip*nx)] = -1.0*x[it+nt*2];
dxdotdp[(6+ip*nx)] = -1.0*x[it+nt*6];
dxdotdp[(7+ip*nx)] = (k[0]*x[it+nt*6] - 2.0*(pow(k[0],2))*x[it+nt*7])/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = -1.0*x[it+nt*8];
dxdotdp[(10+ip*nx)] = -1.0*x[it+nt*10];
dxdotdp[(18+ip*nx)] = -1.0*x[it+nt*18];

  } break;

  }
  }
  
  return;
}
