#include "enhancer_667.h"
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


 int xdot_enhancer_667(realtype t, N_Vector x, N_Vector xdot, void *user_data)
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
xdot_tmp[0] = (2.0*(pow(k[0],2))*p[1]*x_tmp[9] - 4.0*(pow(k[0],2))*p[0]*x_tmp[0] + 2.0*(pow(k[0],2))*p[1]*x_tmp[10] + 2.0*k[0]*p[0]*x_tmp[6] + k[0]*p[1]*x_tmp[7] + k[0]*p[1]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[1] = ((pow(k[0],2))*p[0]*x_tmp[9] - 2.0*(pow(k[0],2))*p[1]*x_tmp[1] - 2.0*(pow(k[0],2))*p[0]*x_tmp[1] + (pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[1]*x_tmp[11] + (pow(k[0],2))*p[1]*x_tmp[12])/(pow(k[0],2));
xdot_tmp[2] = ((pow(k[0],2))*p[0]*x_tmp[9] - 2.0*(pow(k[0],2))*p[1]*x_tmp[2] - 2.0*(pow(k[0],2))*p[0]*x_tmp[2] + (pow(k[0],2))*p[0]*x_tmp[10] + (pow(k[0],2))*p[1]*x_tmp[11] + (pow(k[0],2))*p[1]*x_tmp[12])/(pow(k[0],2));
xdot_tmp[3] = (2.0*(pow(k[0],2))*p[0]*x_tmp[9] - 2.0*(pow(k[0],2))*p[1]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[3] + 2.0*(pow(k[0],2))*p[1]*x_tmp[11] + k[0]*p[0]*x_tmp[6] + k[0]*p[0]*x_tmp[7] + k[0]*p[1]*x_tmp[7] + k[0]*p[1]*x_tmp[8])/(pow(k[0],2));
xdot_tmp[4] = (2.0*(pow(k[0],2))*p[0]*x_tmp[10] - 2.0*(pow(k[0],2))*p[1]*x_tmp[4] - 2.0*(pow(k[0],2))*p[0]*x_tmp[4] + 2.0*(pow(k[0],2))*p[1]*x_tmp[12] + k[0]*p[0]*x_tmp[6] + k[0]*p[1]*x_tmp[8] + k[0]*p[0]*x_tmp[13] + k[0]*p[1]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[5] = (2.0*(pow(k[0],2))*p[0]*x_tmp[11] - 4.0*(pow(k[0],2))*p[1]*x_tmp[5] + 2.0*(pow(k[0],2))*p[0]*x_tmp[12] + k[0]*p[0]*x_tmp[7] + 2.0*k[0]*p[1]*x_tmp[8] + k[0]*p[0]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[6] = (k[0]*p[1]*x_tmp[7] - 2.0*k[0]*p[0]*x_tmp[6] + k[0]*p[1]*x_tmp[13])/k[0];
xdot_tmp[7] = (k[0]*p[0]*x_tmp[6] - 1.0*k[0]*p[0]*x_tmp[7] - 1.0*k[0]*p[1]*x_tmp[7] + k[0]*p[1]*x_tmp[8])/k[0];
xdot_tmp[8] = (k[0]*p[0]*x_tmp[7] - 2.0*k[0]*p[1]*x_tmp[8] + k[0]*p[0]*x_tmp[13])/k[0];
xdot_tmp[9] = ((pow(k[0],2))*p[0]*x_tmp[0] + (pow(k[0],2))*p[1]*x_tmp[1] + (pow(k[0],2))*p[1]*x_tmp[2] + (pow(k[0],2))*p[1]*x_tmp[3] - 3.0*(pow(k[0],2))*p[0]*x_tmp[9] - 1.0*(pow(k[0],2))*p[1]*x_tmp[9] - 1.0*k[0]*p[0]*x_tmp[6] - 1.0*k[0]*p[1]*x_tmp[7])/(pow(k[0],2));
xdot_tmp[10] = ((pow(k[0],2))*p[0]*x_tmp[0] + (pow(k[0],2))*p[1]*x_tmp[1] + (pow(k[0],2))*p[1]*x_tmp[2] + (pow(k[0],2))*p[1]*x_tmp[4] - 3.0*(pow(k[0],2))*p[0]*x_tmp[10] - 1.0*(pow(k[0],2))*p[1]*x_tmp[10] - 1.0*k[0]*p[0]*x_tmp[6] - 1.0*k[0]*p[1]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[11] = ((pow(k[0],2))*p[0]*x_tmp[1] + (pow(k[0],2))*p[0]*x_tmp[2] + (pow(k[0],2))*p[0]*x_tmp[3] + (pow(k[0],2))*p[1]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[11] - 3.0*(pow(k[0],2))*p[1]*x_tmp[11] - 1.0*k[0]*p[0]*x_tmp[7] - 1.0*k[0]*p[1]*x_tmp[8])/(pow(k[0],2));
xdot_tmp[12] = ((pow(k[0],2))*p[0]*x_tmp[1] + (pow(k[0],2))*p[0]*x_tmp[2] + (pow(k[0],2))*p[0]*x_tmp[4] + (pow(k[0],2))*p[1]*x_tmp[5] - 1.0*(pow(k[0],2))*p[0]*x_tmp[12] - 3.0*(pow(k[0],2))*p[1]*x_tmp[12] - 1.0*k[0]*p[1]*x_tmp[8] - 1.0*k[0]*p[0]*x_tmp[13])/(pow(k[0],2));
xdot_tmp[13] = (k[0]*p[0]*x_tmp[6] + k[0]*p[1]*x_tmp[8] - 1.0*k[0]*p[0]*x_tmp[13] - 1.0*k[0]*p[1]*x_tmp[13])/k[0];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_enhancer_667(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
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
xBdot_tmp[0] = 4.0*p[0]*xB_tmp[0] - 1.0*p[0]*xB_tmp[9] - 1.0*p[0]*xB_tmp[10];
xBdot_tmp[1] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[1])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[11] - 1.0*p[1]*xB_tmp[10] - 1.0*p[0]*xB_tmp[12] - 1.0*p[1]*xB_tmp[9];
xBdot_tmp[2] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[2])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[11] - 1.0*p[1]*xB_tmp[10] - 1.0*p[0]*xB_tmp[12] - 1.0*p[1]*xB_tmp[9];
xBdot_tmp[3] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[3])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[11] - 1.0*p[1]*xB_tmp[9];
xBdot_tmp[4] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[4])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[12] - 1.0*p[1]*xB_tmp[10];
xBdot_tmp[5] = 4.0*p[1]*xB_tmp[5] - 1.0*p[1]*xB_tmp[11] - 1.0*p[1]*xB_tmp[12];
xBdot_tmp[6] = 2.0*p[0]*xB_tmp[6] - 1.0*p[0]*xB_tmp[7] - 1.0*p[0]*xB_tmp[13] - (2.0*p[0]*xB_tmp[0])/k[0] - (1.0*p[0]*xB_tmp[3])/k[0] - (1.0*p[0]*xB_tmp[4])/k[0] + (p[0]*xB_tmp[9])/k[0] + (p[0]*xB_tmp[10])/k[0];
xBdot_tmp[7] = (p[1]*xB_tmp[9])/k[0] - 1.0*p[0]*xB_tmp[8] - (1.0*p[1]*xB_tmp[0])/k[0] - (1.0*p[0]*xB_tmp[5])/k[0] - 1.0*p[1]*xB_tmp[6] + (p[0]*xB_tmp[11])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[3])/(pow(k[0],2)) + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[7])/k[0];
xBdot_tmp[8] = 2.0*p[1]*xB_tmp[8] - 1.0*p[1]*xB_tmp[7] - 1.0*p[1]*xB_tmp[13] - (1.0*p[1]*xB_tmp[3])/k[0] - (1.0*p[1]*xB_tmp[4])/k[0] - (2.0*p[1]*xB_tmp[5])/k[0] + (p[1]*xB_tmp[11])/k[0] + (p[1]*xB_tmp[12])/k[0];
xBdot_tmp[9] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[9])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[0] - 1.0*p[0]*xB_tmp[2] - 2.0*p[0]*xB_tmp[3] - 1.0*p[0]*xB_tmp[1];
xBdot_tmp[10] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[10])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[0] - 1.0*p[0]*xB_tmp[2] - 2.0*p[0]*xB_tmp[4] - 1.0*p[0]*xB_tmp[1];
xBdot_tmp[11] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[11])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[2] - 2.0*p[1]*xB_tmp[3] - 2.0*p[0]*xB_tmp[5] - 1.0*p[1]*xB_tmp[1];
xBdot_tmp[12] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[12])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[2] - 2.0*p[0]*xB_tmp[5] - 2.0*p[1]*xB_tmp[4] - 1.0*p[1]*xB_tmp[1];
xBdot_tmp[13] = (p[1]*xB_tmp[10])/k[0] - 1.0*p[0]*xB_tmp[8] - (1.0*p[1]*xB_tmp[0])/k[0] - (1.0*p[0]*xB_tmp[5])/k[0] - 1.0*p[1]*xB_tmp[6] + (p[0]*xB_tmp[12])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[4])/(pow(k[0],2)) + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[13])/k[0];
xBdot_tmp[14] = 4.0*p[0]*xB_tmp[14] - 1.0*p[0]*xB_tmp[23] - 1.0*p[0]*xB_tmp[24];
xBdot_tmp[15] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[15])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[25] - 1.0*p[1]*xB_tmp[24] - 1.0*p[0]*xB_tmp[26] - 1.0*p[1]*xB_tmp[23];
xBdot_tmp[16] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[16])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[25] - 1.0*p[1]*xB_tmp[24] - 1.0*p[0]*xB_tmp[26] - 1.0*p[1]*xB_tmp[23];
xBdot_tmp[17] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[17])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[25] - 1.0*p[1]*xB_tmp[23];
xBdot_tmp[18] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[18])/(pow(k[0],2)) - 1.0*p[0]*xB_tmp[26] - 1.0*p[1]*xB_tmp[24];
xBdot_tmp[19] = 4.0*p[1]*xB_tmp[19] - 1.0*p[1]*xB_tmp[25] - 1.0*p[1]*xB_tmp[26];
xBdot_tmp[20] = 2.0*p[0]*xB_tmp[20] - 1.0*p[0]*xB_tmp[21] - 1.0*p[0]*xB_tmp[27] - (2.0*p[0]*xB_tmp[14])/k[0] - (1.0*p[0]*xB_tmp[17])/k[0] - (1.0*p[0]*xB_tmp[18])/k[0] + (p[0]*xB_tmp[23])/k[0] + (p[0]*xB_tmp[24])/k[0];
xBdot_tmp[21] = (p[1]*xB_tmp[23])/k[0] - 1.0*p[0]*xB_tmp[22] - (1.0*p[1]*xB_tmp[14])/k[0] - (1.0*p[0]*xB_tmp[19])/k[0] - 1.0*p[1]*xB_tmp[20] + (p[0]*xB_tmp[25])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[17])/(pow(k[0],2)) + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[21])/k[0];
xBdot_tmp[22] = 2.0*p[1]*xB_tmp[22] - 1.0*p[1]*xB_tmp[21] - 1.0*p[1]*xB_tmp[27] - (1.0*p[1]*xB_tmp[17])/k[0] - (1.0*p[1]*xB_tmp[18])/k[0] - (2.0*p[1]*xB_tmp[19])/k[0] + (p[1]*xB_tmp[25])/k[0] + (p[1]*xB_tmp[26])/k[0];
xBdot_tmp[23] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[23])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[14] - 1.0*p[0]*xB_tmp[16] - 2.0*p[0]*xB_tmp[17] - 1.0*p[0]*xB_tmp[15];
xBdot_tmp[24] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[24])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[14] - 1.0*p[0]*xB_tmp[16] - 2.0*p[0]*xB_tmp[18] - 1.0*p[0]*xB_tmp[15];
xBdot_tmp[25] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[25])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[16] - 2.0*p[1]*xB_tmp[17] - 2.0*p[0]*xB_tmp[19] - 1.0*p[1]*xB_tmp[15];
xBdot_tmp[26] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[26])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[16] - 2.0*p[0]*xB_tmp[19] - 2.0*p[1]*xB_tmp[18] - 1.0*p[1]*xB_tmp[15];
xBdot_tmp[27] = (p[1]*xB_tmp[24])/k[0] - 1.0*p[0]*xB_tmp[22] - (1.0*p[1]*xB_tmp[14])/k[0] - (1.0*p[0]*xB_tmp[19])/k[0] - 1.0*p[1]*xB_tmp[20] + (p[0]*xB_tmp[26])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[18])/(pow(k[0],2)) + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[27])/k[0];

  for (ixB=0; ixB<28; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_enhancer_667(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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
qBdot_tmp[0+ip*ny] = 2.0*x_tmp[6]*xB_tmp[6] - (1.0*(k[0]*x_tmp[6] - 1.0*k[0]*x_tmp[7])*xB_tmp[7])/k[0] - (1.0*(k[0]*x_tmp[6] - 1.0*k[0]*x_tmp[13])*xB_tmp[13])/k[0] - (1.0*xB_tmp[3]*(k[0]*x_tmp[6] + k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[3] + 2.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[4]*(k[0]*x_tmp[6] + k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[4] + 2.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[5]*(k[0]*x_tmp[7] + k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*xB_tmp[1]*((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[2]*((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) + (xB_tmp[9]*(k[0]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) + (xB_tmp[10]*(k[0]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[6] - 4.0*(pow(k[0],2))*x_tmp[0])*xB_tmp[0])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[7] + k[0]*x_tmp[13])*xB_tmp[8])/k[0] - (1.0*xB_tmp[11]*((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) - (1.0*xB_tmp[12]*((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[12]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = 2.0*x_tmp[6]*xB_tmp[20] - (1.0*(k[0]*x_tmp[6] - 1.0*k[0]*x_tmp[7])*xB_tmp[21])/k[0] - (1.0*(k[0]*x_tmp[6] - 1.0*k[0]*x_tmp[13])*xB_tmp[27])/k[0] - (1.0*xB_tmp[17]*(k[0]*x_tmp[6] + k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[3] + 2.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[18]*(k[0]*x_tmp[6] + k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[4] + 2.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[19]*(k[0]*x_tmp[7] + k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*xB_tmp[15]*((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[16]*((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) + (xB_tmp[23]*(k[0]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) + (xB_tmp[24]*(k[0]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[6] - 4.0*(pow(k[0],2))*x_tmp[0])*xB_tmp[14])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[7] + k[0]*x_tmp[13])*xB_tmp[22])/k[0] - (1.0*xB_tmp[25]*((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) - (1.0*xB_tmp[26]*((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[12]))/(pow(k[0],2));

  } break;

  case 1: {
qBdot_tmp[0+ip*ny] = 2.0*x_tmp[8]*xB_tmp[8] + ((k[0]*x_tmp[7] - 1.0*k[0]*x_tmp[8])*xB_tmp[7])/k[0] - (1.0*(k[0]*x_tmp[8] - 1.0*k[0]*x_tmp[13])*xB_tmp[13])/k[0] - (1.0*xB_tmp[3]*(k[0]*x_tmp[7] + k[0]*x_tmp[8] - 2.0*(pow(k[0],2))*x_tmp[3] + 2.0*(pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) - (1.0*xB_tmp[0]*(k[0]*x_tmp[7] + k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[9] + 2.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[4]*(k[0]*x_tmp[8] + k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[4] + 2.0*(pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*xB_tmp[1]*((pow(k[0],2))*x_tmp[11] - 2.0*(pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*xB_tmp[2]*((pow(k[0],2))*x_tmp[11] - 2.0*(pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) + (xB_tmp[11]*(k[0]*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) + (xB_tmp[12]*(k[0]*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[8] - 4.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[5])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[7] + k[0]*x_tmp[13])*xB_tmp[6])/k[0] - (1.0*xB_tmp[9]*((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[10]*((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = 2.0*x_tmp[8]*xB_tmp[22] + ((k[0]*x_tmp[7] - 1.0*k[0]*x_tmp[8])*xB_tmp[21])/k[0] - (1.0*(k[0]*x_tmp[8] - 1.0*k[0]*x_tmp[13])*xB_tmp[27])/k[0] - (1.0*xB_tmp[17]*(k[0]*x_tmp[7] + k[0]*x_tmp[8] - 2.0*(pow(k[0],2))*x_tmp[3] + 2.0*(pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) - (1.0*xB_tmp[14]*(k[0]*x_tmp[7] + k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[9] + 2.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*xB_tmp[18]*(k[0]*x_tmp[8] + k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[4] + 2.0*(pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*xB_tmp[15]*((pow(k[0],2))*x_tmp[11] - 2.0*(pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*xB_tmp[16]*((pow(k[0],2))*x_tmp[11] - 2.0*(pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) + (xB_tmp[25]*(k[0]*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) + (xB_tmp[26]*(k[0]*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[8] - 4.0*(pow(k[0],2))*x_tmp[5])*xB_tmp[19])/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[7] + k[0]*x_tmp[13])*xB_tmp[20])/k[0] - (1.0*xB_tmp[23]*((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) - (1.0*xB_tmp[24]*((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2));

  } break;

  }
  }

  for (iyp=0; iyp<2*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_enhancer_667(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*14);
x0_tmp[0] = (k[5]*k[19])/(pow(k[0],2));
x0_tmp[1] = (k[10]*k[24])/(pow(k[0],2));
x0_tmp[2] = (k[8]*k[22])/(pow(k[0],2));
x0_tmp[3] = (k[9]*k[23])/(pow(k[0],2));
x0_tmp[4] = (k[12]*k[26])/(pow(k[0],2));
x0_tmp[5] = (k[14]*k[28])/(pow(k[0],2));
x0_tmp[6] = (k[1]*k[15] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[7] = (k[2]*k[16])/k[0];
x0_tmp[8] = (k[4]*k[18])/k[0];
x0_tmp[9] = (k[6]*k[20])/(pow(k[0],2));
x0_tmp[10] = (k[7]*k[21])/(pow(k[0],2));
x0_tmp[11] = (k[11]*k[25])/(pow(k[0],2));
x0_tmp[12] = (k[13]*k[27])/(pow(k[0],2));
x0_tmp[13] = (k[3]*k[17])/k[0];
  
  
  return;
}


 int Jv_enhancer_667(N_Vector v, N_Vector Jv, realtype t,
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
Jv_tmp[0] = 2.0*p[1]*v_tmp[9] - 4.0*p[0]*v_tmp[0] + 2.0*p[1]*v_tmp[10] + (2.0*p[0]*v_tmp[6])/k[0] + (p[1]*v_tmp[7])/k[0] + (p[1]*v_tmp[13])/k[0];
Jv_tmp[1] = p[0]*v_tmp[9] + p[0]*v_tmp[10] + p[1]*v_tmp[11] + p[1]*v_tmp[12] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[1])/(pow(k[0],2));
Jv_tmp[2] = p[0]*v_tmp[9] + p[0]*v_tmp[10] + p[1]*v_tmp[11] + p[1]*v_tmp[12] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[2])/(pow(k[0],2));
Jv_tmp[3] = 2.0*p[0]*v_tmp[9] + 2.0*p[1]*v_tmp[11] + (p[0]*v_tmp[6])/k[0] + (p[1]*v_tmp[8])/k[0] + ((k[0]*p[0] + k[0]*p[1])*v_tmp[7])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[3])/(pow(k[0],2));
Jv_tmp[4] = 2.0*p[0]*v_tmp[10] + 2.0*p[1]*v_tmp[12] + (p[0]*v_tmp[6])/k[0] + (p[1]*v_tmp[8])/k[0] + ((k[0]*p[0] + k[0]*p[1])*v_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[4])/(pow(k[0],2));
Jv_tmp[5] = 2.0*p[0]*v_tmp[11] - 4.0*p[1]*v_tmp[5] + 2.0*p[0]*v_tmp[12] + (p[0]*v_tmp[7])/k[0] + (2.0*p[1]*v_tmp[8])/k[0] + (p[0]*v_tmp[13])/k[0];
Jv_tmp[6] = p[1]*v_tmp[7] - 2.0*p[0]*v_tmp[6] + p[1]*v_tmp[13];
Jv_tmp[7] = p[0]*v_tmp[6] + p[1]*v_tmp[8] - (1.0*(k[0]*p[0] + k[0]*p[1])*v_tmp[7])/k[0];
Jv_tmp[8] = p[0]*v_tmp[7] - 2.0*p[1]*v_tmp[8] + p[0]*v_tmp[13];
Jv_tmp[9] = p[0]*v_tmp[0] + p[1]*v_tmp[1] + p[1]*v_tmp[2] + p[1]*v_tmp[3] - (1.0*p[0]*v_tmp[6])/k[0] - (1.0*p[1]*v_tmp[7])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*v_tmp[9])/(pow(k[0],2));
Jv_tmp[10] = p[0]*v_tmp[0] + p[1]*v_tmp[1] + p[1]*v_tmp[2] + p[1]*v_tmp[4] - (1.0*p[0]*v_tmp[6])/k[0] - (1.0*p[1]*v_tmp[13])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*v_tmp[10])/(pow(k[0],2));
Jv_tmp[11] = p[0]*v_tmp[1] + p[0]*v_tmp[2] + p[0]*v_tmp[3] + p[1]*v_tmp[5] - (1.0*p[0]*v_tmp[7])/k[0] - (1.0*p[1]*v_tmp[8])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*v_tmp[11])/(pow(k[0],2));
Jv_tmp[12] = p[0]*v_tmp[1] + p[0]*v_tmp[2] + p[0]*v_tmp[4] + p[1]*v_tmp[5] - (1.0*p[1]*v_tmp[8])/k[0] - (1.0*p[0]*v_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*v_tmp[12])/(pow(k[0],2));
Jv_tmp[13] = p[0]*v_tmp[6] + p[1]*v_tmp[8] - (1.0*(k[0]*p[0] + k[0]*p[1])*v_tmp[13])/k[0];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_enhancer_667(N_Vector vB, N_Vector JvB, realtype t,
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
JvB_tmp[0] = 4.0*p[0]*vB_tmp[0] - 1.0*p[0]*vB_tmp[9] - 1.0*p[0]*vB_tmp[10];
JvB_tmp[1] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[1])/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[11] - 1.0*p[1]*vB_tmp[10] - 1.0*p[0]*vB_tmp[12] - 1.0*p[1]*vB_tmp[9];
JvB_tmp[2] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[2])/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[11] - 1.0*p[1]*vB_tmp[10] - 1.0*p[0]*vB_tmp[12] - 1.0*p[1]*vB_tmp[9];
JvB_tmp[3] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[3])/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[11] - 1.0*p[1]*vB_tmp[9];
JvB_tmp[4] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[4])/(pow(k[0],2)) - 1.0*p[0]*vB_tmp[12] - 1.0*p[1]*vB_tmp[10];
JvB_tmp[5] = 4.0*p[1]*vB_tmp[5] - 1.0*p[1]*vB_tmp[11] - 1.0*p[1]*vB_tmp[12];
JvB_tmp[6] = 2.0*p[0]*vB_tmp[6] - 1.0*p[0]*vB_tmp[7] - 1.0*p[0]*vB_tmp[13] - (2.0*p[0]*vB_tmp[0])/k[0] - (1.0*p[0]*vB_tmp[3])/k[0] - (1.0*p[0]*vB_tmp[4])/k[0] + (p[0]*vB_tmp[9])/k[0] + (p[0]*vB_tmp[10])/k[0];
JvB_tmp[7] = (p[1]*vB_tmp[9])/k[0] - 1.0*p[0]*vB_tmp[8] - (1.0*p[1]*vB_tmp[0])/k[0] - (1.0*p[0]*vB_tmp[5])/k[0] - 1.0*p[1]*vB_tmp[6] + (p[0]*vB_tmp[11])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*vB_tmp[3])/(pow(k[0],2)) + ((k[0]*p[0] + k[0]*p[1])*vB_tmp[7])/k[0];
JvB_tmp[8] = 2.0*p[1]*vB_tmp[8] - 1.0*p[1]*vB_tmp[7] - 1.0*p[1]*vB_tmp[13] - (1.0*p[1]*vB_tmp[3])/k[0] - (1.0*p[1]*vB_tmp[4])/k[0] - (2.0*p[1]*vB_tmp[5])/k[0] + (p[1]*vB_tmp[11])/k[0] + (p[1]*vB_tmp[12])/k[0];
JvB_tmp[9] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*vB_tmp[9])/(pow(k[0],2)) - 2.0*p[1]*vB_tmp[0] - 1.0*p[0]*vB_tmp[2] - 2.0*p[0]*vB_tmp[3] - 1.0*p[0]*vB_tmp[1];
JvB_tmp[10] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*vB_tmp[10])/(pow(k[0],2)) - 2.0*p[1]*vB_tmp[0] - 1.0*p[0]*vB_tmp[2] - 2.0*p[0]*vB_tmp[4] - 1.0*p[0]*vB_tmp[1];
JvB_tmp[11] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*vB_tmp[11])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[2] - 2.0*p[1]*vB_tmp[3] - 2.0*p[0]*vB_tmp[5] - 1.0*p[1]*vB_tmp[1];
JvB_tmp[12] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*vB_tmp[12])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[2] - 2.0*p[0]*vB_tmp[5] - 2.0*p[1]*vB_tmp[4] - 1.0*p[1]*vB_tmp[1];
JvB_tmp[13] = (p[1]*vB_tmp[10])/k[0] - 1.0*p[0]*vB_tmp[8] - (1.0*p[1]*vB_tmp[0])/k[0] - (1.0*p[0]*vB_tmp[5])/k[0] - 1.0*p[1]*vB_tmp[6] + (p[0]*vB_tmp[12])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*vB_tmp[4])/(pow(k[0],2)) + ((k[0]*p[0] + k[0]*p[1])*vB_tmp[13])/k[0];

  for (ix=0; ix<14; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_enhancer_667(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_enhancer_667(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_enhancer_667(long int N, realtype t, N_Vector x,
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
J->data[0] = -4.0*p[0];
J->data[9] = p[0];
J->data[10] = p[0];
J->data[15] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[23] = p[1];
J->data[24] = p[1];
J->data[25] = p[0];
J->data[26] = p[0];
J->data[30] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[37] = p[1];
J->data[38] = p[1];
J->data[39] = p[0];
J->data[40] = p[0];
J->data[45] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[51] = p[1];
J->data[53] = p[0];
J->data[60] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[66] = p[1];
J->data[68] = p[0];
J->data[75] = -4.0*p[1];
J->data[81] = p[1];
J->data[82] = p[1];
J->data[84] = (2.0*p[0])/k[0];
J->data[87] = p[0]/k[0];
J->data[88] = p[0]/k[0];
J->data[90] = -2.0*p[0];
J->data[91] = p[0];
J->data[93] = -(1.0*p[0])/k[0];
J->data[94] = -(1.0*p[0])/k[0];
J->data[97] = p[0];
J->data[98] = p[1]/k[0];
J->data[101] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[103] = p[0]/k[0];
J->data[104] = p[1];
J->data[105] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];
J->data[106] = p[0];
J->data[107] = -(1.0*p[1])/k[0];
J->data[109] = -(1.0*p[0])/k[0];
J->data[115] = p[1]/k[0];
J->data[116] = p[1]/k[0];
J->data[117] = (2.0*p[1])/k[0];
J->data[119] = p[1];
J->data[120] = -2.0*p[1];
J->data[123] = -(1.0*p[1])/k[0];
J->data[124] = -(1.0*p[1])/k[0];
J->data[125] = p[1];
J->data[126] = 2.0*p[1];
J->data[127] = p[0];
J->data[128] = p[0];
J->data[129] = 2.0*p[0];
J->data[135] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[140] = 2.0*p[1];
J->data[141] = p[0];
J->data[142] = p[0];
J->data[144] = 2.0*p[0];
J->data[150] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[155] = p[1];
J->data[156] = p[1];
J->data[157] = 2.0*p[1];
J->data[159] = 2.0*p[0];
J->data[165] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[169] = p[1];
J->data[170] = p[1];
J->data[172] = 2.0*p[1];
J->data[173] = 2.0*p[0];
J->data[180] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[182] = p[1]/k[0];
J->data[186] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[187] = p[0]/k[0];
J->data[188] = p[1];
J->data[190] = p[0];
J->data[192] = -(1.0*p[1])/k[0];
J->data[194] = -(1.0*p[0])/k[0];
J->data[195] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];

  for (iJ=0; iJ<196; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_enhancer_667(realtype t, N_Vector x,
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
  J->rowvals[1] = 9;
  J->rowvals[2] = 10;
  J->rowvals[3] = 1;
  J->rowvals[4] = 9;
  J->rowvals[5] = 10;
  J->rowvals[6] = 11;
  J->rowvals[7] = 12;
  J->rowvals[8] = 2;
  J->rowvals[9] = 9;
  J->rowvals[10] = 10;
  J->rowvals[11] = 11;
  J->rowvals[12] = 12;
  J->rowvals[13] = 3;
  J->rowvals[14] = 9;
  J->rowvals[15] = 11;
  J->rowvals[16] = 4;
  J->rowvals[17] = 10;
  J->rowvals[18] = 12;
  J->rowvals[19] = 5;
  J->rowvals[20] = 11;
  J->rowvals[21] = 12;
  J->rowvals[22] = 0;
  J->rowvals[23] = 3;
  J->rowvals[24] = 4;
  J->rowvals[25] = 6;
  J->rowvals[26] = 7;
  J->rowvals[27] = 9;
  J->rowvals[28] = 10;
  J->rowvals[29] = 13;
  J->rowvals[30] = 0;
  J->rowvals[31] = 3;
  J->rowvals[32] = 5;
  J->rowvals[33] = 6;
  J->rowvals[34] = 7;
  J->rowvals[35] = 8;
  J->rowvals[36] = 9;
  J->rowvals[37] = 11;
  J->rowvals[38] = 3;
  J->rowvals[39] = 4;
  J->rowvals[40] = 5;
  J->rowvals[41] = 7;
  J->rowvals[42] = 8;
  J->rowvals[43] = 11;
  J->rowvals[44] = 12;
  J->rowvals[45] = 13;
  J->rowvals[46] = 0;
  J->rowvals[47] = 1;
  J->rowvals[48] = 2;
  J->rowvals[49] = 3;
  J->rowvals[50] = 9;
  J->rowvals[51] = 0;
  J->rowvals[52] = 1;
  J->rowvals[53] = 2;
  J->rowvals[54] = 4;
  J->rowvals[55] = 10;
  J->rowvals[56] = 1;
  J->rowvals[57] = 2;
  J->rowvals[58] = 3;
  J->rowvals[59] = 5;
  J->rowvals[60] = 11;
  J->rowvals[61] = 1;
  J->rowvals[62] = 2;
  J->rowvals[63] = 4;
  J->rowvals[64] = 5;
  J->rowvals[65] = 12;
  J->rowvals[66] = 0;
  J->rowvals[67] = 4;
  J->rowvals[68] = 5;
  J->rowvals[69] = 6;
  J->rowvals[70] = 8;
  J->rowvals[71] = 10;
  J->rowvals[72] = 12;
  J->rowvals[73] = 13;
  J->colptrs[0] = 0;
  J->colptrs[1] = 3;
  J->colptrs[2] = 8;
  J->colptrs[3] = 13;
  J->colptrs[4] = 16;
  J->colptrs[5] = 19;
  J->colptrs[6] = 22;
  J->colptrs[7] = 30;
  J->colptrs[8] = 38;
  J->colptrs[9] = 46;
  J->colptrs[10] = 51;
  J->colptrs[11] = 56;
  J->colptrs[12] = 61;
  J->colptrs[13] = 66;
  J->colptrs[14] = 74;
J->data[0] = -4.0*p[0];
J->data[1] = p[0];
J->data[2] = p[0];
J->data[3] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[4] = p[1];
J->data[5] = p[1];
J->data[6] = p[0];
J->data[7] = p[0];
J->data[8] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[9] = p[1];
J->data[10] = p[1];
J->data[11] = p[0];
J->data[12] = p[0];
J->data[13] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[14] = p[1];
J->data[15] = p[0];
J->data[16] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[17] = p[1];
J->data[18] = p[0];
J->data[19] = -4.0*p[1];
J->data[20] = p[1];
J->data[21] = p[1];
J->data[22] = (2.0*p[0])/k[0];
J->data[23] = p[0]/k[0];
J->data[24] = p[0]/k[0];
J->data[25] = -2.0*p[0];
J->data[26] = p[0];
J->data[27] = -(1.0*p[0])/k[0];
J->data[28] = -(1.0*p[0])/k[0];
J->data[29] = p[0];
J->data[30] = p[1]/k[0];
J->data[31] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[32] = p[0]/k[0];
J->data[33] = p[1];
J->data[34] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];
J->data[35] = p[0];
J->data[36] = -(1.0*p[1])/k[0];
J->data[37] = -(1.0*p[0])/k[0];
J->data[38] = p[1]/k[0];
J->data[39] = p[1]/k[0];
J->data[40] = (2.0*p[1])/k[0];
J->data[41] = p[1];
J->data[42] = -2.0*p[1];
J->data[43] = -(1.0*p[1])/k[0];
J->data[44] = -(1.0*p[1])/k[0];
J->data[45] = p[1];
J->data[46] = 2.0*p[1];
J->data[47] = p[0];
J->data[48] = p[0];
J->data[49] = 2.0*p[0];
J->data[50] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[51] = 2.0*p[1];
J->data[52] = p[0];
J->data[53] = p[0];
J->data[54] = 2.0*p[0];
J->data[55] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[56] = p[1];
J->data[57] = p[1];
J->data[58] = 2.0*p[1];
J->data[59] = 2.0*p[0];
J->data[60] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[61] = p[1];
J->data[62] = p[1];
J->data[63] = 2.0*p[1];
J->data[64] = 2.0*p[0];
J->data[65] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[66] = p[1]/k[0];
J->data[67] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[68] = p[0]/k[0];
J->data[69] = p[1];
J->data[70] = p[0];
J->data[71] = -(1.0*p[1])/k[0];
J->data[72] = -(1.0*p[0])/k[0];
J->data[73] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];
  return(0);
}


 int JBBand_enhancer_667(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_enhancer_667(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_enhancer_667(long int N, realtype t, N_Vector x,
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
JB->data[6] = -(2.0*p[0])/k[0];
JB->data[7] = -(1.0*p[1])/k[0];
JB->data[9] = -2.0*p[1];
JB->data[10] = -2.0*p[1];
JB->data[13] = -(1.0*p[1])/k[0];
JB->data[15] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[23] = -1.0*p[0];
JB->data[24] = -1.0*p[0];
JB->data[25] = -1.0*p[1];
JB->data[26] = -1.0*p[1];
JB->data[30] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[37] = -1.0*p[0];
JB->data[38] = -1.0*p[0];
JB->data[39] = -1.0*p[1];
JB->data[40] = -1.0*p[1];
JB->data[45] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[48] = -(1.0*p[0])/k[0];
JB->data[49] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/(pow(k[0],2));
JB->data[50] = -(1.0*p[1])/k[0];
JB->data[51] = -2.0*p[0];
JB->data[53] = -2.0*p[1];
JB->data[60] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[62] = -(1.0*p[0])/k[0];
JB->data[64] = -(1.0*p[1])/k[0];
JB->data[66] = -2.0*p[0];
JB->data[68] = -2.0*p[1];
JB->data[69] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/(pow(k[0],2));
JB->data[75] = 4.0*p[1];
JB->data[77] = -(1.0*p[0])/k[0];
JB->data[78] = -(2.0*p[1])/k[0];
JB->data[81] = -2.0*p[0];
JB->data[82] = -2.0*p[0];
JB->data[83] = -(1.0*p[0])/k[0];
JB->data[90] = 2.0*p[0];
JB->data[91] = -1.0*p[1];
JB->data[97] = -1.0*p[1];
JB->data[104] = -1.0*p[0];
JB->data[105] = (k[0]*p[0] + k[0]*p[1])/k[0];
JB->data[106] = -1.0*p[1];
JB->data[119] = -1.0*p[0];
JB->data[120] = 2.0*p[1];
JB->data[125] = -1.0*p[0];
JB->data[126] = -1.0*p[0];
JB->data[127] = -1.0*p[1];
JB->data[128] = -1.0*p[1];
JB->data[129] = -1.0*p[1];
JB->data[132] = p[0]/k[0];
JB->data[133] = p[1]/k[0];
JB->data[135] = (3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[140] = -1.0*p[0];
JB->data[141] = -1.0*p[1];
JB->data[142] = -1.0*p[1];
JB->data[144] = -1.0*p[1];
JB->data[146] = p[0]/k[0];
JB->data[150] = (3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[153] = p[1]/k[0];
JB->data[155] = -1.0*p[0];
JB->data[156] = -1.0*p[0];
JB->data[157] = -1.0*p[0];
JB->data[159] = -1.0*p[1];
JB->data[161] = p[0]/k[0];
JB->data[162] = p[1]/k[0];
JB->data[165] = ((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[169] = -1.0*p[0];
JB->data[170] = -1.0*p[0];
JB->data[172] = -1.0*p[0];
JB->data[173] = -1.0*p[1];
JB->data[176] = p[1]/k[0];
JB->data[180] = ((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[181] = p[0]/k[0];
JB->data[188] = -1.0*p[0];
JB->data[190] = -1.0*p[1];
JB->data[195] = (k[0]*p[0] + k[0]*p[1])/k[0];

  for (iJ=0; iJ<196; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_enhancer_667(realtype t, N_Vector x,
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
  JB->rowvals[1] = 6;
  JB->rowvals[2] = 7;
  JB->rowvals[3] = 9;
  JB->rowvals[4] = 10;
  JB->rowvals[5] = 13;
  JB->rowvals[6] = 1;
  JB->rowvals[7] = 9;
  JB->rowvals[8] = 10;
  JB->rowvals[9] = 11;
  JB->rowvals[10] = 12;
  JB->rowvals[11] = 2;
  JB->rowvals[12] = 9;
  JB->rowvals[13] = 10;
  JB->rowvals[14] = 11;
  JB->rowvals[15] = 12;
  JB->rowvals[16] = 3;
  JB->rowvals[17] = 6;
  JB->rowvals[18] = 7;
  JB->rowvals[19] = 8;
  JB->rowvals[20] = 9;
  JB->rowvals[21] = 11;
  JB->rowvals[22] = 4;
  JB->rowvals[23] = 6;
  JB->rowvals[24] = 8;
  JB->rowvals[25] = 10;
  JB->rowvals[26] = 12;
  JB->rowvals[27] = 13;
  JB->rowvals[28] = 5;
  JB->rowvals[29] = 7;
  JB->rowvals[30] = 8;
  JB->rowvals[31] = 11;
  JB->rowvals[32] = 12;
  JB->rowvals[33] = 13;
  JB->rowvals[34] = 6;
  JB->rowvals[35] = 7;
  JB->rowvals[36] = 13;
  JB->rowvals[37] = 6;
  JB->rowvals[38] = 7;
  JB->rowvals[39] = 8;
  JB->rowvals[40] = 7;
  JB->rowvals[41] = 8;
  JB->rowvals[42] = 13;
  JB->rowvals[43] = 0;
  JB->rowvals[44] = 1;
  JB->rowvals[45] = 2;
  JB->rowvals[46] = 3;
  JB->rowvals[47] = 6;
  JB->rowvals[48] = 7;
  JB->rowvals[49] = 9;
  JB->rowvals[50] = 0;
  JB->rowvals[51] = 1;
  JB->rowvals[52] = 2;
  JB->rowvals[53] = 4;
  JB->rowvals[54] = 6;
  JB->rowvals[55] = 10;
  JB->rowvals[56] = 13;
  JB->rowvals[57] = 1;
  JB->rowvals[58] = 2;
  JB->rowvals[59] = 3;
  JB->rowvals[60] = 5;
  JB->rowvals[61] = 7;
  JB->rowvals[62] = 8;
  JB->rowvals[63] = 11;
  JB->rowvals[64] = 1;
  JB->rowvals[65] = 2;
  JB->rowvals[66] = 4;
  JB->rowvals[67] = 5;
  JB->rowvals[68] = 8;
  JB->rowvals[69] = 12;
  JB->rowvals[70] = 13;
  JB->rowvals[71] = 6;
  JB->rowvals[72] = 8;
  JB->rowvals[73] = 13;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 6;
  JB->colptrs[2] = 11;
  JB->colptrs[3] = 16;
  JB->colptrs[4] = 22;
  JB->colptrs[5] = 28;
  JB->colptrs[6] = 34;
  JB->colptrs[7] = 37;
  JB->colptrs[8] = 40;
  JB->colptrs[9] = 43;
  JB->colptrs[10] = 50;
  JB->colptrs[11] = 57;
  JB->colptrs[12] = 64;
  JB->colptrs[13] = 71;
  JB->colptrs[14] = 74;
  return(0);
}


 int sx_enhancer_667(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[9] - 4.0*p[0]*sx_tmp[0] + 2.0*p[1]*sx_tmp[10] + (2.0*k[0]*x_tmp[6] - 4.0*(pow(k[0],2))*x_tmp[0])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[6])/k[0] + (p[1]*sx_tmp[7])/k[0] + (p[1]*sx_tmp[13])/k[0];
sxdot_tmp[1] = p[0]*sx_tmp[9] + p[0]*sx_tmp[10] + p[1]*sx_tmp[11] + p[1]*sx_tmp[12] + ((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[10])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[1])/(pow(k[0],2));
sxdot_tmp[2] = p[0]*sx_tmp[9] + p[0]*sx_tmp[10] + p[1]*sx_tmp[11] + p[1]*sx_tmp[12] + ((pow(k[0],2))*x_tmp[9] - 2.0*(pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[10])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[2])/(pow(k[0],2));
sxdot_tmp[3] = 2.0*p[0]*sx_tmp[9] + 2.0*p[1]*sx_tmp[11] + (k[0]*x_tmp[6] + k[0]*x_tmp[7] - 2.0*(pow(k[0],2))*x_tmp[3] + 2.0*(pow(k[0],2))*x_tmp[9])/(pow(k[0],2)) + (p[0]*sx_tmp[6])/k[0] + (p[1]*sx_tmp[8])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = 2.0*p[0]*sx_tmp[10] + 2.0*p[1]*sx_tmp[12] + (k[0]*x_tmp[6] + k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[4] + 2.0*(pow(k[0],2))*x_tmp[10])/(pow(k[0],2)) + (p[0]*sx_tmp[6])/k[0] + (p[1]*sx_tmp[8])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = 2.0*p[0]*sx_tmp[11] - 4.0*p[1]*sx_tmp[5] + 2.0*p[0]*sx_tmp[12] + (k[0]*x_tmp[7] + k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[11] + 2.0*(pow(k[0],2))*x_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[7])/k[0] + (2.0*p[1]*sx_tmp[8])/k[0] + (p[0]*sx_tmp[13])/k[0];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 2.0*p[0]*sx_tmp[6] + p[1]*sx_tmp[13] - 2.0*x_tmp[6];
sxdot_tmp[7] = p[0]*sx_tmp[6] + p[1]*sx_tmp[8] + (k[0]*x_tmp[6] - 1.0*k[0]*x_tmp[7])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[7])/k[0];
sxdot_tmp[8] = p[0]*sx_tmp[7] - 2.0*p[1]*sx_tmp[8] + p[0]*sx_tmp[13] + (k[0]*x_tmp[7] + k[0]*x_tmp[13])/k[0];
sxdot_tmp[9] = p[0]*sx_tmp[0] + p[1]*sx_tmp[1] + p[1]*sx_tmp[2] + p[1]*sx_tmp[3] - (1.0*(k[0]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[9]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[6])/k[0] - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[9])/(pow(k[0],2));
sxdot_tmp[10] = p[0]*sx_tmp[0] + p[1]*sx_tmp[1] + p[1]*sx_tmp[2] + p[1]*sx_tmp[4] - (1.0*(k[0]*x_tmp[6] - 1.0*(pow(k[0],2))*x_tmp[0] + 3.0*(pow(k[0],2))*x_tmp[10]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[6])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[10])/(pow(k[0],2));
sxdot_tmp[11] = p[0]*sx_tmp[1] + p[0]*sx_tmp[2] + p[0]*sx_tmp[3] + p[1]*sx_tmp[5] + ((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[11])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7])/k[0] - (1.0*p[1]*sx_tmp[8])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[11])/(pow(k[0],2));
sxdot_tmp[12] = p[0]*sx_tmp[1] + p[0]*sx_tmp[2] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] + ((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[12])/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[8])/k[0] - (1.0*p[0]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[12])/(pow(k[0],2));
sxdot_tmp[13] = p[0]*sx_tmp[6] + p[1]*sx_tmp[8] + (k[0]*x_tmp[6] - 1.0*k[0]*x_tmp[13])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[13])/k[0];

  } break;

  case 1: {
sxdot_tmp[0] = 2.0*p[1]*sx_tmp[9] - 4.0*p[0]*sx_tmp[0] + 2.0*p[1]*sx_tmp[10] + (k[0]*x_tmp[7] + k[0]*x_tmp[13] + 2.0*(pow(k[0],2))*x_tmp[9] + 2.0*(pow(k[0],2))*x_tmp[10])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[6])/k[0] + (p[1]*sx_tmp[7])/k[0] + (p[1]*sx_tmp[13])/k[0];
sxdot_tmp[1] = p[0]*sx_tmp[9] + p[0]*sx_tmp[10] + p[1]*sx_tmp[11] + p[1]*sx_tmp[12] + ((pow(k[0],2))*x_tmp[11] - 2.0*(pow(k[0],2))*x_tmp[1] + (pow(k[0],2))*x_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[1])/(pow(k[0],2));
sxdot_tmp[2] = p[0]*sx_tmp[9] + p[0]*sx_tmp[10] + p[1]*sx_tmp[11] + p[1]*sx_tmp[12] + ((pow(k[0],2))*x_tmp[11] - 2.0*(pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[12])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[2])/(pow(k[0],2));
sxdot_tmp[3] = 2.0*p[0]*sx_tmp[9] + 2.0*p[1]*sx_tmp[11] + (k[0]*x_tmp[7] + k[0]*x_tmp[8] - 2.0*(pow(k[0],2))*x_tmp[3] + 2.0*(pow(k[0],2))*x_tmp[11])/(pow(k[0],2)) + (p[0]*sx_tmp[6])/k[0] + (p[1]*sx_tmp[8])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[7])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = 2.0*p[0]*sx_tmp[10] + 2.0*p[1]*sx_tmp[12] + (k[0]*x_tmp[8] + k[0]*x_tmp[13] - 2.0*(pow(k[0],2))*x_tmp[4] + 2.0*(pow(k[0],2))*x_tmp[12])/(pow(k[0],2)) + (p[0]*sx_tmp[6])/k[0] + (p[1]*sx_tmp[8])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[13])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[4])/(pow(k[0],2));
sxdot_tmp[5] = 2.0*p[0]*sx_tmp[11] - 4.0*p[1]*sx_tmp[5] + 2.0*p[0]*sx_tmp[12] + (2.0*k[0]*x_tmp[8] - 4.0*(pow(k[0],2))*x_tmp[5])/(pow(k[0],2)) + (p[0]*sx_tmp[7])/k[0] + (2.0*p[1]*sx_tmp[8])/k[0] + (p[0]*sx_tmp[13])/k[0];
sxdot_tmp[6] = p[1]*sx_tmp[7] - 2.0*p[0]*sx_tmp[6] + p[1]*sx_tmp[13] + (k[0]*x_tmp[7] + k[0]*x_tmp[13])/k[0];
sxdot_tmp[7] = p[0]*sx_tmp[6] + p[1]*sx_tmp[8] - (1.0*(k[0]*x_tmp[7] - 1.0*k[0]*x_tmp[8]))/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[7])/k[0];
sxdot_tmp[8] = p[0]*sx_tmp[7] - 2.0*p[1]*sx_tmp[8] + p[0]*sx_tmp[13] - 2.0*x_tmp[8];
sxdot_tmp[9] = p[0]*sx_tmp[0] + p[1]*sx_tmp[1] + p[1]*sx_tmp[2] + p[1]*sx_tmp[3] + ((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[7] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[9])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[6])/k[0] - (1.0*p[1]*sx_tmp[7])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[9])/(pow(k[0],2));
sxdot_tmp[10] = p[0]*sx_tmp[0] + p[1]*sx_tmp[1] + p[1]*sx_tmp[2] + p[1]*sx_tmp[4] + ((pow(k[0],2))*x_tmp[1] - 1.0*k[0]*x_tmp[13] + (pow(k[0],2))*x_tmp[2] + (pow(k[0],2))*x_tmp[4] - 1.0*(pow(k[0],2))*x_tmp[10])/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[6])/k[0] - (1.0*p[1]*sx_tmp[13])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[10])/(pow(k[0],2));
sxdot_tmp[11] = p[0]*sx_tmp[1] + p[0]*sx_tmp[2] + p[0]*sx_tmp[3] + p[1]*sx_tmp[5] - (1.0*(k[0]*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[11]))/(pow(k[0],2)) - (1.0*p[0]*sx_tmp[7])/k[0] - (1.0*p[1]*sx_tmp[8])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[11])/(pow(k[0],2));
sxdot_tmp[12] = p[0]*sx_tmp[1] + p[0]*sx_tmp[2] + p[0]*sx_tmp[4] + p[1]*sx_tmp[5] - (1.0*(k[0]*x_tmp[8] - 1.0*(pow(k[0],2))*x_tmp[5] + 3.0*(pow(k[0],2))*x_tmp[12]))/(pow(k[0],2)) - (1.0*p[1]*sx_tmp[8])/k[0] - (1.0*p[0]*sx_tmp[13])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[12])/(pow(k[0],2));
sxdot_tmp[13] = p[0]*sx_tmp[6] + p[1]*sx_tmp[8] + (k[0]*x_tmp[8] - 1.0*k[0]*x_tmp[13])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[13])/k[0];

  } break;

  }
 for (ix=0; ix<14; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_enhancer_667(int ip, N_Vector sx0, void *user_data)
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


void y_enhancer_667(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*6];
y[it+nt*1] = x[it+nt*0];
    
    return;
}


void dydp_enhancer_667(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }
  
  return;
}


void dydx_enhancer_667(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*28);
dydx[1] = 1.0;
dydx[12] = 1.0;
  
  return;
}


void sy_enhancer_667(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(0+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(6+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(0+np*nx)];

  } break;

  }
  
  return;
}
int root_enhancer_667(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_enhancer_667(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root_enhancer_667(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval_enhancer_667(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval_enhancer_667(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
void deltadisc_enhancer_667(double t, int idisc, N_Vector x, void *user_data){
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
void sdeltadisc_enhancer_667(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data){
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


void dxdotdp_enhancer_667(double t, int nt, int it, double *dxdotdp, double *p, double *k, double *u, double *x, int *plist, int np, int nx){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
dxdotdp[(0+ip*nx)] = (2.0*k[0]*x[it+nt*6] - 4.0*(pow(k[0],2))*x[it+nt*0])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = ((pow(k[0],2))*x[it+nt*9] - 2.0*(pow(k[0],2))*x[it+nt*1] + (pow(k[0],2))*x[it+nt*10])/(pow(k[0],2));
dxdotdp[(2+ip*nx)] = ((pow(k[0],2))*x[it+nt*9] - 2.0*(pow(k[0],2))*x[it+nt*2] + (pow(k[0],2))*x[it+nt*10])/(pow(k[0],2));
dxdotdp[(3+ip*nx)] = (k[0]*x[it+nt*6] + k[0]*x[it+nt*7] - 2.0*(pow(k[0],2))*x[it+nt*3] + 2.0*(pow(k[0],2))*x[it+nt*9])/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = (k[0]*x[it+nt*6] + k[0]*x[it+nt*13] - 2.0*(pow(k[0],2))*x[it+nt*4] + 2.0*(pow(k[0],2))*x[it+nt*10])/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = (k[0]*x[it+nt*7] + k[0]*x[it+nt*13] + 2.0*(pow(k[0],2))*x[it+nt*11] + 2.0*(pow(k[0],2))*x[it+nt*12])/(pow(k[0],2));
dxdotdp[(6+ip*nx)] = -2.0*x[it+nt*6];
dxdotdp[(7+ip*nx)] = (k[0]*x[it+nt*6] - 1.0*k[0]*x[it+nt*7])/k[0];
dxdotdp[(8+ip*nx)] = (k[0]*x[it+nt*7] + k[0]*x[it+nt*13])/k[0];
dxdotdp[(9+ip*nx)] = -(1.0*(k[0]*x[it+nt*6] - 1.0*(pow(k[0],2))*x[it+nt*0] + 3.0*(pow(k[0],2))*x[it+nt*9]))/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = -(1.0*(k[0]*x[it+nt*6] - 1.0*(pow(k[0],2))*x[it+nt*0] + 3.0*(pow(k[0],2))*x[it+nt*10]))/(pow(k[0],2));
dxdotdp[(11+ip*nx)] = ((pow(k[0],2))*x[it+nt*1] - 1.0*k[0]*x[it+nt*7] + (pow(k[0],2))*x[it+nt*2] + (pow(k[0],2))*x[it+nt*3] - 1.0*(pow(k[0],2))*x[it+nt*11])/(pow(k[0],2));
dxdotdp[(12+ip*nx)] = ((pow(k[0],2))*x[it+nt*1] - 1.0*k[0]*x[it+nt*13] + (pow(k[0],2))*x[it+nt*2] + (pow(k[0],2))*x[it+nt*4] - 1.0*(pow(k[0],2))*x[it+nt*12])/(pow(k[0],2));
dxdotdp[(13+ip*nx)] = (k[0]*x[it+nt*6] - 1.0*k[0]*x[it+nt*13])/k[0];

  } break;

  case 1: {
dxdotdp[(0+ip*nx)] = (k[0]*x[it+nt*7] + k[0]*x[it+nt*13] + 2.0*(pow(k[0],2))*x[it+nt*9] + 2.0*(pow(k[0],2))*x[it+nt*10])/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = ((pow(k[0],2))*x[it+nt*11] - 2.0*(pow(k[0],2))*x[it+nt*1] + (pow(k[0],2))*x[it+nt*12])/(pow(k[0],2));
dxdotdp[(2+ip*nx)] = ((pow(k[0],2))*x[it+nt*11] - 2.0*(pow(k[0],2))*x[it+nt*2] + (pow(k[0],2))*x[it+nt*12])/(pow(k[0],2));
dxdotdp[(3+ip*nx)] = (k[0]*x[it+nt*7] + k[0]*x[it+nt*8] - 2.0*(pow(k[0],2))*x[it+nt*3] + 2.0*(pow(k[0],2))*x[it+nt*11])/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = (k[0]*x[it+nt*8] + k[0]*x[it+nt*13] - 2.0*(pow(k[0],2))*x[it+nt*4] + 2.0*(pow(k[0],2))*x[it+nt*12])/(pow(k[0],2));
dxdotdp[(5+ip*nx)] = (2.0*k[0]*x[it+nt*8] - 4.0*(pow(k[0],2))*x[it+nt*5])/(pow(k[0],2));
dxdotdp[(6+ip*nx)] = (k[0]*x[it+nt*7] + k[0]*x[it+nt*13])/k[0];
dxdotdp[(7+ip*nx)] = -(1.0*(k[0]*x[it+nt*7] - 1.0*k[0]*x[it+nt*8]))/k[0];
dxdotdp[(8+ip*nx)] = -2.0*x[it+nt*8];
dxdotdp[(9+ip*nx)] = ((pow(k[0],2))*x[it+nt*1] - 1.0*k[0]*x[it+nt*7] + (pow(k[0],2))*x[it+nt*2] + (pow(k[0],2))*x[it+nt*3] - 1.0*(pow(k[0],2))*x[it+nt*9])/(pow(k[0],2));
dxdotdp[(10+ip*nx)] = ((pow(k[0],2))*x[it+nt*1] - 1.0*k[0]*x[it+nt*13] + (pow(k[0],2))*x[it+nt*2] + (pow(k[0],2))*x[it+nt*4] - 1.0*(pow(k[0],2))*x[it+nt*10])/(pow(k[0],2));
dxdotdp[(11+ip*nx)] = -(1.0*(k[0]*x[it+nt*8] - 1.0*(pow(k[0],2))*x[it+nt*5] + 3.0*(pow(k[0],2))*x[it+nt*11]))/(pow(k[0],2));
dxdotdp[(12+ip*nx)] = -(1.0*(k[0]*x[it+nt*8] - 1.0*(pow(k[0],2))*x[it+nt*5] + 3.0*(pow(k[0],2))*x[it+nt*12]))/(pow(k[0],2));
dxdotdp[(13+ip*nx)] = (k[0]*x[it+nt*8] - 1.0*k[0]*x[it+nt*13])/k[0];

  } break;

  }
  }
  
  return;
}
