#include "enhancer_666.h"
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


 int xdot_enhancer_666(realtype t, N_Vector x, N_Vector xdot, void *user_data)
{
  int ix;
  UserData data = (UserData) user_data;
  double *qpositivex = data->qpositivex;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xdot_tmp = N_VGetArrayPointer(xdot);
  memset(xdot_tmp,0,sizeof(double)*9);
xdot_tmp[0] = -(1.0*(2.0*(pow(k[0],2))*p[0]*x_tmp[0] + 2.0*(pow(k[0],2))*p[1]*x_tmp[0] - 1.0*(pow(k[0],2))*p[1]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[8]))/(pow(k[0],2));
xdot_tmp[1] = (2.0*(pow(k[0],2))*p[0]*x_tmp[3] - 4.0*(pow(k[0],2))*p[1]*x_tmp[1] + 2.0*k[0]*p[1]*x_tmp[2] + k[0]*p[0]*x_tmp[5])/(pow(k[0],2));
xdot_tmp[2] = -(1.0*(2.0*k[0]*p[1]*x_tmp[2] - 1.0*k[0]*p[0]*x_tmp[5]))/k[0];
xdot_tmp[3] = -(1.0*((pow(k[0],2))*p[0]*x_tmp[3] - 2.0*(pow(k[0],2))*p[1]*x_tmp[1] - 2.0*(pow(k[0],2))*p[0]*x_tmp[0] + 3.0*(pow(k[0],2))*p[1]*x_tmp[3] - 1.0*(pow(k[0],2))*p[0]*x_tmp[7] + 2.0*k[0]*p[1]*x_tmp[2] + k[0]*p[0]*x_tmp[5]))/(pow(k[0],2));
xdot_tmp[4] = -(1.0*(2.0*k[0]*p[0]*x_tmp[4] - 1.0*k[0]*p[1]*x_tmp[5]))/k[0];
xdot_tmp[5] = (2.0*k[0]*p[1]*x_tmp[2] + 2.0*k[0]*p[0]*x_tmp[4] - 1.0*k[0]*p[0]*x_tmp[5] - 1.0*k[0]*p[1]*x_tmp[5])/k[0];
xdot_tmp[6] = (2.0*(pow(k[0],2))*p[1]*x_tmp[8] - 4.0*(pow(k[0],2))*p[0]*x_tmp[6] + 2.0*k[0]*p[0]*x_tmp[4] + k[0]*p[1]*x_tmp[5])/(pow(k[0],2));
xdot_tmp[7] = (4.0*(pow(k[0],2))*p[1]*x_tmp[3] - 2.0*(pow(k[0],2))*p[0]*x_tmp[7] + 4.0*(pow(k[0],2))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[1]*x_tmp[7] + 2.0*k[0]*p[1]*x_tmp[2] + 2.0*k[0]*p[0]*x_tmp[4] + k[0]*p[0]*x_tmp[5] + k[0]*p[1]*x_tmp[5])/(pow(k[0],2));
xdot_tmp[8] = -(1.0*(3.0*(pow(k[0],2))*p[0]*x_tmp[8] - 2.0*(pow(k[0],2))*p[0]*x_tmp[6] - 2.0*(pow(k[0],2))*p[1]*x_tmp[0] - 1.0*(pow(k[0],2))*p[1]*x_tmp[7] + (pow(k[0],2))*p[1]*x_tmp[8] + 2.0*k[0]*p[0]*x_tmp[4] + k[0]*p[1]*x_tmp[5]))/(pow(k[0],2));

  for (ix=0; ix<9; ix++) {
    if(mxIsNaN(xdot_tmp[ix])) xdot_tmp[ix] = 0.0;
    if(qpositivex[ix]>0.5 && x_tmp[ix]<0.0 && xdot_tmp[ix]<0.0) xdot_tmp[ix] = -xdot_tmp[ix];
  }

  return(0);
}


 int xBdot_enhancer_666(realtype t, N_Vector x, N_Vector xB, N_Vector xBdot, void *user_data)
{
  int ixB;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double *xB_tmp = N_VGetArrayPointer(xB);
  double *xBdot_tmp = N_VGetArrayPointer(xBdot);
  memset(xBdot_tmp,0,sizeof(double)*18);
xBdot_tmp[0] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[0])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[8] - 2.0*p[0]*xB_tmp[3];
xBdot_tmp[1] = 4.0*p[1]*xB_tmp[1] - 2.0*p[1]*xB_tmp[3];
xBdot_tmp[2] = 2.0*p[1]*xB_tmp[2] - 2.0*p[1]*xB_tmp[5] - (2.0*p[1]*xB_tmp[1])/k[0] + (2.0*p[1]*xB_tmp[3])/k[0] - (2.0*p[1]*xB_tmp[7])/k[0];
xBdot_tmp[3] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[3])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[0] - 4.0*p[1]*xB_tmp[7] - 2.0*p[0]*xB_tmp[1];
xBdot_tmp[4] = 2.0*p[0]*xB_tmp[4] - 2.0*p[0]*xB_tmp[5] - (2.0*p[0]*xB_tmp[6])/k[0] - (2.0*p[0]*xB_tmp[7])/k[0] + (2.0*p[0]*xB_tmp[8])/k[0];
xBdot_tmp[5] = (p[0]*xB_tmp[3])/k[0] - 1.0*p[1]*xB_tmp[4] - (1.0*p[0]*xB_tmp[1])/k[0] - 1.0*p[0]*xB_tmp[2] - (1.0*p[1]*xB_tmp[6])/k[0] + (p[1]*xB_tmp[8])/k[0] + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[5])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[7])/(pow(k[0],2));
xBdot_tmp[6] = 4.0*p[0]*xB_tmp[6] - 2.0*p[0]*xB_tmp[8];
xBdot_tmp[7] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[7])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[8] - 1.0*p[0]*xB_tmp[3];
xBdot_tmp[8] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[8])/(pow(k[0],2)) - 4.0*p[0]*xB_tmp[7] - 2.0*p[1]*xB_tmp[6] - 1.0*p[0]*xB_tmp[0];
xBdot_tmp[9] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[9])/(pow(k[0],2)) - 2.0*p[1]*xB_tmp[17] - 2.0*p[0]*xB_tmp[12];
xBdot_tmp[10] = 4.0*p[1]*xB_tmp[10] - 2.0*p[1]*xB_tmp[12];
xBdot_tmp[11] = 2.0*p[1]*xB_tmp[11] - 2.0*p[1]*xB_tmp[14] - (2.0*p[1]*xB_tmp[10])/k[0] + (2.0*p[1]*xB_tmp[12])/k[0] - (2.0*p[1]*xB_tmp[16])/k[0];
xBdot_tmp[12] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*xB_tmp[12])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[9] - 4.0*p[1]*xB_tmp[16] - 2.0*p[0]*xB_tmp[10];
xBdot_tmp[13] = 2.0*p[0]*xB_tmp[13] - 2.0*p[0]*xB_tmp[14] - (2.0*p[0]*xB_tmp[15])/k[0] - (2.0*p[0]*xB_tmp[16])/k[0] + (2.0*p[0]*xB_tmp[17])/k[0];
xBdot_tmp[14] = (p[0]*xB_tmp[12])/k[0] - 1.0*p[1]*xB_tmp[13] - (1.0*p[0]*xB_tmp[10])/k[0] - 1.0*p[0]*xB_tmp[11] - (1.0*p[1]*xB_tmp[15])/k[0] + (p[1]*xB_tmp[17])/k[0] + ((k[0]*p[0] + k[0]*p[1])*xB_tmp[14])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*xB_tmp[16])/(pow(k[0],2));
xBdot_tmp[15] = 4.0*p[0]*xB_tmp[15] - 2.0*p[0]*xB_tmp[17];
xBdot_tmp[16] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*xB_tmp[16])/(pow(k[0],2)) - 1.0*p[1]*xB_tmp[17] - 1.0*p[0]*xB_tmp[12];
xBdot_tmp[17] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*xB_tmp[17])/(pow(k[0],2)) - 4.0*p[0]*xB_tmp[16] - 2.0*p[1]*xB_tmp[15] - 1.0*p[0]*xB_tmp[9];

  for (ixB=0; ixB<18; ixB++) {
    if(mxIsNaN(xBdot_tmp[ixB])) xBdot_tmp[ixB] = 0.0;
  }

  return(0);
}


 int xQB_enhancer_666(realtype t, N_Vector x, N_Vector xB, N_Vector qBdot, void *user_data)
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
qBdot_tmp[0+ip*ny] = 2.0*x_tmp[4]*xB_tmp[4] - 1.0*x_tmp[5]*xB_tmp[2] + ((2.0*(pow(k[0],2))*x_tmp[0] - 1.0*(pow(k[0],2))*x_tmp[8])*xB_tmp[0])/(pow(k[0],2)) - (1.0*xB_tmp[5]*(2.0*k[0]*x_tmp[4] - 1.0*k[0]*x_tmp[5]))/k[0] - (1.0*xB_tmp[7]*(2.0*k[0]*x_tmp[4] + k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[7] + 4.0*(pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[1])/(pow(k[0],2)) + (xB_tmp[8]*(2.0*k[0]*x_tmp[4] - 2.0*(pow(k[0],2))*x_tmp[6] + 3.0*(pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[3]*(k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[7]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[4] - 4.0*(pow(k[0],2))*x_tmp[6])*xB_tmp[6])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = 2.0*x_tmp[4]*xB_tmp[13] - 1.0*x_tmp[5]*xB_tmp[11] + ((2.0*(pow(k[0],2))*x_tmp[0] - 1.0*(pow(k[0],2))*x_tmp[8])*xB_tmp[9])/(pow(k[0],2)) - (1.0*xB_tmp[14]*(2.0*k[0]*x_tmp[4] - 1.0*k[0]*x_tmp[5]))/k[0] - (1.0*xB_tmp[16]*(2.0*k[0]*x_tmp[4] + k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[7] + 4.0*(pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[10])/(pow(k[0],2)) + (xB_tmp[17]*(2.0*k[0]*x_tmp[4] - 2.0*(pow(k[0],2))*x_tmp[6] + 3.0*(pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) + (xB_tmp[12]*(k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[7]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[4] - 4.0*(pow(k[0],2))*x_tmp[6])*xB_tmp[15])/(pow(k[0],2));

  } break;

  case 1: {
qBdot_tmp[0+ip*ny] = 2.0*x_tmp[2]*xB_tmp[2] - 1.0*x_tmp[5]*xB_tmp[4] + ((2.0*(pow(k[0],2))*x_tmp[0] - 1.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[0])/(pow(k[0],2)) - (1.0*xB_tmp[5]*(2.0*k[0]*x_tmp[2] - 1.0*k[0]*x_tmp[5]))/k[0] - (1.0*xB_tmp[7]*(2.0*k[0]*x_tmp[2] + k[0]*x_tmp[5] + 4.0*(pow(k[0],2))*x_tmp[3] - 2.0*(pow(k[0],2))*x_tmp[7]))/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[8])*xB_tmp[6])/(pow(k[0],2)) + (xB_tmp[3]*(2.0*k[0]*x_tmp[2] - 2.0*(pow(k[0],2))*x_tmp[1] + 3.0*(pow(k[0],2))*x_tmp[3]))/(pow(k[0],2)) + (xB_tmp[8]*(k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[0] - 1.0*(pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[2] - 4.0*(pow(k[0],2))*x_tmp[1])*xB_tmp[1])/(pow(k[0],2));
qBdot_tmp[1+ip*ny] = 2.0*x_tmp[2]*xB_tmp[11] - 1.0*x_tmp[5]*xB_tmp[13] + ((2.0*(pow(k[0],2))*x_tmp[0] - 1.0*(pow(k[0],2))*x_tmp[3])*xB_tmp[9])/(pow(k[0],2)) - (1.0*xB_tmp[14]*(2.0*k[0]*x_tmp[2] - 1.0*k[0]*x_tmp[5]))/k[0] - (1.0*xB_tmp[16]*(2.0*k[0]*x_tmp[2] + k[0]*x_tmp[5] + 4.0*(pow(k[0],2))*x_tmp[3] - 2.0*(pow(k[0],2))*x_tmp[7]))/(pow(k[0],2)) - (1.0*(k[0]*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[8])*xB_tmp[15])/(pow(k[0],2)) + (xB_tmp[12]*(2.0*k[0]*x_tmp[2] - 2.0*(pow(k[0],2))*x_tmp[1] + 3.0*(pow(k[0],2))*x_tmp[3]))/(pow(k[0],2)) + (xB_tmp[17]*(k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[0] - 1.0*(pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*(2.0*k[0]*x_tmp[2] - 4.0*(pow(k[0],2))*x_tmp[1])*xB_tmp[10])/(pow(k[0],2));

  } break;

  }
  }

  for (iyp=0; iyp<2*np; iyp++) {
    if(mxIsNaN(qBdot_tmp[iyp])) qBdot_tmp[iyp] = 0.0;
  }

  return(0);
}


 void x0_enhancer_666(N_Vector x0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x0_tmp = N_VGetArrayPointer(x0);
  memset(x0_tmp,0,sizeof(double)*9);
x0_tmp[0] = (k[6]*k[15])/(pow(k[0],2));
x0_tmp[1] = (k[9]*k[18])/(pow(k[0],2));
x0_tmp[2] = (k[3]*k[12])/k[0];
x0_tmp[3] = (k[8]*k[17])/(pow(k[0],2));
x0_tmp[4] = (k[1]*k[10] - 1.0*k[1] + 1.0)/k[0];
x0_tmp[5] = (k[2]*k[11])/k[0];
x0_tmp[6] = (k[4]*k[13])/(pow(k[0],2));
x0_tmp[7] = (k[7]*k[16])/(pow(k[0],2));
x0_tmp[8] = (k[5]*k[14])/(pow(k[0],2));
  
  
  return;
}


 int Jv_enhancer_666(N_Vector v, N_Vector Jv, realtype t,
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
  memset(Jv_tmp,0,sizeof(double)*9);
Jv_tmp[0] = p[1]*v_tmp[3] + p[0]*v_tmp[8] - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[0])/(pow(k[0],2));
Jv_tmp[1] = 2.0*p[0]*v_tmp[3] - 4.0*p[1]*v_tmp[1] + (2.0*p[1]*v_tmp[2])/k[0] + (p[0]*v_tmp[5])/k[0];
Jv_tmp[2] = p[0]*v_tmp[5] - 2.0*p[1]*v_tmp[2];
Jv_tmp[3] = 2.0*p[0]*v_tmp[0] + 2.0*p[1]*v_tmp[1] + p[0]*v_tmp[7] - (2.0*p[1]*v_tmp[2])/k[0] - (1.0*p[0]*v_tmp[5])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*v_tmp[3])/(pow(k[0],2));
Jv_tmp[4] = p[1]*v_tmp[5] - 2.0*p[0]*v_tmp[4];
Jv_tmp[5] = 2.0*p[1]*v_tmp[2] + 2.0*p[0]*v_tmp[4] - (1.0*(k[0]*p[0] + k[0]*p[1])*v_tmp[5])/k[0];
Jv_tmp[6] = 2.0*p[1]*v_tmp[8] - 4.0*p[0]*v_tmp[6] + (2.0*p[0]*v_tmp[4])/k[0] + (p[1]*v_tmp[5])/k[0];
Jv_tmp[7] = 4.0*p[1]*v_tmp[3] + 4.0*p[0]*v_tmp[8] + (2.0*p[1]*v_tmp[2])/k[0] + (2.0*p[0]*v_tmp[4])/k[0] + ((k[0]*p[0] + k[0]*p[1])*v_tmp[5])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*v_tmp[7])/(pow(k[0],2));
Jv_tmp[8] = 2.0*p[1]*v_tmp[0] + 2.0*p[0]*v_tmp[6] + p[1]*v_tmp[7] - (2.0*p[0]*v_tmp[4])/k[0] - (1.0*p[1]*v_tmp[5])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*v_tmp[8])/(pow(k[0],2));

  for (ix=0; ix<9; ix++) {
    if(mxIsNaN(Jv_tmp[ix])) Jv_tmp[ix] = 0.0;
  }

  return(0);
}
 int JvB_enhancer_666(N_Vector vB, N_Vector JvB, realtype t,
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
  memset(JvB_tmp,0,sizeof(double)*9);
JvB_tmp[0] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[0])/(pow(k[0],2)) - 2.0*p[1]*vB_tmp[8] - 2.0*p[0]*vB_tmp[3];
JvB_tmp[1] = 4.0*p[1]*vB_tmp[1] - 2.0*p[1]*vB_tmp[3];
JvB_tmp[2] = 2.0*p[1]*vB_tmp[2] - 2.0*p[1]*vB_tmp[5] - (2.0*p[1]*vB_tmp[1])/k[0] + (2.0*p[1]*vB_tmp[3])/k[0] - (2.0*p[1]*vB_tmp[7])/k[0];
JvB_tmp[3] = (((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*vB_tmp[3])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[0] - 4.0*p[1]*vB_tmp[7] - 2.0*p[0]*vB_tmp[1];
JvB_tmp[4] = 2.0*p[0]*vB_tmp[4] - 2.0*p[0]*vB_tmp[5] - (2.0*p[0]*vB_tmp[6])/k[0] - (2.0*p[0]*vB_tmp[7])/k[0] + (2.0*p[0]*vB_tmp[8])/k[0];
JvB_tmp[5] = (p[0]*vB_tmp[3])/k[0] - 1.0*p[1]*vB_tmp[4] - (1.0*p[0]*vB_tmp[1])/k[0] - 1.0*p[0]*vB_tmp[2] - (1.0*p[1]*vB_tmp[6])/k[0] + (p[1]*vB_tmp[8])/k[0] + ((k[0]*p[0] + k[0]*p[1])*vB_tmp[5])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*vB_tmp[7])/(pow(k[0],2));
JvB_tmp[6] = 4.0*p[0]*vB_tmp[6] - 2.0*p[0]*vB_tmp[8];
JvB_tmp[7] = ((2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*vB_tmp[7])/(pow(k[0],2)) - 1.0*p[1]*vB_tmp[8] - 1.0*p[0]*vB_tmp[3];
JvB_tmp[8] = ((3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*vB_tmp[8])/(pow(k[0],2)) - 4.0*p[0]*vB_tmp[7] - 2.0*p[1]*vB_tmp[6] - 1.0*p[0]*vB_tmp[0];

  for (ix=0; ix<9; ix++) {
    if(mxIsNaN(JvB_tmp[ix])) JvB_tmp[ix] = 0.0;
  }

  return(0);
}


 int JBand_enhancer_666(long int N, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xdot,
  	DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(J_enhancer_666(N,t,x,xdot,J,user_data,tmp1,tmp2,tmp3));
}


 int J_enhancer_666(long int N, realtype t, N_Vector x,
  	N_Vector xdot, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  int iJ;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  memset(J->data,0,sizeof(double)*81);
J->data[0] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[3] = 2.0*p[0];
J->data[8] = 2.0*p[1];
J->data[10] = -4.0*p[1];
J->data[12] = 2.0*p[1];
J->data[19] = (2.0*p[1])/k[0];
J->data[20] = -2.0*p[1];
J->data[21] = -(2.0*p[1])/k[0];
J->data[23] = 2.0*p[1];
J->data[25] = (2.0*p[1])/k[0];
J->data[27] = p[1];
J->data[28] = 2.0*p[0];
J->data[30] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[34] = 4.0*p[1];
J->data[40] = -2.0*p[0];
J->data[41] = 2.0*p[0];
J->data[42] = (2.0*p[0])/k[0];
J->data[43] = (2.0*p[0])/k[0];
J->data[44] = -(2.0*p[0])/k[0];
J->data[46] = p[0]/k[0];
J->data[47] = p[0];
J->data[48] = -(1.0*p[0])/k[0];
J->data[49] = p[1];
J->data[50] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];
J->data[51] = p[1]/k[0];
J->data[52] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[53] = -(1.0*p[1])/k[0];
J->data[60] = -4.0*p[0];
J->data[62] = 2.0*p[0];
J->data[66] = p[0];
J->data[70] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[71] = p[1];
J->data[72] = p[0];
J->data[78] = 2.0*p[1];
J->data[79] = 4.0*p[0];
J->data[80] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));

  for (iJ=0; iJ<81; iJ++) {
    if(mxIsNaN(J->data[iJ])) J->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparse_enhancer_666(realtype t, N_Vector x,
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
  J->rowvals[1] = 3;
  J->rowvals[2] = 8;
  J->rowvals[3] = 1;
  J->rowvals[4] = 3;
  J->rowvals[5] = 1;
  J->rowvals[6] = 2;
  J->rowvals[7] = 3;
  J->rowvals[8] = 5;
  J->rowvals[9] = 7;
  J->rowvals[10] = 0;
  J->rowvals[11] = 1;
  J->rowvals[12] = 3;
  J->rowvals[13] = 7;
  J->rowvals[14] = 4;
  J->rowvals[15] = 5;
  J->rowvals[16] = 6;
  J->rowvals[17] = 7;
  J->rowvals[18] = 8;
  J->rowvals[19] = 1;
  J->rowvals[20] = 2;
  J->rowvals[21] = 3;
  J->rowvals[22] = 4;
  J->rowvals[23] = 5;
  J->rowvals[24] = 6;
  J->rowvals[25] = 7;
  J->rowvals[26] = 8;
  J->rowvals[27] = 6;
  J->rowvals[28] = 8;
  J->rowvals[29] = 3;
  J->rowvals[30] = 7;
  J->rowvals[31] = 8;
  J->rowvals[32] = 0;
  J->rowvals[33] = 6;
  J->rowvals[34] = 7;
  J->rowvals[35] = 8;
  J->colptrs[0] = 0;
  J->colptrs[1] = 3;
  J->colptrs[2] = 5;
  J->colptrs[3] = 10;
  J->colptrs[4] = 14;
  J->colptrs[5] = 19;
  J->colptrs[6] = 27;
  J->colptrs[7] = 29;
  J->colptrs[8] = 32;
  J->colptrs[9] = 36;
J->data[0] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[1] = 2.0*p[0];
J->data[2] = 2.0*p[1];
J->data[3] = -4.0*p[1];
J->data[4] = 2.0*p[1];
J->data[5] = (2.0*p[1])/k[0];
J->data[6] = -2.0*p[1];
J->data[7] = -(2.0*p[1])/k[0];
J->data[8] = 2.0*p[1];
J->data[9] = (2.0*p[1])/k[0];
J->data[10] = p[1];
J->data[11] = 2.0*p[0];
J->data[12] = -(1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[13] = 4.0*p[1];
J->data[14] = -2.0*p[0];
J->data[15] = 2.0*p[0];
J->data[16] = (2.0*p[0])/k[0];
J->data[17] = (2.0*p[0])/k[0];
J->data[18] = -(2.0*p[0])/k[0];
J->data[19] = p[0]/k[0];
J->data[20] = p[0];
J->data[21] = -(1.0*p[0])/k[0];
J->data[22] = p[1];
J->data[23] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/k[0];
J->data[24] = p[1]/k[0];
J->data[25] = (k[0]*p[0] + k[0]*p[1])/(pow(k[0],2));
J->data[26] = -(1.0*p[1])/k[0];
J->data[27] = -4.0*p[0];
J->data[28] = 2.0*p[0];
J->data[29] = p[0];
J->data[30] = -(1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1]))/(pow(k[0],2));
J->data[31] = p[1];
J->data[32] = p[0];
J->data[33] = 2.0*p[1];
J->data[34] = 4.0*p[0];
J->data[35] = -(1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1]))/(pow(k[0],2));
  return(0);
}


 int JBBand_enhancer_666(long int NeqB, long int mupper, long int mlower,   realtype t, N_Vector x, N_Vector xB,
  	N_Vector xdotB, DlsMat J, void *user_data, 
  	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  return(JB_enhancer_666(NeqB,t,x,xB,xdotB,J,user_data,tmp1,tmp2,tmp3));
}
 int JB_enhancer_666(long int N, realtype t, N_Vector x,
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
JB->data[0] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[3] = -1.0*p[1];
JB->data[8] = -1.0*p[0];
JB->data[10] = 4.0*p[1];
JB->data[11] = -(2.0*p[1])/k[0];
JB->data[12] = -2.0*p[0];
JB->data[14] = -(1.0*p[0])/k[0];
JB->data[20] = 2.0*p[1];
JB->data[23] = -1.0*p[0];
JB->data[27] = -2.0*p[0];
JB->data[28] = -2.0*p[1];
JB->data[29] = (2.0*p[1])/k[0];
JB->data[30] = ((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[32] = p[0]/k[0];
JB->data[34] = -1.0*p[0];
JB->data[40] = 2.0*p[0];
JB->data[41] = -1.0*p[1];
JB->data[47] = -2.0*p[1];
JB->data[49] = -2.0*p[0];
JB->data[50] = (k[0]*p[0] + k[0]*p[1])/k[0];
JB->data[58] = -(2.0*p[0])/k[0];
JB->data[59] = -(1.0*p[1])/k[0];
JB->data[60] = 4.0*p[0];
JB->data[62] = -2.0*p[1];
JB->data[65] = -(2.0*p[1])/k[0];
JB->data[66] = -4.0*p[1];
JB->data[67] = -(2.0*p[0])/k[0];
JB->data[68] = -(1.0*(k[0]*p[0] + k[0]*p[1]))/(pow(k[0],2));
JB->data[70] = (2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])/(pow(k[0],2));
JB->data[71] = -4.0*p[0];
JB->data[72] = -2.0*p[1];
JB->data[76] = (2.0*p[0])/k[0];
JB->data[77] = p[1]/k[0];
JB->data[78] = -2.0*p[0];
JB->data[79] = -1.0*p[1];
JB->data[80] = (3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])/(pow(k[0],2));

  for (iJ=0; iJ<81; iJ++) {
    if(mxIsNaN(JB->data[iJ])) JB->data[iJ] = 0.0;
  }

  return(0);
}


 int JSparseB_enhancer_666(realtype t, N_Vector x,
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
  JB->rowvals[1] = 3;
  JB->rowvals[2] = 8;
  JB->rowvals[3] = 1;
  JB->rowvals[4] = 2;
  JB->rowvals[5] = 3;
  JB->rowvals[6] = 5;
  JB->rowvals[7] = 2;
  JB->rowvals[8] = 5;
  JB->rowvals[9] = 0;
  JB->rowvals[10] = 1;
  JB->rowvals[11] = 2;
  JB->rowvals[12] = 3;
  JB->rowvals[13] = 5;
  JB->rowvals[14] = 7;
  JB->rowvals[15] = 4;
  JB->rowvals[16] = 5;
  JB->rowvals[17] = 2;
  JB->rowvals[18] = 4;
  JB->rowvals[19] = 5;
  JB->rowvals[20] = 4;
  JB->rowvals[21] = 5;
  JB->rowvals[22] = 6;
  JB->rowvals[23] = 8;
  JB->rowvals[24] = 2;
  JB->rowvals[25] = 3;
  JB->rowvals[26] = 4;
  JB->rowvals[27] = 5;
  JB->rowvals[28] = 7;
  JB->rowvals[29] = 8;
  JB->rowvals[30] = 0;
  JB->rowvals[31] = 4;
  JB->rowvals[32] = 5;
  JB->rowvals[33] = 6;
  JB->rowvals[34] = 7;
  JB->rowvals[35] = 8;
  JB->colptrs[0] = 0;
  JB->colptrs[1] = 3;
  JB->colptrs[2] = 7;
  JB->colptrs[3] = 9;
  JB->colptrs[4] = 15;
  JB->colptrs[5] = 17;
  JB->colptrs[6] = 20;
  JB->colptrs[7] = 24;
  JB->colptrs[8] = 30;
  JB->colptrs[9] = 36;
  return(0);
}


 int sx_enhancer_666(int Ns, realtype t, N_Vector x, N_Vector xdot,
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
  memset(sxdot_tmp,0,sizeof(double)*9);
  switch (plist[ip]) {
  case 0: {
sxdot_tmp[0] = p[1]*sx_tmp[3] + p[0]*sx_tmp[8] - (1.0*(2.0*(pow(k[0],2))*x_tmp[0] - 1.0*(pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[0])/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[0]*sx_tmp[3] - 4.0*p[1]*sx_tmp[1] + (k[0]*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[3])/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[2])/k[0] + (p[0]*sx_tmp[5])/k[0];
sxdot_tmp[2] = p[0]*sx_tmp[5] - 2.0*p[1]*sx_tmp[2] + x_tmp[5];
sxdot_tmp[3] = 2.0*p[0]*sx_tmp[0] + 2.0*p[1]*sx_tmp[1] + p[0]*sx_tmp[7] - (1.0*(k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[0] + (pow(k[0],2))*x_tmp[3] - 1.0*(pow(k[0],2))*x_tmp[7]))/(pow(k[0],2)) - (2.0*p[1]*sx_tmp[2])/k[0] - (1.0*p[0]*sx_tmp[5])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[5] - 2.0*p[0]*sx_tmp[4] - 2.0*x_tmp[4];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[2] + 2.0*p[0]*sx_tmp[4] + (2.0*k[0]*x_tmp[4] - 1.0*k[0]*x_tmp[5])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[5])/k[0];
sxdot_tmp[6] = 2.0*p[1]*sx_tmp[8] - 4.0*p[0]*sx_tmp[6] + (2.0*k[0]*x_tmp[4] - 4.0*(pow(k[0],2))*x_tmp[6])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[4])/k[0] + (p[1]*sx_tmp[5])/k[0];
sxdot_tmp[7] = 4.0*p[1]*sx_tmp[3] + 4.0*p[0]*sx_tmp[8] + (2.0*k[0]*x_tmp[4] + k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[7] + 4.0*(pow(k[0],2))*x_tmp[8])/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[2])/k[0] + (2.0*p[0]*sx_tmp[4])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[5])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[7])/(pow(k[0],2));
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[0] + 2.0*p[0]*sx_tmp[6] + p[1]*sx_tmp[7] - (1.0*(2.0*k[0]*x_tmp[4] - 2.0*(pow(k[0],2))*x_tmp[6] + 3.0*(pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (2.0*p[0]*sx_tmp[4])/k[0] - (1.0*p[1]*sx_tmp[5])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[8])/(pow(k[0],2));

  } break;

  case 1: {
sxdot_tmp[0] = p[1]*sx_tmp[3] + p[0]*sx_tmp[8] - (1.0*(2.0*(pow(k[0],2))*x_tmp[0] - 1.0*(pow(k[0],2))*x_tmp[3]))/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[0])/(pow(k[0],2));
sxdot_tmp[1] = 2.0*p[0]*sx_tmp[3] - 4.0*p[1]*sx_tmp[1] + (2.0*k[0]*x_tmp[2] - 4.0*(pow(k[0],2))*x_tmp[1])/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[2])/k[0] + (p[0]*sx_tmp[5])/k[0];
sxdot_tmp[2] = p[0]*sx_tmp[5] - 2.0*p[1]*sx_tmp[2] - 2.0*x_tmp[2];
sxdot_tmp[3] = 2.0*p[0]*sx_tmp[0] + 2.0*p[1]*sx_tmp[1] + p[0]*sx_tmp[7] - (1.0*(2.0*k[0]*x_tmp[2] - 2.0*(pow(k[0],2))*x_tmp[1] + 3.0*(pow(k[0],2))*x_tmp[3]))/(pow(k[0],2)) - (2.0*p[1]*sx_tmp[2])/k[0] - (1.0*p[0]*sx_tmp[5])/k[0] - (1.0*((pow(k[0],2))*p[0] + 3.0*(pow(k[0],2))*p[1])*sx_tmp[3])/(pow(k[0],2));
sxdot_tmp[4] = p[1]*sx_tmp[5] - 2.0*p[0]*sx_tmp[4] + x_tmp[5];
sxdot_tmp[5] = 2.0*p[1]*sx_tmp[2] + 2.0*p[0]*sx_tmp[4] + (2.0*k[0]*x_tmp[2] - 1.0*k[0]*x_tmp[5])/k[0] - (1.0*(k[0]*p[0] + k[0]*p[1])*sx_tmp[5])/k[0];
sxdot_tmp[6] = 2.0*p[1]*sx_tmp[8] - 4.0*p[0]*sx_tmp[6] + (k[0]*x_tmp[5] + 2.0*(pow(k[0],2))*x_tmp[8])/(pow(k[0],2)) + (2.0*p[0]*sx_tmp[4])/k[0] + (p[1]*sx_tmp[5])/k[0];
sxdot_tmp[7] = 4.0*p[1]*sx_tmp[3] + 4.0*p[0]*sx_tmp[8] + (2.0*k[0]*x_tmp[2] + k[0]*x_tmp[5] + 4.0*(pow(k[0],2))*x_tmp[3] - 2.0*(pow(k[0],2))*x_tmp[7])/(pow(k[0],2)) + (2.0*p[1]*sx_tmp[2])/k[0] + (2.0*p[0]*sx_tmp[4])/k[0] + ((k[0]*p[0] + k[0]*p[1])*sx_tmp[5])/(pow(k[0],2)) - (1.0*(2.0*(pow(k[0],2))*p[0] + 2.0*(pow(k[0],2))*p[1])*sx_tmp[7])/(pow(k[0],2));
sxdot_tmp[8] = 2.0*p[1]*sx_tmp[0] + 2.0*p[0]*sx_tmp[6] + p[1]*sx_tmp[7] - (1.0*(k[0]*x_tmp[5] - 2.0*(pow(k[0],2))*x_tmp[0] - 1.0*(pow(k[0],2))*x_tmp[7] + (pow(k[0],2))*x_tmp[8]))/(pow(k[0],2)) - (2.0*p[0]*sx_tmp[4])/k[0] - (1.0*p[1]*sx_tmp[5])/k[0] - (1.0*(3.0*(pow(k[0],2))*p[0] + (pow(k[0],2))*p[1])*sx_tmp[8])/(pow(k[0],2));

  } break;

  }
 for (ix=0; ix<9; ix++) {
    if(mxIsNaN(sxdot_tmp[ix])) sxdot_tmp[ix] = 0.0;
  }

  return(0);
}


 void sx0_enhancer_666(int ip, N_Vector sx0, void *user_data)
{
  UserData data = (UserData) user_data;
  double *p = data->p;
  int *plist = data->plist;
  double *k = data->k;
  double *u = data->u;
  double *sx0_tmp = N_VGetArrayPointer(sx0);
  memset(sx0_tmp,0,sizeof(double)*9);
  switch (ip) {
  }

  return;
}


void y_enhancer_666(double t, int nt, int it, double *y, double *p, double *k, double *u, double *x){
y[it+nt*0] = x[it+nt*4];
y[it+nt*1] = x[it+nt*6];
    
    return;
}


void dydp_enhancer_666(double t, int nt, int it, double *dydp, double *y, double *p, double *k, double *u, double *x, int *plist, int np, int ny){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  }
  }
  
  return;
}


void dydx_enhancer_666(double t,double *dydx, double *y, double *p, double *k, double *x){
  memset(dydx,0,sizeof(double)*18);
dydx[8] = 1.0;
dydx[13] = 1.0;
  
  return;
}


void sy_enhancer_666(double t, int nt, int it, int ip, int np, int nx, int ny, double *sy, double *p, double *k, double *x, double *sx){
  switch (ip) {
  case 0: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(6+np*nx)];

  } break;

  case 1: {
sy[it+nt*(0+np*ny)] = sx[it+nt*(4+np*nx)];
sy[it+nt*(1+np*ny)] = sx[it+nt*(6+np*nx)];

  } break;

  }
  
  return;
}
int root_enhancer_666(double t, N_Vector x, realtype *gout, void *user_data){
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  return(0);
}
double sroot_enhancer_666(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2root_enhancer_666(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double srootval_enhancer_666(double t, int ip, int ir, N_Vector x, N_Vector sx, void *user_data){
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
double s2rootval_enhancer_666(double t, int ip, int jp, int ir, N_Vector x, N_Vector sx, void *user_data){
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
void deltadisc_enhancer_666(double t, int idisc, N_Vector x, void *user_data){
  int ix;
  UserData data = (UserData) user_data;
  double *p = data->p;
  double *k = data->k;
  double *u = data->u;
  double *x_tmp = N_VGetArrayPointer(x);
  double deltadisc[9];
  memset(deltadisc,0,sizeof(double)*9);
  for(ix = 0; ix<9;ix++){;
  x_tmp[ix] += deltadisc[ix];
  };
}
void sdeltadisc_enhancer_666(double t, int idisc, N_Vector x, N_Vector *sx, void *user_data){
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
  double deltadisc[9];
  double *sdeltadisc;
  memset(deltadisc,0,sizeof(double)*9);
  sdeltadisc = mxMalloc(sizeof(double)*9*np);
  memset(sdeltadisc,0,sizeof(double)*9*np);
  for (ip=0; ip<np; ip++) {
  sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
     switch (plist[ip]) {
     }
  }
  for(ip = 0; ip<np;ip++){
      sx_tmp = N_VGetArrayPointer(sx[plist[ip]]);
      for(ix = 0; ix<9;ix++){
      sx_tmp[ix] += sdeltadisc[plist[ip]+np*ix];
     }
  }
  for(ix = 0; ix<9;ix++){
  x_tmp[ix] += deltadisc[ix];
  };
 mxFree(sdeltadisc);
}


void dxdotdp_enhancer_666(double t, int nt, int it, double *dxdotdp, double *p, double *k, double *u, double *x, int *plist, int np, int nx){
  
  int ip;
  for(ip=0; ip<np; ip++) {
  switch (plist[ip]) {
  case 0: {
dxdotdp[(0+ip*nx)] = -(1.0*(2.0*(pow(k[0],2))*x[it+nt*0] - 1.0*(pow(k[0],2))*x[it+nt*8]))/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = (k[0]*x[it+nt*5] + 2.0*(pow(k[0],2))*x[it+nt*3])/(pow(k[0],2));
dxdotdp[(2+ip*nx)] = x[it+nt*5];
dxdotdp[(3+ip*nx)] = -(1.0*(k[0]*x[it+nt*5] - 2.0*(pow(k[0],2))*x[it+nt*0] + (pow(k[0],2))*x[it+nt*3] - 1.0*(pow(k[0],2))*x[it+nt*7]))/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = -2.0*x[it+nt*4];
dxdotdp[(5+ip*nx)] = (2.0*k[0]*x[it+nt*4] - 1.0*k[0]*x[it+nt*5])/k[0];
dxdotdp[(6+ip*nx)] = (2.0*k[0]*x[it+nt*4] - 4.0*(pow(k[0],2))*x[it+nt*6])/(pow(k[0],2));
dxdotdp[(7+ip*nx)] = (2.0*k[0]*x[it+nt*4] + k[0]*x[it+nt*5] - 2.0*(pow(k[0],2))*x[it+nt*7] + 4.0*(pow(k[0],2))*x[it+nt*8])/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = -(1.0*(2.0*k[0]*x[it+nt*4] - 2.0*(pow(k[0],2))*x[it+nt*6] + 3.0*(pow(k[0],2))*x[it+nt*8]))/(pow(k[0],2));

  } break;

  case 1: {
dxdotdp[(0+ip*nx)] = -(1.0*(2.0*(pow(k[0],2))*x[it+nt*0] - 1.0*(pow(k[0],2))*x[it+nt*3]))/(pow(k[0],2));
dxdotdp[(1+ip*nx)] = (2.0*k[0]*x[it+nt*2] - 4.0*(pow(k[0],2))*x[it+nt*1])/(pow(k[0],2));
dxdotdp[(2+ip*nx)] = -2.0*x[it+nt*2];
dxdotdp[(3+ip*nx)] = -(1.0*(2.0*k[0]*x[it+nt*2] - 2.0*(pow(k[0],2))*x[it+nt*1] + 3.0*(pow(k[0],2))*x[it+nt*3]))/(pow(k[0],2));
dxdotdp[(4+ip*nx)] = x[it+nt*5];
dxdotdp[(5+ip*nx)] = (2.0*k[0]*x[it+nt*2] - 1.0*k[0]*x[it+nt*5])/k[0];
dxdotdp[(6+ip*nx)] = (k[0]*x[it+nt*5] + 2.0*(pow(k[0],2))*x[it+nt*8])/(pow(k[0],2));
dxdotdp[(7+ip*nx)] = (2.0*k[0]*x[it+nt*2] + k[0]*x[it+nt*5] + 4.0*(pow(k[0],2))*x[it+nt*3] - 2.0*(pow(k[0],2))*x[it+nt*7])/(pow(k[0],2));
dxdotdp[(8+ip*nx)] = -(1.0*(k[0]*x[it+nt*5] - 2.0*(pow(k[0],2))*x[it+nt*0] - 1.0*(pow(k[0],2))*x[it+nt*7] + (pow(k[0],2))*x[it+nt*8]))/(pow(k[0],2));

  } break;

  }
  }
  
  return;
}
