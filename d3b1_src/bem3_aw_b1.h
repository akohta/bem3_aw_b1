/*
 * d3b1_const.h
 *
 *  Created on: Nov 10, 2021
 *      Author: ohta
 */
#ifndef BEM3_AW_B1_H_
#define BEM3_AW_B1_H_

#include "d3b1_elem.h"

// -- bem3_aw_b1.c --
void read_dmda(int argc,char **argv,DMDA *ad); // read datafile
void print_dmda(DMDA *ad);                     // print datafile 
void initialize_dmda(DMDA *ad);                // initialize coefficient
void finalize_dmda(DMDA *ad);                  // free allocated memory
int domain_id_m_dmda(double *rt,DMDA *ad);     // return domain id of point rt, return the main domain id on boundary
int domain_id_s_dmda(double *rt,DMDA *ad);     // return domain id of point rt, return the  sub domain id on boundary
void dat_write_dmda(char *filename,DMDA *ad);  // write datafile
void dat_read_dmda(char *filename,DMDA *ad);   // reat datafile 

// -- bem3_aw_b1_solve_bieq.c --
void solve_bieq_dmda(DMDA *ad); // solve boundary integral equations


// -- bem3_aw_b1_force.c --
// calculation of net radiation force and torque
void force_FN(double *F,double *N,double *rc,int type,DMDA *ad);
// outputs
// F : radiation force,  F[0]=F_x, F[1]=F_y, F[2]=F_z.
// N : radiation torque, N[0]=N_x, N[1]=N_y, N[2]=N_z.
// inputs
// rc : coordinate of rotation center, rc[0]=rc_x, rc[1]=rc_y, rc[2]=rc_z.
// type : setting of numerical integration, type=0:4-point GL, type!=0:9-point or 7-point GL.
// ad : pointer of DMDA object.

void surface_area(double *S,int did,int type,DMDA *ad);
// output
// S : surface area of domain id "did".
// did : domain id.
// type : setting of numerical integration, type=0:4-point GL, type!=0:9-point or 7-point GL.
// ad : pointer of DMDA object.


// -- bem3_aw_b1_force.c --
// velocity potential
int phi_s_dmda(double complex *phi,double *rt,int type,DMDA *ad); // scattered or internal field
int phi_t_dmda(double complex *phi,double *rt,int type,DMDA *ad); // total(incident + scattered) field
int phi_i_dmda(double complex *phi,double *rt,int type,DMDA *ad); // incident field
// output
// phi : velocity potential.
// intputs
// rt : coordinate of calcluation point, rt[0]=x, rt[1]=y, rt[2]=z.
// type : setting of numerical integration, type=0:4pGL,1:9pGL,2:GLN-point GL,3:GHN-point GL,4:DE.
// ad : pointer of DMDA object.
// return domain id.

// sound pressure
int p_s_dmda(double complex *p,double *rt,int type,DMDA *ad); // scattered or internal field
int p_t_dmda(double complex *p,double *rt,int type,DMDA *ad); // total field
int p_i_dmda(double complex *p,double *rt,int type,DMDA *ad); // incident field
// output
// p : sound pressure.
// others are the same as phi_*_dmda().

// sound pressure and particle velocity
int pv_s_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // scattered or internal field
int pv_t_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // total field
int pv_i_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // incident field
int pv_s_dbieq_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // for far-field (using derivative boundary integral eq.)
int pv_t_dbieq_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // 
int pv_i_dbieq_dmda(double complex *p,double complex *pv,double *rt,int type,DMDA *ad); // same as pv_i_dmda()
// outputs
// p : sound pressure.
// pv : particle velocity.
// others are the same as p_*_dmda().

// did: domain id, t:element id of each domain, zeta_t,eta_t: parameter on surface (main domain coordinate)
// velocity potential on the boundary
void phi_s_bd_dmda(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad);
void phi_t_bd_dmda(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad);
void phi_i_bd_dmda(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad);
// output
// phi : velocity potential.
// inputs 
// did : domain id.
// t : element id.
// zeta_t, eta_t : parameter of the calculation point on the boundary (main domain coordinate), -1 < zeta_t < 1, -1 < eta_t < 1
// type : setting of numerical integration, type=0:4pGL,1:9pGL,2:GLN-point GL,3:GHN-point GL,4:DE.

// sound pressure and particle velocity on the boundary
void pv_s_bd_dmda(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad);
void pv_t_bd_dmda(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad);
void pv_i_bd_dmda(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,DMDA *ad);
// outputs
// p : sound pressure.
// pv : particle velocity.
// others are the same as phi_*_bd_dmda(). 

#endif
