#include "bem3_aw_b1.h"


void force_FN(double *F,double *N,double *rc,int type,DMDA *ad)
{
  void bil_force_4p(double *F,double *N,double *rc,int s,DMDA *ad);
  void bil_force_9p(double *F,double *N,double *rc,int s,DMDA *ad);
  void lit_force_4p(double *F,double *N,double *rc,int s,DMDA *ad);
  void lit_force_7p(double *F,double *N,double *rc,int s,DMDA *ad);
  
  double tf[3],tn[3],Fx,Fy,Fz,Nx,Ny,Nz;
  int t,td;

  Fx=0.0;    Fy=0.0;    Fz=0.0;
  Nx=0.0;    Ny=0.0;    Nz=0.0;
  
  #pragma omp parallel for schedule(dynamic) reduction(+:Fx,Fy,Fz,Nx,Ny,Nz) private(td,tf,tn)
  for(t=1;t<=ad->bd.sb[0].Ne;t++){
    td=ad->bd.sb[0].sid[t];
    if( ELT4==check_element_type(td,&(ad->bd)) ){
      if(type==0) bil_force_4p(tf,tn,rc,t,ad);
      else        bil_force_9p(tf,tn,rc,t,ad);
    }
    else {
      if(type==0) lit_force_4p(tf,tn,rc,t,ad);
      else        lit_force_7p(tf,tn,rc,t,ad);
    }
    Fx+=tf[0];
    Fy+=tf[1];
    Fz+=tf[2];
    Nx+=tn[0];
    Ny+=tn[1];
    Nz+=tn[2];
  }

  F[0]=Fx;
  F[1]=Fy;
  F[2]=Fz;
  N[0]=Nx;
  N[1]=Ny;
  N[2]=Nz;  
}

void surface_area(double *S,int did,int type,DMDA *ad)
{
  void bil_sa_4p(double *S,int d,int s,DMDA *ad);
  void bil_sa_9p(double *S,int d,int s,DMDA *ad);
  void lit_sa_4p(double *S,int d,int s,DMDA *ad);
  void lit_sa_7p(double *S,int d,int s,DMDA *ad);

  double ts,sa;
  int t,td;

  sa=0.0;
  for(t=1;t<=ad->bd.sb[did].Ne;t++){
    td=ad->bd.sb[did].sid[t];

    if( ELT4==check_element_type(td,&(ad->bd)) ){
      if(type==0) bil_sa_4p(&ts,did,t,ad);
      else bil_sa_9p(&ts,did,t,ad);
    }
    else {
      if(type==0) lit_sa_4p(&ts,did,t,ad);
      else lit_sa_7p(&ts,did,t,ad);
    }
    sa+=ts;
  }
  *S=sa;
  
}

////////////////////////////////////////////////////////////////////
void bil_force_4p(double *F,double *N,double *rc,int s,DMDA *ad)
{
  double complex p,v[3];
  double r[3],w[3],T[9],c1,c2,c3,t1;
  int sd,asd,i,j,l,m;

  c1=0.25/(ad->aw.c0*ad->aw.c0*ad->aw.rho0); // 1/(4*K0)
  c2=0.25*ad->aw.rho0; // rho0/4
  c3=0.50*ad->aw.rho0; // rho0/2

  sd=ad->bd.sb[0].sid[s];
  asd=abs(sd);
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<4;i++){
    for(j=0;j<3;j++){
      r[j]=ad->bd.ren[asd][i][j];
      w[j]=ad->bd.wen[asd][i][j];
    }
    calc_maw_pv(&p,v,r,&(ad->aw)); // multi_aw.h
    p+=-ad->bd.sb[0].P[s][i]*ad->k2[0];
    for(l=0;l<3;l++) v[l]+=ad->bd.sb[0].pv[s][i][l];
    // radiation stress tensor
    t1=c1*creal(p*conj(p))-c2*creal(v[0]*conj(v[0])+v[1]*conj(v[1])+v[2]*conj(v[2]));
    for(l=0;l<3;l++) for(m=0;m<3;m++) T[l*3+m]=c3*creal(v[l]*conj(v[m]));
    T[0*3+0]+=t1;
    T[1*3+1]+=t1;
    T[2*3+2]+=t1; 
    // force
    F[0]+=(T[0*3+0]*w[0]+T[0*3+1]*w[1]+T[0*3+2]*w[2])*ad->bd.wt_44[i];
    F[1]+=(T[1*3+0]*w[0]+T[1*3+1]*w[1]+T[1*3+2]*w[2])*ad->bd.wt_44[i];
    F[2]+=(T[2*3+0]*w[0]+T[2*3+1]*w[1]+T[2*3+2]*w[2])*ad->bd.wt_44[i]; 
    // torque
    N[0]+=( ((r[1]-rc[1])*T[2*3+0]-(r[2]-rc[2])*T[1*3+0])*w[0]
           +((r[1]-rc[1])*T[2*3+1]-(r[2]-rc[2])*T[1*3+1])*w[1]
           +((r[1]-rc[1])*T[2*3+2]-(r[2]-rc[2])*T[1*3+2])*w[2])*ad->bd.wt_44[i];
    N[1]+=( ((r[2]-rc[2])*T[0*3+0]-(r[0]-rc[0])*T[2*3+0])*w[0]
           +((r[2]-rc[2])*T[0*3+1]-(r[0]-rc[0])*T[2*3+1])*w[1]
           +((r[2]-rc[2])*T[0*3+2]-(r[0]-rc[0])*T[2*3+2])*w[2])*ad->bd.wt_44[i];
    N[2]+=( ((r[0]-rc[0])*T[1*3+0]-(r[1]-rc[1])*T[0*3+0])*w[0]
           +((r[0]-rc[0])*T[1*3+1]-(r[1]-rc[1])*T[0*3+1])*w[1]
           +((r[0]-rc[0])*T[1*3+2]-(r[1]-rc[1])*T[0*3+2])*w[2])*ad->bd.wt_44[i];
  }

}

void bil_force_9p(double *F,double *N,double *rc,int s,DMDA *ad)
{
   double complex p,v[3];
  double r[3],w[3],T[9],c1,c2,c3,t1,cr[3][4],cw[3][3];
  int sd,asd,i,j,l,m;

  c1=0.25/(ad->aw.c0*ad->aw.c0*ad->aw.rho0); // 1/(4*K0)
  c2=0.25*ad->aw.rho0; // rho0/4
  c3=0.50*ad->aw.rho0; // rho0/2

  sd=ad->bd.sb[0].sid[s];
  asd=abs(sd);
  bil_copy_elem_const_rw(cr,cw,asd,&(ad->bd));
  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<9;i++){
    bil_rw_zeta_eta(r,w,ad->bd.zt_49[i],ad->bd.et_49[i],cr,cw);
    pv_t_bd_dmda(&p,v,0,s,ad->bd.zt_49[i],ad->bd.et_49[i],1,ad);
    
    // radiation stress tensor
    t1=c1*creal(p*conj(p))-c2*creal(v[0]*conj(v[0])+v[1]*conj(v[1])+v[2]*conj(v[2]));
    for(l=0;l<3;l++) for(m=0;m<3;m++) T[l*3+m]=c3*creal(v[l]*conj(v[m]));
    T[0*3+0]+=t1;
    T[1*3+1]+=t1;
    T[2*3+2]+=t1; 
    // force
    F[0]+=(T[0*3+0]*w[0]+T[0*3+1]*w[1]+T[0*3+2]*w[2])*ad->bd.wt_49[i];
    F[1]+=(T[1*3+0]*w[0]+T[1*3+1]*w[1]+T[1*3+2]*w[2])*ad->bd.wt_49[i];
    F[2]+=(T[2*3+0]*w[0]+T[2*3+1]*w[1]+T[2*3+2]*w[2])*ad->bd.wt_49[i]; 
    // torque
    N[0]+=( ((r[1]-rc[1])*T[2*3+0]-(r[2]-rc[2])*T[1*3+0])*w[0]
           +((r[1]-rc[1])*T[2*3+1]-(r[2]-rc[2])*T[1*3+1])*w[1]
           +((r[1]-rc[1])*T[2*3+2]-(r[2]-rc[2])*T[1*3+2])*w[2])*ad->bd.wt_49[i];
    N[1]+=( ((r[2]-rc[2])*T[0*3+0]-(r[0]-rc[0])*T[2*3+0])*w[0]
           +((r[2]-rc[2])*T[0*3+1]-(r[0]-rc[0])*T[2*3+1])*w[1]
           +((r[2]-rc[2])*T[0*3+2]-(r[0]-rc[0])*T[2*3+2])*w[2])*ad->bd.wt_49[i];
    N[2]+=( ((r[0]-rc[0])*T[1*3+0]-(r[1]-rc[1])*T[0*3+0])*w[0]
           +((r[0]-rc[0])*T[1*3+1]-(r[1]-rc[1])*T[0*3+1])*w[1]
           +((r[0]-rc[0])*T[1*3+2]-(r[1]-rc[1])*T[0*3+2])*w[2])*ad->bd.wt_49[i];
  }
}

void lit_force_4p(double *F,double *N,double *rc,int s,DMDA *ad)
{
  double complex p,v[3];
  double r[3],w[3],cr[3][4],cw[3][3],t1,T[9],c1,c2,c3;
  int i,j,sd,asd,l,m;
  
  c1=0.25/(ad->aw.c0*ad->aw.c0*ad->aw.rho0); // 1/(4*K0)
  c2=0.25*ad->aw.rho0; // rho0/4
  c3=0.50*ad->aw.rho0; // rho0/2
  
  sd=ad->bd.sb[0].sid[s];
  asd=abs(sd);
  lit_copy_elem_const_rw(cr,cw,asd,&(ad->bd));

  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<4;i++){
    if(i<3){
      for(j=0;j<3;j++){
        r[j]=ad->bd.ren[asd][i][j];
        w[j]=ad->bd.wen[asd][i][j];
      }
      calc_maw_pv(&p,v,r,&(ad->aw)); // multi_aw.h
      p+=-ad->bd.sb[0].P[s][i]*ad->k2[0];
      for(j=0;j<3;j++) v[j]+=ad->bd.sb[0].pv[s][i][j];
    }
    else {
      lit_rw_zeta_eta(r,w,ad->bd.zt_34[i],ad->bd.et_34[i],cr,cw);
      pv_t_bd_dmda(&p,v,0,s,ad->bd.zt_34[i],ad->bd.et_34[i],0,ad);
    }
    
    // radiation stress tensor
    t1=c1*creal(p*conj(p))-c2*creal(v[0]*conj(v[0])+v[1]*conj(v[1])+v[2]*conj(v[2]));
    for(l=0;l<3;l++) for(m=0;m<3;m++) T[l*3+m]=c3*creal(v[l]*conj(v[m]));
    T[0*3+0]+=t1;
    T[1*3+1]+=t1;
    T[2*3+2]+=t1; 
    // force
    F[0]+=(T[0*3+0]*w[0]+T[0*3+1]*w[1]+T[0*3+2]*w[2])*ad->bd.wt_34[i];
    F[1]+=(T[1*3+0]*w[0]+T[1*3+1]*w[1]+T[1*3+2]*w[2])*ad->bd.wt_34[i];
    F[2]+=(T[2*3+0]*w[0]+T[2*3+1]*w[1]+T[2*3+2]*w[2])*ad->bd.wt_34[i]; 
    // torque
    N[0]+=( ((r[1]-rc[1])*T[2*3+0]-(r[2]-rc[2])*T[1*3+0])*w[0]
           +((r[1]-rc[1])*T[2*3+1]-(r[2]-rc[2])*T[1*3+1])*w[1]
           +((r[1]-rc[1])*T[2*3+2]-(r[2]-rc[2])*T[1*3+2])*w[2])*ad->bd.wt_34[i];
    N[1]+=( ((r[2]-rc[2])*T[0*3+0]-(r[0]-rc[0])*T[2*3+0])*w[0]
           +((r[2]-rc[2])*T[0*3+1]-(r[0]-rc[0])*T[2*3+1])*w[1]
           +((r[2]-rc[2])*T[0*3+2]-(r[0]-rc[0])*T[2*3+2])*w[2])*ad->bd.wt_34[i];
    N[2]+=( ((r[0]-rc[0])*T[1*3+0]-(r[1]-rc[1])*T[0*3+0])*w[0]
           +((r[0]-rc[0])*T[1*3+1]-(r[1]-rc[1])*T[0*3+1])*w[1]
           +((r[0]-rc[0])*T[1*3+2]-(r[1]-rc[1])*T[0*3+2])*w[2])*ad->bd.wt_34[i];
  }

  for(j=0;j<3;j++){
    F[j]*=0.5;
    N[j]*=0.5;
  }
}

void lit_force_7p(double *F,double *N,double *rc,int s,DMDA *ad)
{
  double complex p,v[3];
  double r[3],w[3],cr[3][4],cw[3][3],c1,c2,c3,t1,T[9];
  int i,j,asd,sd,l,m;
  
  c1=0.25/(ad->aw.c0*ad->aw.c0*ad->aw.rho0); // 1/(4*K0)
  c2=0.25*ad->aw.rho0; // rho0/4
  c3=0.50*ad->aw.rho0; // rho0/2
  
  sd=ad->bd.sb[0].sid[s];
  asd=abs(sd);
  lit_copy_elem_const_rw(cr,cw,asd,&(ad->bd));

  for(j=0;j<3;j++){
    F[j]=0.0;
    N[j]=0.0;
  }

  for(i=0;i<7;i++){
    lit_rw_zeta_eta(r,w,ad->bd.zt_37[i],ad->bd.et_37[i],cr,cw);
    pv_t_bd_dmda(&p,v,0,s,ad->bd.zt_37[i],ad->bd.et_37[i],1,ad);
    
    // radiation stress tensor
    t1=c1*creal(p*conj(p))-c2*creal(v[0]*conj(v[0])+v[1]*conj(v[1])+v[2]*conj(v[2]));
    for(l=0;l<3;l++) for(m=0;m<3;m++) T[l*3+m]=c3*creal(v[l]*conj(v[m]));
    T[0*3+0]+=t1;
    T[1*3+1]+=t1;
    T[2*3+2]+=t1; 
    // force
    F[0]+=(T[0*3+0]*w[0]+T[0*3+1]*w[1]+T[0*3+2]*w[2])*ad->bd.wt_37[i];
    F[1]+=(T[1*3+0]*w[0]+T[1*3+1]*w[1]+T[1*3+2]*w[2])*ad->bd.wt_37[i];
    F[2]+=(T[2*3+0]*w[0]+T[2*3+1]*w[1]+T[2*3+2]*w[2])*ad->bd.wt_37[i]; 
    // torque
    N[0]+=( ((r[1]-rc[1])*T[2*3+0]-(r[2]-rc[2])*T[1*3+0])*w[0]
           +((r[1]-rc[1])*T[2*3+1]-(r[2]-rc[2])*T[1*3+1])*w[1]
           +((r[1]-rc[1])*T[2*3+2]-(r[2]-rc[2])*T[1*3+2])*w[2])*ad->bd.wt_37[i];
    N[1]+=( ((r[2]-rc[2])*T[0*3+0]-(r[0]-rc[0])*T[2*3+0])*w[0]
           +((r[2]-rc[2])*T[0*3+1]-(r[0]-rc[0])*T[2*3+1])*w[1]
           +((r[2]-rc[2])*T[0*3+2]-(r[0]-rc[0])*T[2*3+2])*w[2])*ad->bd.wt_37[i];
    N[2]+=( ((r[0]-rc[0])*T[1*3+0]-(r[1]-rc[1])*T[0*3+0])*w[0]
           +((r[0]-rc[0])*T[1*3+1]-(r[1]-rc[1])*T[0*3+1])*w[1]
           +((r[0]-rc[0])*T[1*3+2]-(r[1]-rc[1])*T[0*3+2])*w[2])*ad->bd.wt_37[i];
  }

  for(j=0;j<3;j++){
    F[j]*=0.5;
    N[j]*=0.5;
  }  
}

void bil_sa_4p(double *S,int d,int s,DMDA *ad)
{
  double w[3];
  int sd,asd,i,j;

  sd=ad->bd.sb[d].sid[s];
  asd=abs(sd);
  *S=0.0;
  for(i=0;i<4;i++){
    for(j=0;j<3;j++) w[j]=ad->bd.wen[asd][i][j];
    *S+=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*ad->bd.wt_44[i];
  }
}

void bil_sa_9p(double *S,int d,int s,DMDA *ad)
{
  double r[3],w[3],cr[3][4],cw[3][3];
  int sd,asd,i;

  sd=ad->bd.sb[d].sid[s];
  asd=abs(sd);
  bil_copy_elem_const_rw(cr,cw,asd,&(ad->bd));
  *S=0.0;

  for(i=0;i<9;i++){
    bil_rw_zeta_eta(r,w,ad->bd.zt_49[i],ad->bd.et_49[i],cr,cw);
    *S+=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*ad->bd.wt_49[i];
  }
}

void lit_sa_4p(double *S,int d,int s,DMDA *ad)
{
  double r[3],w[3],cr[3][4],cw[3][3];
  int i,j,sd,asd;

  sd=ad->bd.sb[d].sid[s];
  asd=abs(sd);
  lit_copy_elem_const_rw(cr,cw,asd,&(ad->bd));

  *S=0.0;
  for(i=0;i<4;i++){
    if(i<3){
      for(j=0;j<3;j++){
        r[j]=ad->bd.ren[asd][i][j];
        w[j]=ad->bd.wen[asd][i][j];
      }
    }
    else {
      lit_rw_zeta_eta(r,w,ad->bd.zt_34[i],ad->bd.et_34[i],cr,cw);
    }
    *S+=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*ad->bd.wt_34[i];
  }
  *S*=0.5;
}

void lit_sa_7p(double *S,int d,int s,DMDA *ad)
{
  double r[3],w[3],cr[3][4],cw[3][3];
  int i,asd,sd;

  sd=ad->bd.sb[d].sid[s];
  asd=abs(sd);
  lit_copy_elem_const_rw(cr,cw,asd,&(ad->bd));

  *S=0.0;
  for(i=0;i<7;i++){
    lit_rw_zeta_eta(r,w,ad->bd.zt_37[i],ad->bd.et_37[i],cr,cw);
    *S+=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2])*ad->bd.wt_37[i];
  }
  *S*=0.5;
}
