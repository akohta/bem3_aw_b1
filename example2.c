// analysis sample of intensity distribution.
#include "bem3_aw_b1.h"

int main(int argc,char *argv[]) 
{
  DMDA ad;
  FILE *fp1;
  double complex p,v[3];
  double rang,dr,r[3],*ip,mf;
  int max,i,j,type,sn;

  if(argc!=2 && argc!=5){
    printf("Usage : %s datafile_name [sampling_number multplier_factor type](optional)\n",argv[0]);
    printf("default sampling number 200, multiplier factor 3 (range is -3*lambda0 to 3*lambda0), type 1 (9 or 7 point GaussLegendre)\n");
    exit(0);
  }
  else if(argc==5){
    sn=atoi(argv[2]);
    mf=atof(argv[3]);
    type=atoi(argv[4]);
  }
  else{
    sn=200;
    mf=3.0;
    type=1;
  }

  dat_read_dmda(argv[1],&ad); // read datafile 
  print_dmda(&ad);            // print data 

  max=sn;
  rang=mf*ad.aw.lambda0;
  dr=rang*2.0/(double)(max-1);
  
  ip=(double *)m_alloc2(max,sizeof(double),"exampl2.c,ip");
  
  // x=0 plane
  if((fp1=fopen("Ip_yz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# y z sound_pressure_intensity");
  r[0]=0.0;
  for(i=0;i<max;i++){
    r[1]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(p,v) // omp parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      pv_t_dmda(&p,v,r,type,&ad); // total field
      ip[j]=creal(p*conj(p));     // |p|^2
    } // end parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[1],r[2],ip[j]);
    }
    fprintf(fp1,"\n");
  }
  fclose(fp1);

  // y=0 plane
  if((fp1=fopen("Ip_xz.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  }
  fprintf(fp1,"%s\n","# x z sound_pressure_intensity");
  r[1]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(p,v) // omp parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      pv_t_dmda(&p,v,r,type,&ad); // total field
      ip[j]=creal(p*conj(p));
    }// end parallel
    for(j=0;j<max;j++){
      r[2]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[0],r[2],ip[j]);
    }
    fprintf(fp1,"\n");
  }
  fclose(fp1);

  // z=0 plane  
  if((fp1=fopen("Ip_xy.txt","wt"))==NULL){    printf("Can not open the file.\n");    exit(1);  } 
  fprintf(fp1,"%s\n","# x y sound_pressure_intensity");
  r[2]=0.0;
  for(i=0;i<max;i++){
    r[0]=-rang+(double)i*dr;
    #pragma omp parallel for schedule(dynamic) firstprivate(r) private(p,v) // omp parallel
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      pv_t_dmda(&p,v,r,type,&ad); // total field
      ip[j]=creal(p*conj(p));
    }
    for(j=0;j<max;j++){
      r[1]=-rang+(double)j*dr;
      fprintf(fp1,"%g %g %15.14e\n",r[0],r[1],ip[j]);
    }
    fprintf(fp1,"\n");
  }
  fclose(fp1);
  
  printf("Intensity plot is finished\n");
  
  free(ip);
  finalize_dmda(&ad);

  return 0;
}
