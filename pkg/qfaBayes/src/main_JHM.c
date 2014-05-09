#include "headers_JHM.h"
#include "datain_JHM.h"
#include "functions_JHM.h"
#include "print_JHM.h"

int main_JHM(int *arga,int *argb,int *argc,int *argd,double *OUT,
  char **HEADER,int *QFAIA,double *QFADyA,double *QFADxA,int *QFADNoORFA,
  int *QFANoTIMEA,int *QFAIB,double *QFADyB,double *QFADxB,
  int *QFADNoORFB,int *QFANoTIMEB,double *PRIORS,double *TUNING)
{
  GetRNGstate();

  struct_data_JHM *data= malloc(sizeof(struct_data_JHM));
  struct_para_JHM *para= malloc(sizeof(struct_para_JHM));
  struct_priors_JHM *priors= malloc(sizeof(struct_priors_JHM));
  struct_tuning_JHM *tuning = malloc(sizeof(struct_tuning_JHM));
  struct_adaptive_JHM *adaptive = malloc(sizeof(struct_tuning_JHM));

  int burn,iters,thin,adaptive_phase;
  
  burn=*arga;   /*burn in*/
  iters=*argb;    /*iterations*/
  thin=*argc;        /*thinning*/
  adaptive_phase=*argd; /*adaptive phase length within burn in*/
  
  inzstruct_data_JHM(data,QFAIA,QFADyA,QFADxA,QFADNoORFA,QFANoTIMEA,QFAIB,QFADyB,QFADxB,QFADNoORFB,QFANoTIMEB);
  inzstruct_tuning_JHM(tuning,TUNING);
  inzstruct_adaptive_JHM(adaptive);
  inzstruct_priors_JHM(priors,PRIORS);
  inzstruct_para_JHM(para,data,priors);

  gibbsandMHloop_JHM(burn,1,data,para,priors,tuning,adaptive,0,adaptive_phase,OUT,HEADER);
  gibbsandMHloop_JHM(iters,thin,data,para,priors,tuning,adaptive,1,0,OUT,HEADER);

  PutRNGstate();
  
  return 0;
}
