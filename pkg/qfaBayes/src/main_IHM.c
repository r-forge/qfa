#include "headers_IHM.h"
#include "datain_IHM.h"
#include "functions_IHM.h"
#include "print_IHM.h"

int main_IHM(int *arga,int *argb,int *argc,int *argd,int *arge,double *OUT,
  char **HEADER,int *QFAIA,double *QFADyA,int *QFADNoORFA,int *QFAIB,
  double *QFADyB,int *QFADNoORFB,double *PRIORS,double *TUNING)
{
  GetRNGstate();

  struct_data_IHM *data= malloc(sizeof(struct_data_IHM));
  struct_para_IHM *para= malloc(sizeof(struct_para_IHM));
  struct_priors_IHM *priors= malloc(sizeof(struct_priors_IHM));
  struct_tuning_IHM *tuning = malloc(sizeof(struct_tuning_IHM));
  struct_adaptive_IHM *adaptive = malloc(sizeof(struct_tuning_IHM));
  
  int burn,iters,thin,adaptive_phase,adaptive_period;
  
  burn=*arga;   /*burn in*/
  iters=*argb;    /*iterations*/
  thin=*argc;        /*thinning*/
  adaptive_phase=*argd; /*adaptive phase length within burn in*/
  adaptive_period=*arge; /*period to calculate adaptive statistic over*/
  
  inzstruct_data_IHM(data,QFAIA,QFADyA,QFADNoORFA,QFAIB,QFADyB,QFADNoORFB);
  inzstruct_tuning_JHM(tuning,TUNING);
  inzstruct_adaptive_JHM(adaptive);
  inzstruct_priors_IHM(priors,PRIORS);
  inzstruct_para_IHM(para,data,priors);

  gibbsandMHloop_IHM(burn,1,data,para,priors,tuning,adaptive,0,adaptive_phase,adaptive_period,OUT,HEADER);
  gibbsandMHloop_IHM(iters,thin,data,para,priors,tuning,adaptive,1,0,0,OUT,HEADER);
  
  PutRNGstate();
  
  return 0;
}
