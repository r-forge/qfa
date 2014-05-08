#include "headers_SHM.h"
#include "datain_SHM.h"
#include "functions_SHM.h"
#include "print_SHM.h"

int main_SHM(int *arga,int *argb,int *argc,int *argd,int *arge,double *OUT,
  char **HEADER,int *QFAI,double *QFADy,double *QFADx,int *QFADNoORF,
  int *QFADNoTIME,double *PRIORS,double *TUNING)
{
  GetRNGstate();
  
  struct_data *data= malloc(sizeof(struct_data));
  struct_para *para= malloc(sizeof(struct_para));
  struct_priors *priors= malloc(sizeof(struct_priors));
  struct_tuning *tuning = malloc(sizeof(struct_tuning));
  struct_adaptive *adaptive = malloc(sizeof(struct_tuning));
  
  int burn,iters,thin,adaptive_phase,adaptive_period;
  
  burn=*arga;    /*burn in*/
  iters=*argb;    /*iterations*/
  thin=*argc;         /*thinning*/
  adaptive_phase=*argd; /*adaptive phase length within burn in*/
  adaptive_period=*arge; /*period to calculate adaptive statistic over*/

  inzstruct_data(data,QFAI,QFADy,QFADx,QFADNoORF,QFADNoTIME);
  inzstruct_tuning(tuning,TUNING);
  inzstruct_adaptive(adaptive);
  inzstruct_priors(priors,PRIORS);
  inzstruct_para(para,data,priors);

  gibbsandMHloop(burn,1,data,para,priors,tuning,adaptive,0,adaptive_phase,adaptive_period,OUT,HEADER);
  gibbsandMHloop(iters,thin,data,para,priors,tuning,adaptive,1,0,0,OUT,HEADER);
  
  PutRNGstate();
  
  return 0;
}
