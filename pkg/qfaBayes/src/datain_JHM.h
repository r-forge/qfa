#ifndef _datain_h
#define _datain_h

int testargc_JHM(int argc);
int testsame_JHM(int a,int b);

int datadouble_JHM(char filename[], char filename2[], double datavec[],int length);
int dataint_JHM(char filename[],char filename2[], int datavec[] ,int lengtha,int lengthb );
int dataLMN_JHM(char filename[],char filename2[], 
int *datavecL,int *datavecM,int *datavecN,int *datavecmaxy,int *datavecmaxTIMEa,int *datavecmaxTIMEb);

int inzstruct_priors_JHM(struct_priors_JHM *D_priors,double *PRIORS);
int inzstruct_data_JHM(struct_data_JHM *data,int *QFAIA,double *QFADyA,double *QFADxA,int *QFADNoORFA,int *QFADNoTIMEA,int *QFAIB,double *QFADyB,double *QFADxB,int *QFADNoORFB,int *QFADNoTIMEB);
int inzstruct_para_JHM(struct_para_JHM *para,struct_data_JHM *data,struct_priors_JHM *priors);

int filldata_JHM(struct_data_JHM *D);
int fillpara_JHM(struct_para_JHM *D_para,struct_data_JHM *D,struct_priors_JHM *D_priors);
int fillpriors_JHM(struct_priors_JHM *D_priors);

int inzstruct_tuning_JHM(struct_tuning_JHM *tuning,double *TUNING);
int inzstruct_adaptive_JHM(struct_adaptive_JHM *adaptive);
#endif
