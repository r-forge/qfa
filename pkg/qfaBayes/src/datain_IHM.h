#ifndef _datain_h
#define _datain_h

int testargc_IHM(int argc);

int datadouble_IHM(struct_data_IHM *D,double *QFADyA,double *QFADyB);

int inzstruct_priors_IHM(struct_priors_IHM *D_priors,double *PRIORS);
int inzstruct_data_IHM(struct_data_IHM *data,int *QFAIA,double *QFADyA,int *QFADNoORFA,int *QFAIB,double *QFADyB,int *QFADNoORFB);
int inzstruct_para_IHM(struct_para_IHM *para,struct_data_IHM *data,struct_priors_IHM *priors);

int filldata_IHM(struct_data_IHM *D);
int fillpara_IHM(struct_para_IHM *D_para,struct_data_IHM *D,struct_priors_IHM *D_priors);
int fillpriors_IHM(struct_priors_IHM *D_priors);

int inzstruct_tuning_IHM(struct_tuning_IHM *tuning,double *TUNING);
int inzstruct_adaptive_IHM(struct_adaptive_IHM *adaptive);
#endif
