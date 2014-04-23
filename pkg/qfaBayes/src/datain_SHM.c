#include "headers_SHM.h"
#include "datain_SHM.h"

/*INZ*/

int inzstruct_MH(struct_MH *MH,double *TUNING)
{
  fillMH(MH,TUNING);
  return 0;
}

int inzstruct_priors(struct_priors *D_priors,double *PRIORS)
{
   /*K*/
    D_priors->tau_K_mu=PRIORS[0];
   
    D_priors->eta_tau_K_p=PRIORS[1];          

    D_priors->eta_K_o=PRIORS[2]; 

    D_priors->psi_K_o=PRIORS[3];            
    /*r*/

    D_priors->tau_r_mu=PRIORS[4];              

    D_priors->eta_tau_r_p=PRIORS[5];          

    D_priors->eta_r_o=PRIORS[6];              

    D_priors->psi_r_o=PRIORS[7];       
    /*nu*/

    D_priors->eta_nu=PRIORS[8];         

    D_priors->psi_nu=PRIORS[9];      
    /*K*//*r*//*nu*//*P*/

    D_priors->K_mu=PRIORS[10];     

    D_priors->eta_K_p=PRIORS[11];      /*Normal  LMean; Precisions */

    D_priors->r_mu=PRIORS[12];   

    D_priors->eta_r_p=PRIORS[13];      /*Normal  LMean; Precisions */

    D_priors->nu_mu=PRIORS[14];      

    D_priors->eta_nu_p=PRIORS[15];     /*Normal  LMean; Precisions */

    D_priors->P_mu=PRIORS[16];   

    D_priors->eta_P=PRIORS[17];   /*Normal  LMean; Precisions */
    D_priors->df=3;
    D_priors->eta_tau_K=D_priors->eta_tau_K_p;  D_priors->psi_tau_K=D_priors ->eta_tau_K_p;
    D_priors->eta_tau_r=D_priors->eta_tau_r_p;  D_priors->psi_tau_r=D_priors->eta_tau_r_p;
  /*fillpriors(priors);*/
  return 0;
}

int inzstruct_data(struct_data *data,int *QFAI,double *QFADy,double *QFADx,int *QFADNoORF,int *QFADNoTIME)
{
	int i;
  long size;
	data->L=QFAI[0];
	data->M=QFAI[1];
	data->N=QFAI[2];
	data->maxy=QFAI[3];
	data->maxNoTIME=QFAI[4];
	
  size=data->L*data->M*data->N; /*input from file*/ 
  data->y=malloc(size*sizeof(double));        /*Cycle with SHIFTlmn*/
  data->x=malloc(size*sizeof(double));        /*Cycle with SHIFTlmn*/
  size=data->L;
  data->NoORF=malloc(size*sizeof(double));    /*Cycle with data->L*/
  data->NoSUM=malloc(size*sizeof(double));    /*Cycle with data->L*/

 for (i=0;i<(data->L*data->M*data->N);i++){
data->y[i]=QFADy[i];
data->x[i]=QFADx[i];
}
	for (i=0;i<(data->L*data->M*data->N);i++){
          if(data->y[i]<0){ data->y[i]=0;}
	  if(data->x[i]<0){ data->x[i]=0;}
        }


 for (i=0;i<(data->L);i++){
data->NoORF[i]=QFADNoORF[i];
}
 
  filldata(data);
  size=data->SHIFTlmn;/*inputfromfile*/
  data->NoTIME=malloc(size*sizeof(double));   /*Cycle with SHIFTlm*/

 for (i=0;i<(data->SHIFTlmn);i++){
data->NoTIME[i]=QFADNoTIME[i];
}

  return 0;
}

int inzstruct_para(struct_para *para,struct_data *data,struct_priors *priors)
{
  long size;

  size=data->SHIFTlmn;
  para->K_lm=malloc(size*sizeof(double));
  para->r_lm=malloc(size*sizeof(double));

  size=data->L;
  para->tau_K_l=malloc(size*sizeof(double));
  para->tau_r_l=malloc(size*sizeof(double));
  para->K_o_l=malloc(size*sizeof(double));
  para->r_o_l=malloc(size*sizeof(double));
  para->nu_l=malloc(size*sizeof(double));

  fillpara(para,data,priors);
  return 0;
}

/*FILL*/

int fillMH(struct_MH *MH,double *TUNING)
{
  MH->hK=TUNING[0];
  MH->hr=TUNING[1];
  MH->hnu=TUNING[2];
  MH->hP=TUNING[3];
  MH->accept_K=0;
  MH->accept_r=0;
  MH->accept_nu=0;
  MH->accept_P=0; 

/*HARDCODED VER
  MH->hK=0.1;MH->accept_K=0;
  MH->hr=0.1;MH->accept_r=0;
  MH->hnu=0.1;MH->accept_nu=0;
  MH->hP=0.2;MH->accept_P=0;  */
  return 0;
}

int filldata(struct_data *D)
{
  int l;

  D->NoSUM[0]=0;
  for (l=1;l<(D->L);l++){
    D->NoSUM[l]=D->NoSUM[l-1]+D->NoORF[l-1];
  }
  D->SHIFTlmn=D->NoSUM[D->L-1]+D->NoORF[D->L-1];/*create mnSHIFT*/
  return 0;
}

int fillpara(struct_para *D_para, struct_data *D,struct_priors *D_priors)
{
  int l,m,mm;
  double SUM=0,SUMa=0;
  /*initials*/
  /*K*/

 for (l=0;l<D->L;l++){
    for (m=0;m<D->NoORF[l];m++){
      mm=D->NoSUM[l]+m;
      if(D->y[l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]<=0){D_para->K_lm[mm]=D_priors->P_mu;SUM+=D_para->K_lm[mm];}
	else{     
	  D_para->K_lm[mm]=log(D->y[l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]);SUM+=D_para->K_lm[mm];
	}
    }
    D_para->K_o_l[l]=SUM/D->NoORF[l];
    SUM=0;
    SUMa+=D_para->K_o_l[l];
  }
 D_para->K_p=SUMa/D->L;       /*LMean*/

  for (l=0;l<D->L;l++)          {D_para->tau_K_l[l]=D_priors->tau_K_mu;}                  /*Precision*/
  
  D_para->sigma_K_o=D_para->sigma_K_o_b=D_priors->eta_K_o;               /*Precision*/
  /*r*/
  for (l=0;l<D->L;l++){
    for (m=0;m<D->NoORF[l];m++){
      mm=D->NoSUM[l]+m;
      D_para->r_lm[mm]=D_priors->r_mu;
    }
  }                          /*LMean*/

  for (l=0;l<D->L;l++)          {D_para->tau_r_l[l]=D_priors->tau_r_mu;}                  /*Precision*/

  for (l=0;l<D->L;l++)          {D_para->r_o_l[l]=D_priors->r_mu;}        /*LMean*/
  D_para->sigma_r_o=D_priors->eta_r_o;               /*Precision*/

  D_para->r_p=D_priors->r_mu;       /*LMean*/
  /*nu*/
  for (l=0;l<D->L;l++)          {D_para->nu_l[l]=D_priors->nu_mu;}                      /*LMean*/
  D_para->sigma_nu=D_priors->eta_nu;   /*Precision for lMean*/

  D_para->nu_p=D_priors->nu_mu;   /*LMean*/
  /*P*/
  D_para->P=D_priors->P_mu;      /*LMean*/


  D_para->tau_K_p=D_priors->tau_K_mu;
  D_para->sigma_tau_K=D_priors->eta_tau_K;
  D_para->tau_r_p=D_priors->tau_r_mu;
  D_para->sigma_tau_r=D_priors->eta_tau_r;

  return 0;
}

/* eof  */