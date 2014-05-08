#include "headers_SHM.h"
#include "datain_SHM.h"

/*INZ*/

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

int inzstruct_tuning(struct_tuning *tuning,double *TUNING)
{
    tuning->K_lm=TUNING[0],        tuning->tau_K_l=TUNING[1],
    tuning->r_lm=TUNING[2],        tuning->tau_r_l=TUNING[3],

    tuning->K_o_l=TUNING[4],        tuning->sigma_K_o=TUNING[5],
    tuning->r_o_l=TUNING[6],        tuning->sigma_r_o=TUNING[7],
    tuning->nu_l=TUNING[8],         tuning->sigma_nu=TUNING[9],

	tuning->tau_K_p=TUNING[16],  	tuning->sigma_tau_K=TUNING[17],
	tuning->tau_r_p=TUNING[18],     tuning->sigma_tau_r=TUNING[19],
	
    tuning->K_p=TUNING[20],			tuning->r_p=TUNING[21],
    tuning->P=TUNING[22];

return 0;
}

int inzstruct_adaptive(struct_adaptive *adaptive)
{
    adaptive->K_lm=0,            adaptive->tau_K_l=0,
    adaptive->r_lm=0,            adaptive->tau_r_l=0,

    adaptive->K_o_l=0,            adaptive->sigma_K_o=0,
    adaptive->r_o_l=0,            adaptive->sigma_r_o=0,
    adaptive->nu_l=0,             adaptive->sigma_nu=0,
   
	adaptive->tau_K_p=0,		  adaptive->sigma_tau_K=0,
	adaptive->tau_r_p=0,		  adaptive->sigma_tau_r=0,
	
    adaptive->K_p=0,			  adaptive->r_p=0,
    adaptive->P=0;
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
  
  D_para->sigma_K_o=D_priors->eta_K_o;               /*Precision*/
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
