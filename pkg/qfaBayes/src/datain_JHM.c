#include "headers_JHM.h"
#include "datain_JHM.h"

/*INZ*/

int inzstruct_priors_JHM(struct_priors_JHM *D_priors,double *PRIORS)
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

   D_priors->eta_K_p=PRIORS[11];     

 D_priors->r_mu=PRIORS[12];       

     D_priors->eta_r_p=PRIORS[13];      

	D_priors->nu_mu=PRIORS[14];           

 D_priors->eta_nu_p=PRIORS[15];    

	D_priors->P_mu=PRIORS[16];     

   D_priors->eta_P=PRIORS[17];   

	/*data2.c*/       

	D_priors->alpha_mu=PRIORS[18];     

     D_priors->eta_alpha=PRIORS[19];

	D_priors->beta_mu=PRIORS[20];        

   D_priors->eta_beta=PRIORS[21];

	D_priors->p=PRIORS[22];    

	D_priors->eta_gamma=PRIORS[23];

	D_priors->psi_gamma=PRIORS[24];

	D_priors->eta_omega=PRIORS[25];	

D_priors->psi_omega=PRIORS[26];
   
   	D_priors->eta_tau_K=PRIORS[27];  
	D_priors->psi_tau_K=PRIORS[28];
	D_priors->eta_tau_r=PRIORS[29];  
	D_priors->psi_tau_r=PRIORS[30];
	
	D_priors->df=3;
	/*fillpriors(priors);*/

return 0;
}

int inzstruct_data_JHM(struct_data_JHM *data,int *QFAIA,double *QFADyA,double *QFADxA,int *QFADNoORFA,int *QFADNoTIMEA,int *QFAIB,double *QFADyB,double *QFADxB,int *QFADNoORFB,int *QFADNoTIMEB)
{
	int i;
	long size;

		data->L=QFAIA[0];     

		data->M=QFAIA[1];
		data->N=QFAIA[2];
		data->maxTIMEa=QFAIA[4];
		data->L=QFAIB[0];
		data->M=QFAIB[1];
		data->N=QFAIB[2];
		data->maxy=QFAIA[3]+QFAIB[3];
		data->maxTIMEb=QFAIB[4];

	data->SHIFTlmn=QFAIA[3];
	size=data->maxy;
  	data->y=malloc(size*sizeof(double));  
        data->x=malloc(size*sizeof(double));  
	size=data->L*2;
	data->NoORF=malloc(size*sizeof(double));  
	data->NoSUM=malloc(size*sizeof(double));
	size=data->maxTIMEa+data->maxTIMEb;
	data->NoTIME=malloc(size*sizeof(double));  
/**/
 for (i=0;i<(data->maxy-data->SHIFTlmn);i++){
data->y[i]=QFADyA[i];
data->x[i]=QFADxA[i];
data->y[i+data->SHIFTlmn]=QFADyB[i];
data->x[i+data->SHIFTlmn]=QFADxB[i];
}

	for (i=0;i<data->maxy;i++){
          if(data->y[i]<0){ data->y[i]=0;}
	  if(data->x[i]<0){ data->x[i]=0;}
        }


 for (i=0;i<(data->L);i++){
data->NoORF[i]=QFADNoORFA[i];
data->NoORF[i+data->L]=QFADNoORFB[i];
}
 for (i=0;i<(data->maxTIMEa);i++){
data->NoTIME[i]=QFADNoTIMEA[i];
}
 for (i=0;i<(data->maxTIMEb);i++){
data->NoTIME[i+data->maxTIMEa]=QFADNoTIMEB[i];
}
/**/
	filldata_JHM(data);
return 0;
}

int inzstruct_para_JHM(struct_para_JHM *para,struct_data_JHM *data,struct_priors_JHM *priors)
{
	long size;
	size=data->L*2;
	para->nu_l=malloc(size*sizeof(double));
	para->tau_K_cl=malloc(size*sizeof(double));
	para->tau_r_cl=malloc(size*sizeof(double));
	size=data->maxTIMEa+data->maxTIMEb;/*inputfromfile*/
	para->K_clm=malloc(size*sizeof(double));
	para->r_clm=malloc(size*sizeof(double));
	size=data->L;
	para->delta_l=malloc(size*sizeof(double));
	para->gamma_cl=malloc(size*sizeof(double));
	para->omega_cl=malloc(size*sizeof(double));
	para->K_o_l=malloc(size*sizeof(double));
	para->r_o_l=malloc(size*sizeof(double));

	size=2;
	para->alpha_c=malloc(size*sizeof(double));
	para->beta_c=malloc(size*sizeof(double));

	para->tau_K_p=malloc(size*sizeof(double));
	para->tau_r_p=malloc(size*sizeof(double));
	para->sigma_tau_K=malloc(size*sizeof(double));
	para->sigma_tau_r=malloc(size*sizeof(double));

	fillpara_JHM(para,data,priors);
return 0;
}

/*FILL*/

int  inzstruct_tuning_JHM(struct_tuning_JHM *tuning,double *TUNING)
{
    tuning->K_clm=TUNING[0],        tuning->tau_K_cl=TUNING[1],
    tuning->r_clm=TUNING[2],        tuning->tau_r_cl=TUNING[3],

    tuning->K_o_l=TUNING[4],        tuning->sigma_K_o=TUNING[5],
    tuning->r_o_l=TUNING[6],        tuning->sigma_r_o=TUNING[7],
    tuning->nu_l=TUNING[8],         tuning->sigma_nu=TUNING[9],

    tuning->gamma_cl=TUNING[10],	tuning->sigma_gamma=TUNING[11],
    tuning->omega_cl=TUNING[12],    tuning->sigma_omega=TUNING[13],
	
    tuning->alpha_c=TUNING[14],		tuning->beta_c=TUNING[15],
	
	tuning->tau_K_p=TUNING[16],  	tuning->sigma_tau_K=TUNING[17],
	tuning->tau_r_p=TUNING[18],     tuning->sigma_tau_r=TUNING[19],
	
    tuning->K_p=TUNING[20],			tuning->r_p=TUNING[21],
    tuning->P=TUNING[22];

return 0;
}

int  inzstruct_adaptive_JHM(struct_adaptive_JHM *adaptive)
{
    adaptive->K_clm=0,            adaptive->tau_K_cl=0,
    adaptive->r_clm=0,            adaptive->tau_r_cl=0,

    adaptive->K_o_l=0,            adaptive->sigma_K_o=0,
    adaptive->r_o_l=0,            adaptive->sigma_r_o=0,
    adaptive->nu_l=0,             adaptive->sigma_nu=0,
   
    adaptive->gamma_cl=0,	      adaptive->sigma_gamma=0,
    adaptive->omega_cl=0,         adaptive->sigma_omega=0,
    
	adaptive->alpha_c=0,		  adaptive->beta_c=0,
	
	adaptive->tau_K_p=0,		  adaptive->sigma_tau_K=0,
	adaptive->tau_r_p=0,		  adaptive->sigma_tau_r=0,
	
    adaptive->K_p=0,			  adaptive->r_p=0,
    adaptive->P=0;
return 0;
}

int filldata_JHM(struct_data_JHM *D)
{
	int l;
	D->NoSUM[0]=0;
	for (l=1;l<(D->L);l++){
		D->NoSUM[l]=D->NoSUM[l-1]+D->NoORF[l-1];
	}

	D->NoSUM[D->L]=D->SHIFTmn=D->NoSUM[D->L-1]+D->NoORF[D->L-1];/*create mnSHIFT*/

	for (l=(1+D->L);l<(2*D->L);l++){
		D->NoSUM[l]=D->NoSUM[l-1]+D->NoORF[l-1];
	}


return 0;
}

int fillpara_JHM(struct_para_JHM *D_para, struct_data_JHM *D,struct_priors_JHM *D_priors)
{
int c,l,m,ll,mm;
 double SUM=0,SUMa=0,SUMb=0;
	/*initials*/
  /*K*/
for (c=0;c<2;c++){
    for (l=0;l<D->L;l++){
      ll=c*D->L+l;
      for (m=0;m<D->NoORF[ll];m++){
	mm=D->NoSUM[ll]+m;
	if (D->y[c*D->SHIFTlmn+l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]<=0){D_para->K_clm[mm]=D_priors->P_mu;}
	else{     
		D_para->K_clm[mm]=log(D->y[c*D->SHIFTlmn+l*D->M*D->N + m*D->N + D->NoTIME[mm]-1]);
		}
	if (c==0){SUMa+=exp(D_para->K_clm[mm]);SUM+=D_para->K_clm[mm];}
	  else{ SUMb+=exp(D_para->K_clm[mm]);}
      }
      if(c==0){D_para->K_o_l[l]=SUM/D->NoORF[l];}
      SUM=0;
      }
  }
 D_para->alpha_c[1]=log((SUMb/(D->NoORF[D->L-1]+D->NoSUM[2*D->L-1]-D->NoSUM[D->L]))/(SUMa/D->NoSUM[D->L]));
 D_para->beta_c[1]=D_para->alpha_c[1];
 SUM=0;
 for (l=0;l<(D->L);l++){SUM+=D_para->K_o_l[l];}
 D_para->K_p=SUM/D->L;		
	for (l=0;l<(2*D->L);l++)          {D_para->tau_K_cl[l]=D_priors->tau_K_mu;}                  /*Precision*/
       
	D_para->sigma_K_o=D_priors->eta_K_o;               /*Precision*/

	/*r*/
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
			ll=c*D->L+l;
			for (m=0;m<D->NoORF[ll];m++){
				mm=D->NoSUM[ll]+m;
				D_para->r_clm[mm]=D_priors->r_mu;   
			} 
		}  
	}
               
	for (l=0;l<2*D->L;l++)          {D_para->tau_r_cl[l]=D_priors->tau_r_mu;}                  /*Precision*/

	for (l=0;l<D->L;l++)          {D_para->r_o_l[l]=D_priors->r_mu;}     	   /*LMean*/
	D_para->sigma_r_o=D_priors->eta_r_o;               /*Precision*/

	D_para->r_p=D_priors->r_mu;       /*LMean*/
	
	/*nu*/
	for (l=0;l<2*D->L;l++)          {D_para->nu_l[l]=D_priors->nu_mu;}                      /*LMean*/
	D_para->sigma_nu=D_priors->eta_nu;   /*Precision for lMean*/

	D_para->nu_p=D_priors->nu_mu;   /*LMean*/
	/*P*/
  D_para->P=D_priors->P_mu;      /*LMean*/

	for (l=0;l<D->L;l++)          {D_para->gamma_cl[l]=0;} 

	for (l=0;l<D->L;l++)          {D_para->omega_cl[l]=0;}
	for (l=0;l<D->L;l++)          {D_para->delta_l[l]=1;}/*!*/  
 
	D_para->alpha_c[0]=0;
	D_para->beta_c[0]=0;
	D_para->sigma_gamma=D_priors->eta_gamma;
	D_para->sigma_omega=D_priors->eta_omega;
 
	for (c=0;c<2;c++){
	D_para->tau_K_p[c]=D_priors->tau_K_mu;
	D_para->sigma_tau_K[c]=D_priors->eta_tau_K;
	D_para->tau_r_p[c]=D_priors->tau_r_mu;
	D_para->sigma_tau_r[c]=D_priors->eta_tau_r;
	}
return 0;
}

