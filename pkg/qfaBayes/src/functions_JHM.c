#include "headers_JHM.h"
#include "functions_JHM.h"
#include "print_JHM.h"

/*trun*/
double trun_const_low_JHM(double x,double m, double l)
{
  return(1-0.5-(  (pow(3*l,0.5))*(x-m)/(((l*pow(x-m,2))+3)) +atan((pow(l/3,0.5))*(x-m)))/M_PI);
}
/*trun*/

/*Logistic Growth*/

double logistic_function_JHM(double t,double K, double r, double P){
	double output;
	output=(K*P*exp(r*t))/(K+P*(exp(r*t)-1));
	return(output);
}

double logistic_function_E_JHM(double t,double K, double r, double P){
	double output;
	output=logistic_function_JHM(t,exp(K), exp(r), exp(P));
	return(output);
}

/*Gibbs*/

double gauss_sample_JHM(struct_data_JHM *D ,int start, int N,double x[],double tau,double mu_0,double tau_0){
	double vec,Ndouble,SUM=0;
	int i;
	Ndouble=N;
	for (i=start;i<(start+N);i++){SUM=SUM+x[i];}
	vec=rnorm((tau_0*mu_0+tau*SUM)/(tau_0+Ndouble*tau),1/sqrt(tau_0+Ndouble*tau));
	return(vec);
}

/*MCMC MH*/

double MCMC_base_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data_JHM *D,struct struct_para_JHM *D_para,struct struct_priors_JHM *D_priors,double,int,int,int),int c,int l, int m){
	double logu,logaprob,can;
	can=rnorm(para,*h);
	logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

  logu=log(1-runif(0,1));
	if (logaprob>logu){para=can;*accept=*accept+1;}
	return(para); 
	}

double MCMC_base_truncate_low_JHM(double truncate,struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data_JHM *D,struct struct_para_JHM *D_para,struct struct_priors_JHM *D_priors,double,int,int,int),int c,int l, int m){
  double logu,logaprob,can;

	can=rnorm(para,*h);
    if(can<=(truncate)){
      can=para;
    }
 
  logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

  logu=log(1-runif(0,1));
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para);
}
double MCMC_base_truncate_high_JHM(double truncate,struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data_JHM *D,struct struct_para_JHM *D_para,struct struct_priors_JHM *D_priors,double,int,int,int),int c,int l, int m){
  double logu,logaprob,can;
 
	can=rnorm(para,*h);
    if(can>=(truncate)){
      can=para;
    }
  logaprob=(*foo)(D,D_para,D_priors,can,c,l,m)-(*foo)(D,D_para,D_priors,para,c,l,m);

  logu=log(1-runif(0,1));
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para);
}



double MCMC_K_clm_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll=c*D->L+l;
	mm=D->NoSUM[ll]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N+n;
		F=D->y[nn]-logistic_function_E_JHM(D->x[nn],para, D_para->r_clm[mm],D_para->P);
		SUM+=F*F*exp(D_para->nu_l[ll]);
	}	
	F=para-(D_para->alpha_c[c]+D_para->K_o_l[l]+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
	density=F*F*exp(D_para->tau_K_cl[ll])+SUM; 
	return(-0.5*density);
}

double MCMC_r_clm_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll=c*D->L+l;
	mm=D->NoSUM[ll]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N+n;
		F=D->y[nn]-logistic_function_E_JHM(D->x[nn], D_para->K_clm[mm],para,D_para->P);	
		SUM+=F*F*exp(D_para->nu_l[ll]);
	}	
	F=para-(D_para->beta_c[c]+D_para->r_o_l[l]+c*D_para->delta_l[l]*D_para->omega_cl[l]);
	density=F*F*exp(D_para->tau_r_cl[ll])+SUM; 
	return(-0.5*density); 
}


double MCMC_tau_K_cl_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,SUM=0,F;
	int mm,ll=c*D->L+l;
	for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->K_clm[mm]-(D_para->alpha_c[c]+D_para->K_o_l[l]+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
		SUM+=-para+F*F*exp(para)/**/+2*log(pnorm(0,(D_para->alpha_c[c]+D_para->K_o_l[l]+c*D_para->delta_l[l]*D_para->gamma_cl[l]),1/pow(exp(para),0.5),1,0));
	}	
	F=para-D_para->tau_K_p[c];
	density=F*F*exp(D_para->sigma_tau_K[c])+SUM; 
	return(-0.5*density); 
}

double MCMC_tau_r_cl_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,SUM=0,F;
	int mm,ll=c*D->L+l;
	for (m=0;m<D->NoORF[ll];m++){
		mm=D->NoSUM[ll]+m;
		F=D_para->r_clm[mm]-(D_para->beta_c[c]+D_para->r_o_l[l]+c*D_para->delta_l[l]*D_para->omega_cl[l]);
		SUM+=-para+F*F*exp(para)/**/+2*log(pnorm(3.5,(D_para->beta_c[c]+D_para->r_o_l[l]+c*D_para->delta_l[l]*D_para->omega_cl[l]),1/pow(exp(para),0.5),1,0));
	}	
	F=para-D_para->tau_r_p[c];
	density=F*F*exp(D_para->sigma_tau_r[c])+SUM; 
	return(-0.5*density); 
}



double MCMC_nu_l_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll;
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			for (n=0;n<D->NoTIME[mm];n++){
				nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N+n;
				F=D->y[nn]-logistic_function_E_JHM(D->x[nn],D_para->K_clm[mm], D_para->r_clm[mm],D_para->P);
				SUM+=-para+F*F*exp(para);
			}
		}
       
	F=para-D_para->nu_p;
	density=F*F*exp(D_para->sigma_nu)+SUM; 
	return(-0.5*density); 
}

double MCMC_sigma_nu_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	for (l=0;l<2*D->L;l++){
		F=D_para->nu_l[l]-D_para->nu_p;
		SUM+=-para+F*F*exp(para);
	}	
	F=para-D_priors->eta_nu;
	density=F*F*D_priors->psi_nu+SUM; 
	return(-0.5*density); 
}

double MCMC_alpha_c_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll,mm;
	for (l=0;l<D->L;l++){
		ll=D->L+l;
 		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->K_clm[mm]-(para+D_para->K_o_l[l]+D_para->delta_l[l]*D_para->gamma_cl[l]);
			SUM+=F*F*exp(D_para->tau_K_cl[ll])/**/+2*log(pnorm(0,(para+D_para->K_o_l[l]+D_para->delta_l[l]*D_para->gamma_cl[l]),1/pow(exp(D_para->tau_K_cl[ll]),0.5),1,0));
		}
	}	
	F=para-D_priors->alpha_mu;
	density=F*F*D_priors->eta_alpha+SUM; 

	return(-0.5*density); 
}

double MCMC_beta_c_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int ll,mm;
	for (l=0;l<D->L;l++){
		ll=D->L+l;
 		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->r_clm[mm]-(para+D_para->r_o_l[l]+D_para->delta_l[l]*D_para->omega_cl[l]);
			SUM+=F*F*exp(D_para->tau_r_cl[ll])/**/+2*log(pnorm(3.5,(para+D_para->r_o_l[l]+D_para->delta_l[l]*D_para->omega_cl[l]),1/pow(exp(D_para->tau_r_cl[ll]),0.5),1,0));
		}
	}	
	F=para-D_priors->beta_mu;	
	density=F*F*D_priors->eta_beta+SUM; 
	return(-0.5*density); 
}

double MCMC_P_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn,ll;
	for (c=0;c<2;c++){
		for (l=0;l<D->L;l++){
		ll=c*D->L+l;
        		for (m=0;m<D->NoORF[ll];m++){
				mm=D->NoSUM[ll]+m;
				for (n=0;n<D->NoTIME[mm];n++){
					nn=c*D->SHIFTlmn+l*D->M*D->N + m*D->N + n;
					F=D->y[nn]-logistic_function_E_JHM(D->x[nn],D_para->K_clm[mm], D_para->r_clm[mm],para);
					SUM+=F*F*exp(D_para->nu_l[ll]);
				}
			}
		}
	}
	F=para-D_priors->P_mu;
	density=F*F*D_priors->eta_P+SUM;
	return(-0.5*density); 
}

/*Data Aug*/

double aug_delta_l_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,int l){
  double u,para,PRODa=0,PRODb=0,F,F1;
  int m,ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D_para->K_clm[mm]-(D_para->alpha_c[1]+D_para->K_o_l[l]);
    F1=D_para->r_clm[mm]-(D_para->beta_c[1]+D_para->r_o_l[l]);
    PRODa+=-0.5*F*F*exp(D_para->tau_K_cl[ll])
      -0.5*F1*F1*exp(D_para->tau_r_cl[ll])
      -log(pnorm(0,(D_para->alpha_c[1]+D_para->K_o_l[l]),1/pow(exp(D_para->tau_K_cl[ll]),0.5),1,0))
      -log(pnorm(3.5,(D_para->beta_c[1]+D_para->r_o_l[l]),1/pow(exp(D_para->tau_r_cl[ll]),0.5),1,0));
  }
  PRODa+=log(1-D_priors->p);

  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D_para->K_clm[mm]-(D_para->alpha_c[1]+D_para->K_o_l[l]+D_para->gamma_cl[l]);
    F1=D_para->r_clm[mm]-(D_para->beta_c[1]+D_para->r_o_l[l]+D_para->omega_cl[l]);
    PRODb+=-0.5*F*F*exp(D_para->tau_K_cl[ll])
      -0.5*F1*F1*exp(D_para->tau_r_cl[ll])
      -log(pnorm(0,(D_para->alpha_c[1]+D_para->K_o_l[l]+D_para->gamma_cl[l]),1/pow(exp(D_para->tau_K_cl[l]),0.5),1,0))
      -log(pnorm(3.5,(D_para->beta_c[1]+D_para->r_o_l[l]+D_para->omega_cl[l]),1/pow(exp(D_para->tau_r_cl[ll]),0.5),1,0));
  }
  PRODb+=log(D_priors->p);

	u=runif(0,1);
  if (log(u/(1-u))<(PRODa-PRODb)){para=0;}
  else{para=1;}
  return(para);
}



double MCMC_K_p_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->K_o_l[l])-exp(para);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_K_o)/D_priors->df)/**/+2*log(trun_const_low_JHM(0,exp(para),exp(D_para->sigma_K_o)));
  }
  F=para-(D_priors->K_mu);
  density=F*F*D_priors->eta_K_p+SUM; 
  return(-0.5*density); 
}

double MCMC_K_o_l_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
	int ll,mm;
	for (c=0;c<2;c++){
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->K_clm[mm]-(D_para->alpha_c[c]+log(para)+c*D_para->delta_l[l]*D_para->gamma_cl[l]);
      			SUM+=F*F*exp(D_para->tau_K_cl[ll])/**/+2*log(pnorm(0,(D_para->alpha_c[c]+log(para)+c*D_para->delta_l[l]*D_para->gamma_cl[l]),1/pow(exp(D_para->tau_K_cl[ll]),0.5),1,0));
    		}
	}
	F=(para)-exp(D_para->K_p);
	density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_K_o)/D_priors->df)+SUM;
  return(-0.5*density); 
}

double MCMC_r_p_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->r_o_l[l])-exp(para);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_r_o)/D_priors->df)/**/+2*log(trun_const_low_JHM(0,exp(para),exp(D_para->sigma_r_o)));
  }
  F=para-(D_priors->r_mu);
  density=F*F*D_priors->eta_r_p+SUM; 
  return(-0.5*density); 
}


double MCMC_r_o_l_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
	int ll,mm;
	for (c=0;c<2;c++){
		ll=c*D->L+l;
		for (m=0;m<D->NoORF[ll];m++){
			mm=D->NoSUM[ll]+m;
			F=D_para->r_clm[mm]-(D_para->beta_c[c]+log(para)+c*D_para->delta_l[l]*D_para->omega_cl[l]);
      			SUM+=F*F*exp(D_para->tau_r_cl[ll])/**/+2*log(pnorm(3.5,(D_para->beta_c[c]+log(para)+c*D_para->delta_l[l]*D_para->omega_cl[l]),1/pow(exp(D_para->tau_r_cl[ll]),0.5),1,0));
    		}
	}
	F=(para)-exp(D_para->r_p);
	density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_r_o)/D_priors->df)+SUM;
  return(-0.5*density); 
}


double MCMC_sigma_K_o_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=exp(D_para->K_o_l[l])-exp(D_para->K_p);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para)/**/+2*log(trun_const_low_JHM(0,exp(D_para->K_p),exp(para)));
  }
  F=(para)-(D_priors->eta_K_o);
  density=F*F*D_priors->psi_K_o+SUM; 
  return(-0.5*density); 
}

double MCMC_sigma_r_o_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=exp(D_para->r_o_l[l])-exp(D_para->r_p);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para)/**/+2*log(trun_const_low_JHM(0,exp(D_para->r_p),exp(para))) ;
  }
  F=(para)-(D_priors->eta_r_o);
  density=F*F*D_priors->psi_r_o+SUM; 
  return(-0.5*density); 
}

double MCMC_gamma_cl_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D_para->K_clm[mm]-(D_para->alpha_c[1]+D_para->K_o_l[l]+D_para->delta_l[l]*log(para));
    SUM+=F*F*exp(D_para->tau_K_cl[ll])/**/+2*log(pnorm(0,(D_para->alpha_c[1]+D_para->K_o_l[l]+D_para->delta_l[l]*log(para)),1/pow(exp(D_para->tau_K_cl[ll]),0.5),1,0));
  }
  F=para-1;
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_gamma)/D_priors->df)+SUM;
  return(-0.5*density);
}

double MCMC_omega_cl_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll=D->L+l,mm;
  for (m=0;m<D->NoORF[ll];m++){
    mm=D->NoSUM[ll]+m;
    F=D_para->r_clm[mm]-(D_para->beta_c[1]+D_para->r_o_l[l]+D_para->delta_l[l]*log(para));
    SUM+=F*F*exp(D_para->tau_r_cl[ll])/**/+2*log(pnorm(3.5,(D_para->beta_c[1]+D_para->r_o_l[l]+D_para->delta_l[l]*log(para)),1/pow(exp(D_para->tau_r_cl[ll]),0.5),1,0));
  }
  F=para-1;
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_omega)/D_priors->df)+SUM;
  return(-0.5*density);
}

double MCMC_sigma_gamma_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->gamma_cl[l])-1;
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para)/**/+2*log(trun_const_low_JHM(0,1,exp(para)));
  }
  F=para-D_priors->eta_gamma;
  density=F*F*D_priors->psi_gamma+SUM;
  return(-0.5*density);
}

double MCMC_sigma_omega_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->omega_cl[l])-1;
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para)/**/+2*log(trun_const_low_JHM(0,1,exp(para)));
  }
  F=para-D_priors->eta_omega;
  density=F*F*D_priors->psi_omega+SUM;
  return(-0.5*density);
}

double MCMC_tau_K_p_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll;
  for (l=0;l<D->L;l++){
    ll=l+c*D->L;
    F=D_para->tau_K_cl[ll]-(para);
    SUM+=F*F*exp(D_para->sigma_tau_K[c])/**/+2*log(1-pnorm(0,(para),1/pow(exp(D_para->sigma_tau_K[c]),0.5),1,0));
  }
  F=para-D_priors->tau_K_mu;
  density=F*F*D_priors->eta_tau_K_p+SUM;
  return(-0.5*density);
}

double MCMC_sigma_tau_K_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,SUM=0,F;
  int ll;
  for (l=0;l<D->L;l++){
    ll=l+c*D->L;
    F=D_para->tau_K_cl[ll]-D_para->tau_K_p[c];
    SUM+=-para+F*F*exp(para)/**/+2*log(1-pnorm(0,(D_para->tau_K_p[c]),1/pow(exp(para),0.5),1,0));
  }
  F=para-D_priors->eta_tau_K;
  density=F*F*D_priors->psi_tau_K+SUM;
  return(-0.5*density);
}
double MCMC_tau_r_p_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l, int m){
  double density,F,SUM=0;
  int ll;
 for (l=0;l<D->L;l++){
   ll=l+c*D->L;
   F=D_para->tau_r_cl[ll]-para;
    SUM+=F*F*exp(D_para->sigma_tau_r[c]);
  }
  F=para-D_priors->tau_r_mu;
  density=F*F*D_priors->eta_tau_r_p+SUM;
  return(-0.5*density);
}

double MCMC_sigma_tau_r_JHM(struct_data_JHM *D,struct_para_JHM *D_para,struct_priors_JHM *D_priors,double para,int c,int l,int m){
  double density,SUM=0,F;
  int ll;
  for (l=0;l<D->L;l++){
    ll=l+c*D->L;
    F=D_para->tau_r_cl[ll]-(D_para->tau_r_p[c]);
    SUM+=-para+F*F*exp(para);
  }
  F=para-D_priors->eta_tau_r;
  density=F*F*D_priors->psi_tau_r+SUM;
  return(-0.5*density);
}



/*Gibbs and MH steps*/
int gibbsandMHloop_JHM(int iter,int thin,struct_data_JHM *D,
  struct_para_JHM *D_para,struct_priors_JHM *D_priors,
  struct_tuning_JHM *D_tuning,struct_adaptive_JHM *D_adaptive,
  int print,int adaptive_phase,int adaptive_period,double *OUT,char **HEADER){
int i,j,l,m,mm,c,ll,*T,t;
T=&t;
*T=0;

if (print==0){
	printheader_JHM(D,HEADER);
}

for (i=0;i<iter;i++){
	for (j=0;j<thin;j++){
		D_para->P=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->P,&D_tuning->P,D_para->P,MCMC_P_JHM,-999,-999,-999);
		D_para->alpha_c[1]=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->alpha_c,&D_tuning->alpha_c,D_para->alpha_c[1],MCMC_alpha_c_JHM,-999,-999,-999);
		D_para->beta_c[1]=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->beta_c,&D_tuning->beta_c,D_para->beta_c[1],MCMC_beta_c_JHM,-999,-999,-999);
		D_para->sigma_gamma=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->sigma_gamma,&D_tuning->sigma_gamma,D_para->sigma_gamma,MCMC_sigma_gamma_JHM,-999,-999,-999);
		D_para->sigma_omega=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->sigma_omega,&D_tuning->sigma_omega,D_para->sigma_omega,MCMC_sigma_omega_JHM,-999,-999,-999);
		D_para->sigma_nu=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->sigma_nu,&D_tuning->sigma_nu,D_para->sigma_nu,MCMC_sigma_nu_JHM,-999,-999,-999);
		D_para->sigma_K_o=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->sigma_K_o,&D_tuning->sigma_K_o,D_para->sigma_K_o,MCMC_sigma_K_o_JHM,-999,-999,-999);
		D_para->sigma_r_o=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->sigma_r_o,&D_tuning->sigma_r_o,D_para->sigma_r_o,MCMC_sigma_r_o_JHM,-999,-999,-999);
		D_para->K_p=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->K_p,&D_tuning->K_p,D_para->K_p,MCMC_K_p_JHM,-999,-999,-999);
		D_para->r_p=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->r_p,&D_tuning->r_p,D_para->r_p,MCMC_r_p_JHM,-999,-999,-999);
		D_para->nu_p=gauss_sample_JHM(D,0,2*D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu_p);
		   																     						     
		for (c=0;c<2;c++){
			D_para->tau_K_p[c]=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->tau_K_p,&D_tuning->tau_K_p,D_para->tau_K_p[c],MCMC_tau_K_p_JHM,c,-999,-999);
			D_para->tau_r_p[c]=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->tau_r_p,&D_tuning->tau_r_p,D_para->tau_r_p[c],MCMC_tau_r_p_JHM,c,-999,-999);
			D_para->sigma_tau_K[c]=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->sigma_tau_K,&D_tuning->sigma_tau_K,D_para->sigma_tau_K[c],MCMC_sigma_tau_K_JHM,c,-999,-999);
			D_para->sigma_tau_r[c]=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->sigma_tau_r,&D_tuning->sigma_tau_r,D_para->sigma_tau_r[c],MCMC_sigma_tau_r_JHM,c,-999,-999);
		}

		for (l=0;l<D->L;l++){       
			D_para->K_o_l[l]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_adaptive->K_o_l,&D_tuning->K_o_l,exp(D_para->K_o_l[l]),MCMC_K_o_l_JHM,-999,l,-999);
			D_para->K_o_l[l]=log(D_para->K_o_l[l]);
			D_para->r_o_l[l]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_adaptive->r_o_l,&D_tuning->r_o_l,exp(D_para->r_o_l[l]),MCMC_r_o_l_JHM,-999,l,-999);
			D_para->r_o_l[l]=log(D_para->r_o_l[l]);
			D_para->gamma_cl[l]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_adaptive->gamma_cl,&D_tuning->gamma_cl,exp(D_para->gamma_cl[l]),MCMC_gamma_cl_JHM,-999,l,-999);
			D_para->gamma_cl[l]=log(D_para->gamma_cl[l]);
			D_para->omega_cl[l]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_adaptive->omega_cl,&D_tuning->omega_cl,exp(D_para->omega_cl[l]),MCMC_omega_cl_JHM,-999,l,-999);
			D_para->omega_cl[l]=log(D_para->omega_cl[l]);
			D_para->delta_l[l]=aug_delta_l_JHM(D,D_para,D_priors,l);
		  
			for (c=0;c<2;c++){
				ll=c*D->L+l;
				D_para->nu_l[ll]=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->nu_l,&D_tuning->nu_l,D_para->nu_l[ll],MCMC_nu_l_JHM,c,l,-999);
				D_para->tau_K_cl[ll]=MCMC_base_truncate_low_JHM(0,D,D_para,D_priors,&D_adaptive->tau_K_cl,&D_tuning->tau_K_cl,D_para->tau_K_cl[ll],MCMC_tau_K_cl_JHM,c,l,-999);
				D_para->tau_r_cl[ll]=MCMC_base_JHM(D,D_para,D_priors,&D_adaptive->tau_r_cl,&D_tuning->tau_r_cl,D_para->tau_r_cl[ll],MCMC_tau_r_cl_JHM,c,l,-999);

			    for (m=0;m<D->NoORF[ll];m++){ 
					mm=D->NoSUM[ll]+m;			
					D_para->K_clm[mm]=MCMC_base_truncate_high_JHM(0,D,D_para,D_priors,&D_adaptive->K_clm,&D_tuning->K_clm,D_para->K_clm[mm],MCMC_K_clm_JHM,c,l,m);
					D_para->r_clm[mm]=MCMC_base_truncate_high_JHM(3.5,D,D_para,D_priors,&D_adaptive->r_clm,&D_tuning->r_clm,D_para->r_clm[mm],MCMC_r_clm_JHM,c,l,m);

					}
			}
		}
			
		if (print==0){
			if (adaptive_phase>0){
				adaptive_phase_process_JHM(D_tuning,D_adaptive,
				  adaptive_period,print,iter);	
			}
		}	
		
	}	
	if (print==1){
		printdata_JHM(D,D_para,OUT,T);
	}
}	
return 0;
}

int adaptive_phase_process_JHM(struct_tuning_JHM *tuning,
  struct_adaptive_JHM *adaptive,int adaptive_period,int print,int iter){
  double explore=1.01,exploredown=1/1.01,ideal_accept_rate=0.25;
    
		if (adaptive->K_clm==1){
			tuning->K_clm=tuning->K_clm*explore;/*more difficult to accept*/
		} else {
			tuning->K_clm=tuning->K_clm*exploredown;/*less difficult to accept*/
		}
		
		if (adaptive->tau_K_cl==1){
			tuning->tau_K_cl=tuning->tau_K_cl*explore;
		} else {
			tuning->tau_K_cl=tuning->tau_K_cl*exploredown;
		}
		
		if (adaptive->r_clm==1){
			tuning->r_clm=tuning->r_clm*explore;
		} else {
			tuning->r_clm=tuning->r_clm*exploredown;
		}
		
		if (adaptive->tau_r_cl==1){
			tuning->tau_r_cl=tuning->tau_r_cl*explore;
		} else {
			tuning->tau_r_cl=tuning->tau_r_cl*exploredown;
		}

		if (adaptive->K_o_l==1){
			tuning->K_o_l=tuning->K_o_l*explore;
		} else {
			tuning->K_o_l=tuning->K_o_l*exploredown;
		}
		
		if (adaptive->sigma_K_o==1){
			tuning->sigma_K_o=tuning->sigma_K_o*explore;
		} else {
			tuning->sigma_K_o=tuning->sigma_K_o*exploredown;
		}

		if (adaptive->r_o_l==1){
			tuning->r_o_l=tuning->r_o_l*explore;
		} else {
			tuning->r_o_l=tuning->r_o_l*exploredown;
		}
		
		if (adaptive->sigma_r_o==1){
			tuning->sigma_r_o=tuning->sigma_r_o*explore;
		} else {
			tuning->sigma_r_o=tuning->sigma_r_o*exploredown;
		}		
 
 		if (adaptive->nu_l==1){
			tuning->nu_l=tuning->nu_l*explore;
		} else {
			tuning->nu_l=tuning->nu_l*exploredown;
		}
		
		if (adaptive->sigma_nu==1){
			tuning->sigma_nu=tuning->sigma_nu*explore;
		} else {
			tuning->sigma_nu=tuning->sigma_nu*exploredown;
		}		

 		if (adaptive->gamma_cl==1){
			tuning->gamma_cl=tuning->gamma_cl*explore;
		} else {
			tuning->gamma_cl=tuning->gamma_cl*exploredown;
		}
		
		if (adaptive->sigma_gamma==1){
			tuning->sigma_gamma=tuning->sigma_gamma*explore;
		} else {
			tuning->sigma_gamma=tuning->sigma_gamma*exploredown;
		}		
		
	 	if (adaptive->omega_cl==1){
			tuning->omega_cl=tuning->omega_cl*explore;
		} else {
			tuning->omega_cl=tuning->omega_cl*exploredown;
		}
		
		if (adaptive->sigma_omega==1){
			tuning->sigma_omega=tuning->sigma_omega*explore;
		} else {
			tuning->sigma_omega=tuning->sigma_omega*exploredown;
		}		

		if (adaptive->alpha_c==1){
			tuning->alpha_c=tuning->alpha_c*explore;
		} else {
			tuning->alpha_c=tuning->alpha_c*exploredown;
		}
		
		if (adaptive->beta_c==1){
			tuning->beta_c=tuning->beta_c*explore;
		} else {
			tuning->beta_c=tuning->beta_c*exploredown;
		}

		if (adaptive->tau_K_p==1){
			tuning->tau_K_p=tuning->tau_K_p*explore;
		} else {
			tuning->tau_K_p=tuning->tau_K_p*exploredown;
		}
		
		if (adaptive->sigma_tau_K==1){
			tuning->sigma_tau_K=tuning->sigma_tau_K*explore;
		} else {
			tuning->sigma_tau_K=tuning->sigma_tau_K*exploredown;
		}		
		
		if (adaptive->tau_r_p==1){
			tuning->tau_r_p=tuning->tau_r_p*explore;
		} else {
			tuning->tau_r_p=tuning->tau_r_p*exploredown;
		}
		
		if (adaptive->sigma_tau_r==1){
			tuning->sigma_tau_r=tuning->sigma_tau_r*explore;
		} else {
			tuning->sigma_tau_r=tuning->sigma_tau_r*exploredown;
		}	

		if (adaptive->K_p==1){
			tuning->K_p=tuning->K_p*explore;
		} else {
			tuning->K_p=tuning->K_p*exploredown;
		}

		if (adaptive->r_p==1){
			tuning->r_p=tuning->r_p*explore;
		} else {
			tuning->r_p=tuning->r_p*exploredown;
		}

		if (adaptive->P/adaptive_period<0.5){
			tuning->P=tuning->P*explore;
		} else {
			tuning->P=tuning->P*exploredown;
		}
	
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





