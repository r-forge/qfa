#include "headers_SHM.h"
#include "functions_SHM.h"
#include "print_SHM.h"

/*trun*/
double trun_const_low(double x,double m, double l)
{
  return(1-0.5-(  (pow(3*l,0.5))*(x-m)/(((l*pow(x-m,2))+3)) +atan((pow(l/3,0.5))*(x-m)))/M_PI);
}


/*Logistic Growth*/

double logistic_function(double t,double K, double r, double P){
	double output;
	output=(K*P*exp(r*t))/(K+P*(exp(r*t)-1));
	return(output);
}

double logistic_function_E(double t,double K, double r, double P){
	double output;
	output=logistic_function(t,exp(K), exp(r), exp(P));
	return(output);
}

/*Gibbs*/

double gauss_sample(struct_data *D ,int start, int N,double x[],double tau,double mu_0,double tau_0){
	double vec,Ndouble,SUM=0;
	int i;
	Ndouble=N;
	for (i=start;i<(start+N);i++){SUM=SUM+x[i];}
		vec=rnorm((tau_0*mu_0+tau*SUM)/(tau_0+Ndouble*tau),1/sqrt(tau_0+Ndouble*tau));/*vec=(tau_0*mu_0+tau*SUM)/(tau_0+Ndouble*tau)+gsl_ran_gaussian(RNG,1/sqrt(tau_0+Ndouble*tau));*/
	return(vec);
}

/*MCMC MH*/

double MCMC_base(struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int),int l, int m){
	double logu,logaprob,can;
		can=rnorm(para,*h);/*can=para+gsl_ran_gaussian(RNG,*h);*/
		logaprob=(*foo)(D,D_para,D_priors,can,l,m)-(*foo)(D,D_para,D_priors,para,l,m); /*logaprob=(*foo)(D,D_para,D_priors,can,l,m)-(*foo)(D,D_para,D_priors,para,l,m);*/

		logu=log(1-runif(0,1));/*logu=log(1-gsl_rng_uniform(RNG));*/
	if (logaprob>logu){para=can;*accept=*accept+1;}
	return(para); 
	}


double MCMC_base_truncate_low(double truncate, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int),int l, int m){
  double logu,logaprob,can;
   can=rnorm(para,*h); /*can=para+gsl_ran_gaussian(RNG,*h);*/
    if(can<=(truncate)){
      can=para;
    }
  logaprob=(*foo)(D,D_para,D_priors,can,l,m)-(*foo)(D,D_para,D_priors,para,l,m);

	logu=log(1-runif(0,1));  /*logu=log(1-gsl_rng_uniform(RNG));*/
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para);
}

double MCMC_base_truncate_high(double truncate, struct_data *D,struct_para *D_para,struct_priors *D_priors,double *accept,double *h,double para,double (*foo)(struct struct_data *D,struct struct_para *D_para,struct struct_priors *D_priors,double,int,int),int l, int m){
  double logu,logaprob,can;
   can=rnorm(para,*h); 
    if(can>=(truncate)){
      can=para;
    }
  logaprob=(*foo)(D,D_para,D_priors,can,l,m)-(*foo)(D,D_para,D_priors,para,l,m);

  logu=log(1-runif(0,1)); 
  if (logaprob>logu){para=can;*accept=*accept+1;}
  return(para);
}



double MCMC_K_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	mm=D->NoSUM[l]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=l*D->M*D->N + m*D->N + n;
		F=D->y[nn]-logistic_function_E(D->x[nn],para, D_para->r_lm[mm],D_para->P);
		SUM+=F*F*exp(D_para->nu_l[l]);
	}	
	F=para-D_para->K_o_l[l];
	density=F*F*exp(D_para->tau_K_l[l])+SUM; 
	return(-0.5*density);
}

double MCMC_r_lm(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	mm=D->NoSUM[l]+m;
	for (n=0;n<D->NoTIME[mm];n++){
		nn=l*D->M*D->N + m*D->N + n;
		F=D->y[nn]-logistic_function_E(D->x[nn], D_para->K_lm[mm],para,D_para->P);	
		SUM+=F*F*exp(D_para->nu_l[l]);
	}	
	F=para-D_para->r_o_l[l];
	density=F*F*exp(D_para->tau_r_l[l])+SUM; 
	return(-0.5*density); 
}

double MCMC_nu_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	for (m=0;m<D->NoORF[l];m++){
		mm=D->NoSUM[l]+m;
		for (n=0;n<D->NoTIME[mm];n++){
			nn=l*D->M*D->N + m*D->N + n;
			F=D->y[nn]-logistic_function_E(D->x[nn],D_para->K_lm[mm],D_para->r_lm[mm],D_para->P);
			SUM+=-para+F*F*exp(para);
		}
	}
	F=para-D_para->nu_p;
	density=F*F*exp(D_para->sigma_nu)+SUM; 
	return(-0.5*density); 
}


double MCMC_sigma_nu(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,SUM=0,F;
	for (l=0;l<D->L;l++){
		F=D_para->nu_l[l]-D_para->nu_p;
		SUM+=-para+F*F*exp(para);
	}	
	F=para-D_priors->eta_nu;
	density=F*F*D_priors->psi_nu+SUM; 
	return(-0.5*density); 
}

double MCMC_P(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
	double density,F,SUM=0;
	int n,mm,nn;
	for (l=0;l<D->L;l++){
		for (m=0;m<D->NoORF[l];m++){
			mm=D->NoSUM[l]+m;
			for (n=0;n<D->NoTIME[mm];n++){
				nn=l*D->M*D->N + m*D->N + n;
				F=D->y[nn]-logistic_function_E(D->x[nn],D_para->K_lm[mm],D_para->r_lm[mm],para);
				SUM+=F*F*exp(D_para->nu_l[l]);
			}
		}
	}
	F=para-D_priors->P_mu;
	density=F*F*D_priors->eta_P+SUM;

	return(-0.5*density); 

}


double MCMC_K_p(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->K_o_l[l])-exp(para);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_K_o)/D_priors->df)/**/+2*log(trun_const_low(0,exp(para),exp(D_para->sigma_K_o)));
  }
  F=para-(D_priors->K_mu);
  density=F*F*D_priors->eta_K_p+SUM; 
  return(-0.5*density); 
}

double MCMC_K_o_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,F,SUM=0;
  int mm;
    for (m=0;m<D->NoORF[l];m++){
      mm=D->NoSUM[l]+m;
      F=D_para->K_lm[mm]-log(para);
      SUM+=F*F*exp(D_para->tau_K_l[l])/**/+2*log(pnorm(0,log(para),1/pow(exp(D_para->tau_K_l[l]),0.5),1,0));
    }
    F=para-exp(D_para->K_p);
    density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_K_o)/D_priors->df)+SUM;
  return(-0.5*density); 
}

double MCMC_r_p(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=exp(D_para->r_o_l[l])-exp(para);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_r_o)/D_priors->df)/**/+2*log(trun_const_low(0,exp(para),exp(D_para->sigma_r_o)));
  }
  F=(para)-(D_priors->r_mu);
  density=F*F*D_priors->eta_r_p+SUM; 
  return(-0.5*density); 
}

double MCMC_r_o_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,F,SUM=0;
  int mm;
    for (m=0;m<D->NoORF[l];m++){
      mm=D->NoSUM[l]+m;
      F=D_para->r_lm[mm]-log(para);
      SUM+=F*F*exp(D_para->tau_r_l[l])/**/+2*log(pnorm(3.5,log(para),1/pow(exp(D_para->tau_r_l[l]),0.5),1,0));
    }
    F=(para)-exp(D_para->r_p);
  density=(1+D_priors->df)*log(1+pow(F,2)*exp(D_para->sigma_r_o)/D_priors->df)+SUM;
  return(-0.5*density); 
}

double MCMC_sigma_K_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=exp(D_para->K_o_l[l])-exp(D_para->K_p);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para)/**/+2*log(trun_const_low(0,exp(D_para->K_p),exp(para)));
  }
  F=(para)-D_priors->eta_K_o;
  density=F*F*D_priors->psi_K_o+SUM;
  return(-0.5*density); 
}

double MCMC_sigma_r_o(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=exp(D_para->r_o_l[l])-exp(D_para->r_p);
    SUM+=(1+D_priors->df)*log(1+pow(F,2)*exp(para)/D_priors->df)-(para)/**/+2*log(trun_const_low(0,exp(D_para->r_p),exp(para))) ;
  }
  F=(para)-D_priors->eta_r_o;
  density=F*F*D_priors->psi_r_o+SUM; 
  return(-0.5*density); 
}
double MCMC_tau_K_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,SUM=0,F;
  int mm;
  for (m=0;m<D->NoORF[l];m++){
    mm=D->NoSUM[l]+m;
    F=D_para->K_lm[mm]-D_para->K_o_l[l];
    SUM+=-para+F*F*exp(para)/**/+2*log(pnorm(0,D_para->K_o_l[l],1/pow(exp(para),0.5),1,0));
  }
  F=para-(D_para->tau_K_p);
  density=F*F*exp(D_para->sigma_tau_K)+SUM;
  return(-0.5*density);
}


double MCMC_tau_r_l(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,SUM=0,F;
  int mm;
  for (m=0;m<D->NoORF[l];m++){
    mm=D->NoSUM[l]+m;
    F=D_para->r_lm[mm]-D_para->r_o_l[l];
    SUM+=-para+F*F*exp(para)/**/+2*log(pnorm(3.5,D_para->r_o_l[l],1/pow(exp(para),0.5),1,0));
  }
  F=para-(D_para->tau_r_p);
  density=F*F*exp(D_para->sigma_tau_r)+SUM;
  return(-0.5*density);
}


double MCMC_tau_K_p(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=D_para->tau_K_l[l]-(para);
    SUM+=F*F*exp(D_para->sigma_tau_K)/**/+2*log(1-pnorm(0,(para),1/pow(exp(D_para->sigma_tau_K),0.5),1,0));
  }
  F=D_priors->tau_K_mu-para;
  density=F*F*D_priors->eta_tau_K_p+SUM;
  return(-0.5*density);
}

double MCMC_sigma_tau_K(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=D_para->tau_K_l[l]-D_para->tau_K_p;
    SUM+=-para+F*F*exp(para)/**/+2*log(1-pnorm(0,(D_para->tau_K_p),1/pow(exp(para),0.5),1,0));
  }
  F=para-D_priors->eta_tau_K;
  density=F*F*D_priors->psi_tau_K+SUM;
  return(-0.5*density);
}

double MCMC_tau_r_p(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l, int m){
  double density,F,SUM=0;
  for (l=0;l<D->L;l++){
    F=D_para->tau_r_l[l]-(para);
    SUM+=F*F*exp(D_para->sigma_tau_r);
  }
  F=D_priors->tau_r_mu-para;
  density=F*F*D_priors->eta_tau_r_p+SUM;
  return(-0.5*density);
}



double MCMC_sigma_tau_r(struct_data *D,struct_para *D_para,struct_priors *D_priors,double para,int l,int m){
  double density,SUM=0,F;
  for (l=0;l<D->L;l++){
    F=D_para->tau_r_l[l]-(D_para->tau_r_p);
    SUM+=-para+F*F*exp(para);
  }
  F=para-D_priors->eta_tau_r;
  density=F*F*D_priors->psi_tau_r+SUM;
  return(-0.5*density);
}

/*Gibbs and MH steps*/
int gibbsandMHloop(int iter,int thin,struct_data *D,struct_para *D_para,
  struct_priors *D_priors,struct_tuning *D_tuning,
  struct_adaptive *D_adaptive,int print,int adaptive_phase,
  int adaptive_period,double *OUT,char **HEADER){
int i,j,l,m,mm,*T,t;
T=&t;
*T=0;
	
	if (print==0){
		printheader(D,HEADER);
	}	
for (i=0;i<iter;i++){
	for (j=0;j<thin;j++){
		D_para->P=MCMC_base(D,D_para,D_priors,&D_adaptive->P,&D_tuning->P,D_para->P,MCMC_P,-999,-999);
		D_para->sigma_K_o=MCMC_base(D,D_para,D_priors,&D_adaptive->sigma_K_o,&D_tuning->sigma_K_o,D_para->sigma_K_o,MCMC_sigma_K_o,-999,-999);
		D_para->sigma_r_o=MCMC_base(D,D_para,D_priors,&D_adaptive->sigma_r_o,&D_tuning->sigma_r_o,D_para->sigma_r_o,MCMC_sigma_r_o,-999,-999);
		D_para->sigma_nu=MCMC_base(D,D_para,D_priors,&D_adaptive->sigma_nu,&D_tuning->sigma_nu,D_para->sigma_nu,MCMC_sigma_nu,-999,-999);
		D_para->K_p=MCMC_base(D,D_para,D_priors,&D_adaptive->K_p,&D_tuning->K_p,D_para->K_p,MCMC_K_p,-999,-999);
		D_para->r_p=MCMC_base(D,D_para,D_priors,&D_adaptive->r_p,&D_tuning->r_p,D_para->r_p,MCMC_r_p,-999,-999);	
		D_para->nu_p=gauss_sample(D,0,D->L,D_para->nu_l,exp(D_para->sigma_nu),D_priors->nu_mu,D_priors->eta_nu_p);
		D_para->tau_K_p=MCMC_base(D,D_para,D_priors,&D_adaptive->tau_K_p,&D_tuning->tau_K_p,D_para->tau_K_p,MCMC_tau_K_p,-999,-999);
		D_para->sigma_tau_K=MCMC_base(D,D_para,D_priors,&D_adaptive->sigma_tau_K,&D_tuning->sigma_tau_K,D_para->sigma_tau_K,MCMC_sigma_tau_K,-999,-999);
		D_para->tau_r_p=MCMC_base(D,D_para,D_priors,&D_adaptive->tau_r_p,&D_tuning->tau_r_p,D_para->tau_r_p,MCMC_tau_r_p,-999,-999);
		D_para->sigma_tau_r=MCMC_base(D,D_para,D_priors,&D_adaptive->sigma_tau_r,&D_tuning->sigma_tau_r,D_para->sigma_tau_r,MCMC_sigma_tau_r,-999,-999);
		  
		for (l=0;l<D->L;l++){
			D_para->tau_K_l[l]=MCMC_base_truncate_low(0,D,D_para,D_priors,&D_adaptive->tau_K_l,&D_tuning->tau_K_l,D_para->tau_K_l[l],MCMC_tau_K_l,l,-999);
			D_para->tau_r_l[l]=MCMC_base(D,D_para,D_priors,&D_adaptive->tau_r_l,&D_tuning->tau_r_l,D_para->tau_r_l[l],MCMC_tau_r_l,l,-999);	  
			D_para->K_o_l[l]=MCMC_base_truncate_low(0,D,D_para,D_priors,&D_adaptive->K_o_l,&D_tuning->K_o_l,exp(D_para->K_o_l[l]),MCMC_K_o_l,l,-999);
			D_para->K_o_l[l]=log(D_para->K_o_l[l]);
			D_para->r_o_l[l]=MCMC_base_truncate_low(0,D,D_para,D_priors,&D_adaptive->r_o_l,&D_tuning->r_o_l,exp(D_para->r_o_l[l]),MCMC_r_o_l,l,-999);
			D_para->r_o_l[l]=log(D_para->r_o_l[l]);
			D_para->nu_l[l]=MCMC_base(D,D_para,D_priors,&D_adaptive->nu_l,&D_tuning->nu_l,D_para->nu_l[l],MCMC_nu_l,l,-999);

			for (m=0;m<D->NoORF[l];m++){ 
				mm=D->NoSUM[l]+m;
				D_para->K_lm[mm]=MCMC_base_truncate_high(0,D,D_para,D_priors,&D_adaptive->K_lm,&D_tuning->K_lm,D_para->K_lm[mm],MCMC_K_lm,l,m);
				D_para->r_lm[mm]=MCMC_base_truncate_high(3.5,D,D_para,D_priors,&D_adaptive->r_lm,&D_tuning->r_lm,D_para->r_lm[mm],MCMC_r_lm,l,m);
			}
		}
		
		if (print==0){
			if (adaptive_phase>0){
				if (adaptive_period>0){
					if (i<adaptive_phase){
						if (i % adaptive_period == 0){
							adaptive_phase_process(D_tuning,D_adaptive,adaptive_period,print,iter);
						}
					}
				}	
			}
		}
	}
	if (print==1){
		printdata(D,D_para,OUT,T);
	}
}
return 0;
}

int adaptive_phase_process(struct_tuning *tuning,
  struct_adaptive *adaptive,int adaptive_period,int print,int iter){
  int explore=1.7,exploredown=0.8,ideal_accept_rate=0.25;
  
  explore=explore+runif(-0.15,0.15);
  
		if (adaptive->K_lm/adaptive_period>ideal_accept_rate){
			tuning->K_lm=adaptive->K_lm*explore;
		} else {
			tuning->K_lm=adaptive->K_lm*exploredown;
		}
		
		if (adaptive->tau_K_l/adaptive_period>ideal_accept_rate){
			tuning->tau_K_l=adaptive->tau_K_l*explore;
		} else {
			tuning->tau_K_l=adaptive->tau_K_l*exploredown;
		}
		
		if (adaptive->r_lm/adaptive_period>ideal_accept_rate){
			tuning->r_lm=adaptive->r_lm*explore;
		} else {
			tuning->r_lm=adaptive->r_lm*exploredown;
		}
		
		if (adaptive->tau_r_l/adaptive_period>ideal_accept_rate){
			tuning->tau_r_l=adaptive->tau_r_l*explore;
		} else {
			tuning->tau_r_l=adaptive->tau_r_l*exploredown;
		}

		if (adaptive->K_o_l/adaptive_period>ideal_accept_rate){
			tuning->K_o_l=adaptive->K_o_l*explore;
		} else {
			tuning->K_o_l=adaptive->K_o_l*exploredown;
		}
		
		if (adaptive->sigma_K_o/adaptive_period>ideal_accept_rate){
			tuning->sigma_K_o=adaptive->sigma_K_o*explore;
		} else {
			tuning->sigma_K_o=adaptive->sigma_K_o*exploredown;
		}

		if (adaptive->r_o_l/adaptive_period>ideal_accept_rate){
			tuning->r_o_l=adaptive->r_o_l*explore;
		} else {
			tuning->r_o_l=adaptive->r_o_l*exploredown;
		}
		
		if (adaptive->sigma_r_o/adaptive_period>ideal_accept_rate){
			tuning->sigma_r_o=adaptive->sigma_r_o*explore;
		} else {
			tuning->sigma_r_o=adaptive->sigma_r_o*exploredown;
		}		
 
 		if (adaptive->nu_l/adaptive_period>ideal_accept_rate){
			tuning->nu_l=adaptive->nu_l*explore;
		} else {
			tuning->nu_l=adaptive->nu_l*exploredown;
		}
		
		if (adaptive->sigma_nu/adaptive_period>ideal_accept_rate){
			tuning->sigma_nu=adaptive->sigma_nu*explore;
		} else {
			tuning->sigma_nu=adaptive->sigma_nu*exploredown;
		}		

		if (adaptive->tau_K_p/adaptive_period>ideal_accept_rate){
			tuning->tau_K_p=adaptive->tau_K_p*explore;
		} else {
			tuning->tau_K_p=adaptive->tau_K_p*exploredown;
		}
		
		if (adaptive->sigma_tau_K/adaptive_period>ideal_accept_rate){
			tuning->sigma_tau_K=adaptive->sigma_tau_K*explore;
		} else {
			tuning->sigma_tau_K=adaptive->sigma_tau_K*exploredown;
		}		
		
		if (adaptive->tau_r_p/adaptive_period>ideal_accept_rate){
			tuning->tau_r_p=adaptive->tau_r_p*explore;
		} else {
			tuning->tau_r_p=adaptive->tau_r_p*exploredown;
		}
		
		if (adaptive->sigma_tau_r/adaptive_period>ideal_accept_rate){
			tuning->sigma_tau_r=adaptive->sigma_tau_r*explore;
		} else {
			tuning->sigma_tau_r=adaptive->sigma_tau_r*exploredown;
		}	

		if (adaptive->K_p/adaptive_period>ideal_accept_rate){
			tuning->K_p=adaptive->K_p*explore;
		} else {
			tuning->K_p=adaptive->K_p*exploredown;
		}

		if (adaptive->r_p/adaptive_period>ideal_accept_rate){
			tuning->r_p=adaptive->r_p*explore;
		} else {
			tuning->r_p=adaptive->r_p*exploredown;
		}

		if (adaptive->P/adaptive_period<0.5){
			tuning->P=adaptive->P*explore;
		} else {
			tuning->P=adaptive->P*exploredown;
		}
	
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
