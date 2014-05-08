#include <stdio.h>
#include <math.h>
#include <stdlib.h>    
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

typedef struct struct_data_IHM {
  double *y;
  int L, M,N,*NoORF,*NoSUM,*NoTIME,SHIFTmn,MAXmn;
} struct_data_IHM;

typedef struct struct_para_IHM {
  double
    *alpha_c,
    *delta_l,
    *gamma_cl,	        sigma_gamma,

    *Z_l,              sigma_Z,         
    *nu_cl,             sigma_nu,
    Z_p,
    nu_p;

} struct_para_IHM;

typedef struct struct_priors_IHM {
  double
	Z_mu,			eta_Z_p,
	eta_Z,			psi_Z,

	eta_nu,			psi_nu,
	nu_mu,			eta_nu_p,
	alpha_mu,		eta_alpha,
	p,   
	eta_gamma,		psi_gamma,
    df;
} struct_priors_IHM;

typedef struct struct_tuning_IHM {
  double 
    alpha_c,
    gamma_cl,	        sigma_gamma,

    Z_l,              sigma_Z,         
    nu_cl,             sigma_nu,
    Z_p;

} struct_tuning_IHM;

typedef struct struct_adaptive_IHM {
  double 
    alpha_c,
    delta_l,
    gamma_cl,	        sigma_gamma,

    Z_l,              sigma_Z,         
    nu_cl,             sigma_nu,
    Z_p;

} struct_adaptive_IHM;
