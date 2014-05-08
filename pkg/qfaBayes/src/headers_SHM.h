#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

typedef struct struct_data {
  double *y, *x;
  int L, M,N,*NoORF,*NoSUM,*NoTIME,SHIFTlmn,maxy,maxNoTIME;
} struct_data;

typedef struct struct_para {
  double
    *K_lm,            *tau_K_l,
    *r_lm,            *tau_r_l,

    *K_o_l,            sigma_K_o,
    *r_o_l,            sigma_r_o,
    *nu_l,             sigma_nu,

    K_p,
    r_p,
    nu_p,
    P,
    tau_K_p,tau_r_p,sigma_tau_K,sigma_tau_r;
} struct_para;

typedef struct struct_priors {
  double
    eta_K_o,                psi_K_o,
    eta_r_o,                psi_r_o,
    eta_nu,                 psi_nu,

    K_mu,                   eta_K_p,
    r_mu,                   eta_r_p,
    nu_mu,                  eta_nu_p,
    P_mu,                   eta_P,
    df,
    tau_K_mu,tau_r_mu,
    eta_tau_K,eta_tau_r,
    eta_tau_K_p,eta_tau_r_p,
    psi_tau_K,psi_tau_r;

} struct_priors;

typedef struct struct_tuning {
  double      
    K_lm,            tau_K_l,
    r_lm,            tau_r_l,

    K_o_l,            sigma_K_o,
    r_o_l,            sigma_r_o,
    nu_l,             sigma_nu,

    K_p,
    r_p,
    P,
    tau_K_p,tau_r_p,sigma_tau_K,sigma_tau_r;
} struct_tuning;

typedef struct struct_adaptive {
  double 
   
    K_lm,            tau_K_l,
    r_lm,            tau_r_l,

    K_o_l,            sigma_K_o,
    r_o_l,            sigma_r_o,
    nu_l,             sigma_nu,

    K_p,
    r_p,
    P,
    tau_K_p,tau_r_p,sigma_tau_K,sigma_tau_r;
} struct_adaptive;
