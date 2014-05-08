#########################################################################
# Demonstration of the separate hierarchical model with a small data set.
# Estimating logistic growth parameters for a wild-type deletion (URA3D) data set at 27C
# 8 replicates of a plate (plate 15).
#########################################################################

require(qfaBayes)
require(parallel)

## load control and query qfa datasets
data("URA3_Raw_extratrim_15")#Control
a<-a_15
data("CDC13-1_Raw_extratrim_15")#Query
b<-b_15

## display experiment variables: Treatment, Screen name and Master plate number
qfa.variables(a)
qfa.variables(b)

## choose experimental variables of interest and filter by them
## for the control and query screen
Treatment_a<-27
Screen_a<-unique(a$Screen.Name)
MPlate_a<-15
Treatment_b<-27
Screen_b<-unique(b$Screen.Name)
MPlate_b<-15
remove_row_a<-c(1,16)
remove_col_a<-c(1,24) 
remove_row_b<-c(1,16)
remove_col_b<-c(1,24) 
SHM_a<-SHM_postpro(a=a,Treatment=Treatment_a,Screen=Screen_a,MPlate=MPlate_a,
  remove_row=remove_row_a,remove_col=remove_col_a)
SHM_b<-SHM_postpro(a=b,Treatment=Treatment_b,Screen=Screen_b,MPlate=MPlate_b,
  remove_row=remove_row_b,remove_col=remove_col_b)

## load SHM specific priors and tuning parameters
data("priors_SHM")
PRIORS<-priors_SHM[[1]]
data("tuning_SHM")
TUNING<-tuning_SHM[[1]]

## select lengths for burn-in, posterior sample and thinning 
## for SHM control and query MCMC
burn<-600000
iters<-1000
thin<-100
adaptive_phase<-1000
adaptive_period<-100

## run MCMC code to produce posterior samples from the SHM
## for the control and query
tasks <- list(
  job1 = function() {require(qfaBayes);SHM_main(burn=burn,iters=iters,thin=thin,
    adaptive_phase=adaptive_phase,adaptive_period=adaptive_period,
    QFA.I=SHM_a$QFA.I,QFA.y=SHM_a$QFA.y,QFA.x=SHM_a$QFA.x,
    QFA.NoORF=SHM_a$QFA.NoORF,QFA.NoTIME=SHM_a$QFA.NoTIME,
    PRIORS=PRIORS,TUNING=TUNING)},
  job2 = function() {require(qfaBayes);SHM_main(burn=burn,iters=iters,thin=thin,
    adaptive_phase=adaptive_phase,adaptive_period=adaptive_period,
    QFA.I=SHM_b$QFA.I,QFA.y=SHM_b$QFA.y,QFA.x=SHM_b$QFA.x,
    QFA.NoORF=SHM_b$QFA.NoORF,QFA.NoTIME=SHM_b$QFA.NoTIME,
    PRIORS=PRIORS,TUNING=TUNING)}
)
cl <- makeCluster(2)
clusterExport(cl=cl, varlist=ls())
out <- clusterApply( 
  cl,
  tasks,
  function(f) f()
)
stopCluster(cl)
SHM_output_a=out[[1]]
SHM_output_b=out[[2]]

## Calculate MDRxMDR using posterior means for the control and query screen
IHM_MDRxMDP=IHM_MDRxMDP_postpro(SHM_a,SHM_output_a,SHM_b,SHM_output_b)

## load IHM specific priors and tuning parameters
data("priors_IHM")
PRIORS_IHM=priors_IHM[[1]]
data("tuning_IHM")
TUNING_IHM=tuning_IHM[[1]]

## select lengths for burn-in, posterior sample and thinning for the IHM
burn_IHM<-500000
iters_IHM<-1000
thin_IHM<-100
adaptive_phase_IHM<-1000
adaptive_period_IHM<-100

## run MCMC code to produce posterior samples from the IHM
IHM_output=IHM_main(burn=burn_IHM,iters=iters_IHM,thin=thin_IHM,
  adaptive_phase=adaptive_phase_IHM,adaptive_period=adaptive_period_IHM,
  QFA.IA=SHM_a$QFA.I,QFA.yA=IHM_MDRxMDP$QFA.yA,QFA.NoORFA=SHM_a$QFA.NoORF,
  QFA.IB=SHM_b$QFA.I,QFA.yB=IHM_MDRxMDP$QFA.yB,QFA.NoORFB=SHM_b$QFA.NoORF,
  PRIORS=PRIORS_IHM,TUNING=TUNING_IHM)

## plot simple IHM fitness plot (control vs query)
plot_IHM_simple(IHM_output,SHM_a)