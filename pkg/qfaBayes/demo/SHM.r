#########################################################################
## Demonstration of the separate hierarchical model with a small data set.
## Estimating logistic growth parameters for a wild-type deletion (URA3D) data set at 27C
## 8 replicates of a plate (plate 15).
#########################################################################

require(qfaBayes)

## load qfa dataset
data("URA3_Raw_extratrim_15")
a<-a_15

## display experiment variables: Treatment, Screen name and Master plate number
qfa.variables(a)

## choose experimental variables of interest and filter by them
Treatment<-27
Screen<-unique(a$Screen.Name)
MPlate<-15
remove_row<-c(1,16)
remove_col<-c(1,24) 
SHM<-SHM_postpro(a=a,Treatment=Treatment,Screen=Screen,MPlate=MPlate,
  remove_row=remove_row,remove_col=remove_col)

## load SHM specific priors and tuning parameters
data("priors_SHM")
PRIORS<-priors_SHM[[1]]
data("tuning_SHM")
TUNING<-tuning_SHM[[1]]

## select lengths for burn-in, posterior sample and thinning
burn<-800000
iters<-1000
thin<-100
adaptive_phase<-1000
adaptive_period<-100

## run MCMC code to produce posterior samples from the SHM
SHM_output<-SHM_main(burn=burn,iters=iters,thin=thin,
  adaptive_phase=adaptive_phase,adaptive_period=adaptive_period,
  QFA.I=SHM$QFA.I,QFA.y=SHM$QFA.y,QFA.x=SHM$QFA.x,
  QFA.NoORF=SHM$QFA.NoORF,QFA.NoTIME=SHM$QFA.NoTIME,
  PRIORS=PRIORS,TUNING=TUNING)

## plot simple SHM logistic growth plots
plot_SHM_simple(SHM_output,SHM)