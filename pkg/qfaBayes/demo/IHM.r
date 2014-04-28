#########################################################################
# Demonstration of the separate hierarchical model with a small data set.
# Estimating logistic growth parameters for a wild-type deletion (URA3D) data set at 27C
# 8 replicates of a plate (plate 15).
#########################################################################

require(qfaBayes)

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
Treatment_a=27
Screen_a<-unique(a$Screen.Name)
MPlate_a<-15
Treatment_b=27
Screen_b<-unique(b$Screen.Name)
MPlate_b<-15
SHM_a<-SHM_postpro(a=a,Treatment=Treatment_a,Screen=Screen_a,MPlate=MPlate_a)
SHM_b<-SHM_postpro(a=b,Treatment=Treatment_b,Screen=Screen_b,MPlate=MPlate_b)

## load SHM specific priors and tuning parameters
data("priors_SHM")
PRIORS=priors_SHM[[1]]
data("tuning_SHM")
TUNING=tuning_SHM[[1]]

## select lengths for burn-in, posterior sample and thinning 
## for SHM control and query MCMC
burn<-600000
iters<-1000
thin<-100

## run MCMC code to produce posterior samples from the SHM
## for the control and query
SHM_output_a<-SHM_main(burn=burn,iters=iters,thin=thin,QFA.I=SHM_a$QFA.I,
  QFA.y=SHM_a$QFA.y,QFA.x=SHM_a$QFA.x,QFA.NoORF=SHM_a$QFA.NoORF,
  QFA.NoTIME=SHM_a$QFA.NoTIME,PRIORS=PRIORS,TUNING=TUNING)
SHM_output_b<-SHM_main(burn=burn,iters=iters,thin=thin,QFA.I=SHM_b$QFA.I,
  QFA.y=SHM_b$QFA.y,QFA.x=SHM_b$QFA.x,QFA.NoORF=SHM_b$QFA.NoORF,
  QFA.NoTIME=SHM_b$QFA.NoTIME,PRIORS=PRIORS,TUNING=TUNING)

## Calculate posterior means for the control and query screen
SHM_a$QFA.yA=colMeans(SHM_output_a)
SHM_b$QFA.yB=colMeans(SHM_output_b)

## load IHM specific priors and tuning parameters
data("priors_IHM")
PRIORS_IHM=priors_IHM[[1]]
data("tuning_IHM")
TUNING_IHM=tuning_IHM[[1]]

## select lengths for burn-in, posterior sample and thinning for the IHM
burn_IHM<-500000
iters_IHM<-1000
thin_IHM<-100

## run MCMC code to produce posterior samples from the IHM
IHM_output=IHM_main(burn=burn_IHM,iters=iters_IHM,thin=thin_IHM,
  QFA.IA=SHM_a$QFA.I,QFA.yA=SHM_a$QFA.yA,QFA.NoORFA=SHM_a$QFA.NoORF,
  QFA.IB=SHM_b$QFA.I,QFA.yB=SHM_b$QFA.yB,QFA.NoORFB=SHM_b$QFA.NoORF,
  PRIORS=PRIORS_IHM,TUNING=TUNING_IHM)

## plot simple IHM fitness plot (control vs query)
plot_IHM_simple(IHM_output,SHM_a)