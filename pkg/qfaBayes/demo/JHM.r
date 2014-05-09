#########################################################################
# Demonstration of the separate hierarchical model with a small data set.
# Estimating logistic growth parameters for a wild-type deletion (URA3D) data set at 27C
# 8 replicates of a plate (plate 15).
#########################################################################

require(qfaBayes)

## load control and query qfa datasets
data("URA3_Raw_extratrim_15")#Control
a<-a_15
a$Expt.Time[a$Expt.Time<0]=0
data("CDC13-1_Raw_extratrim_15")#Query
b<-b_15
b$Expt.Time[b$Expt.Time<0]=0

## display experiment variables: Treatment, Screen name and Master plate number
qfa.variables(a)
qfa.variables(b)

## choose experimental variables of interest and filter by them
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
JHM<-JHM_postpro(a=a,Treatment_a=Treatment_a,Screen_a=Screen_a,
  MPlate_a=MPlate_a,b=b,Treatment_b=Treatment_b,Screen_b=Screen_b,
  MPlate_b=MPlate_b,remove_row_a=remove_row_a,remove_col_a=remove_col_a,
  remove_row_b=remove_row_b,remove_col_b=remove_col_b)

## load JHM specific priors and tuning parameters
data("priors_JHM")
PRIORS<-priors_JHM[[1]]
data("tuning_JHM")
TUNING<-tuning_JHM[[1]]

## select lengths for burn-in, posterior sample and thinning
burn<-800000
iters<-1000
thin<-100
adaptive_phase<-1000
adaptive_period<-100

## run MCMC code to produce posterior samples from the JHM
JHM_output<-JHM_main(burn=burn,iters=iters,thin=thin,adaptive_phase=adaptive_phase,
  adaptive_period=adaptive_period,QFA.IA=JHM$QFA.IA,QFA.yA=JHM$QFA.yA,
  QFA.xA=JHM$QFA.xA,QFA.NoORFA=JHM$QFA.NoORFA,
  QFA.NoTIMEA=JHM$QFA.NoTIMEA,QFA.IB=JHM$QFA.IB,QFA.yB=JHM$QFA.yB,
  QFA.xB=JHM$QFA.xB,QFA.NoORFB=JHM$QFA.NoORFB,QFA.NoTIMEB=JHM$QFA.NoTIMEB,
  PRIORS=PRIORS,TUNING=TUNING)
  
## plot simple JHM fitness plots (control vs query)
plot_JHM_simple(JHM_output,JHM)