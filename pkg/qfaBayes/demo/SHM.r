data("URA3_Raw_extratrim_15")
a<-a_15

qfa.variables(a)

data("priors_SHM")
PRIORS<-priors_SHM[[1]]
data("tuning_SHM")
TUNING<-tuning_SHM[[1]]

Treat<-27
Screen<-unique(a$Screen.Name)
MPlate<-15

burn<-1
iters<-1
thin<-1

SHM<-SHM_postpro(a=a,Treat=Treat,Screen=Screen,MPlate=MPlate)

SHM_output<-SHM_main(burn=burn,iters=iters,thin=thin,QFA.I=SHM$QFA.I,
  QFA.y=SHM$QFA.y,QFA.x=SHM$QFA.x,QFA.NoORF=SHM$QFA.NoORF,
  QFA.NoTIME=SHM$QFA.NoTIME,PRIORS=PRIORS,TUNING=TUNING)

ask_plot_simple()
plot_SHM_simple(SHM_output,SHM)

