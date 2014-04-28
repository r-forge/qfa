data("URA3_Raw_extratrim_15")#Control a 
a<-a_15
data("CDC13-1_Raw_extratrim_15")#Query b 
b<-b_15

qfa.variables(a)
qfa.variables(b)

data("priors_JHM")
PRIORS=priors_JHM[[1]]
data("tuning_JHM")
TUNING=tuning_JHM[[1]]

Treat_a=27
Screen_a<-unique(a$Screen.Name)
MPlate_a<-15

Treat_b=27
Screen_b<-unique(b$Screen.Name)
MPlate_b<-15

burn<-1
iters<-1
thin<-1

JHM<-JHM_postpro(a,TreatA=Treat_a,Screen_a=Screen_a,MPlate_a,b,TreatB=Treat_b,
  Screen_b,MPlate_b)
  
JHM_output<-JHM_main(burn=burn,iters=iters,thin=thin,QFA.IA=JHM$QFA.IA,
  QFA.yA=JHM$QFA.yA,QFA.xA=JHM$QFA.xA,QFA.NoORFA=JHM$QFA.NoORFA,
  QFA.NoTIMEA=JHM$QFA.NoTIMEA,QFA.IB=JHM$QFA.IB,QFA.yB=JHM$QFA.yB,
  QFA.xB=JHM$QFA.xB,QFA.NoORFB=JHM$QFA.NoORFB,QFA.NoTIMEB=JHM$QFA.NoTIMEB,
  PRIORS=PRIORS,TUNING=TUNING)

ask_plot_simple()
plot_JHM_simple(JHM_output,JHM)


