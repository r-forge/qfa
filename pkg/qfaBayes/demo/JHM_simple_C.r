data("URA3_Raw_extratrim_15")#Control a 
a<-a_15
data("CDC13-1_Raw_extratrim_15")#Query b 
b<-b_15
qfa.variables(a)
qfa.variables(b)

data("priors_JHM")
PRIORS=as.double((priors_JHM)[[1]])
PRIORS[19]=0##
data("tuning_JHM")
TUNING=as.double((tuning_JHM))[1:5]

TreatA=27
Screen_a<-as.character(unique(a$Screen.Name))
MPlate_a<-as.character(unique(a$MasterPlate.Number))

TreatB=27
Screen_b<-as.character(unique(b$Screen.Name))
MPlate_b<-as.character(unique(b$MasterPlate.Number))

JHM<-JHM_postpro(a,TreatA=TreatA,Screen_a=Screen_a,MPlate_a,b,TreatB=TreatB,
  Screen_b,MPlate_b)
JHM_output<-JHM_main(burn=1,iters=1,thin=1,QFA.IA=JHM$QFA.IA,QFA.yA=JHM$QFA.yA,
  QFA.xA=JHM$QFA.xA,QFA.NoORFA=JHM$QFA.NoORFA,QFA.NoTIMEA=JHM$QFA.NoTIMEA,
  QFA.IB=JHM$QFA.IB,QFA.yB=JHM$QFA.yB,QFA.xB=JHM$QFA.xB,
  QFA.NoORFB=JHM$QFA.NoORFB,QFA.NoTIMEB=JHM$QFA.NoTIMEB,PRIORS,TUNING)

ask_plot_simple()
filename="DEMO"
#pdf(paste("JHM_plot_",filename,".pdf",sep=""),useDingbats=F)
plot_JHM_simple(JHM_output,JHM)
#dev.off()

