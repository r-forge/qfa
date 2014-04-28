data("URA3_Raw_extratrim_15")
a<-a_15
data("CDC13-1_Raw_extratrim_15")
b<-b_15
qfa.variables(a)
qfa.variables(b)

data("priors_SHM")
PRIORS=as.double((priors_SHM)[[1]])
data("tuning_SHM")
TUNING=as.double((tuning_SHM))[1:4]
data("tuning_IHM")
TUNING_IHM=as.double((tuning_IHM))[1:8]

Screen_a<-unique(a$Screen.Name)
SHM_a<-SHM_postpro(a=a,Treat=27,Screen=Screen_a,MPlate=15)
SHM_output_a<-SHM_main(burn=1,iters=1,thin=1,QFA.I=SHM_a$QFA.I,
  QFA.y=SHM_a$QFA.y,QFA.x=SHM_a$QFA.x,QFA.NoORF=SHM_a$QFA.NoORF,
  QFA.NoTIME=SHM_a$QFA.NoTIME,PRIORS=PRIORS,TUNING=TUNING)

Screen_b<-unique(b$Screen.Name)
SHM_b<-SHM_postpro(a=b,Treat=27,Screen=Screen_b,MPlate=15)
SHM_output_b<-SHM_main(burn=1,iters=1,thin=1,QFA.I=SHM_b$QFA.I,
  QFA.y=SHM_b$QFA.y,QFA.x=SHM_b$QFA.x,QFA.NoORF=SHM_b$QFA.NoORF,
  QFA.NoTIME=SHM_b$QFA.NoTIME,PRIORS=PRIORS,TUNING=TUNING)

data("priors_IHM")
PRIORS_IHM=as.double((priors_IHM)[[1]])
data("tuning_SHM")
TUNING=as.double((tuning_SHM))[1:8]

SHM_a$QFA.yA=colMeans(SHM_output_a)
SHM_b$QFA.yB=colMeans(SHM_output_b)

IHM_output=IHM_main(burn=1,iters=1,thin=1,QFA.IA=SHM_a$QFA.I,
  QFA.yA=SHM_a$QFA.yA,QFA.NoORFA=SHM_a$QFA.NoORF,QFA.IB=SHM_b$QFA.I,
  QFA.yB=SHM_b$QFA.yB,QFA.NoORFB=SHM_b$QFA.NoORF,PRIORS=PRIORS_IHM,
  TUNING=TUNING_IHM)

ask_plot_simple()
filename="DEMO"
#pdf(paste("IHM_plot_",filename,".pdf",sep=""),useDingbats=F)
plot_IHM_simple(IHM_output,SHM_a)
#dev.off()
