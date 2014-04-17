data("URA3_Raw_extratrim_15")
a<-a_15

data("priors_SHM")
PRIORS=as.double((priors_SHM)[[1]])[1:18]
data("tuning_SHM")
TUNING=as.double((tuning_SHM))[1:8]
qfa.variables(a)

Screen<-unique(a$Screen.Name)

SHM<-SHM_postpro(a=a,Treat=27,Screen=Screen,MPlate=15)
SHM_output<-SHM_main(burn=1,iters=1,thin=1,CAPL=50,QFA.I=SHM$QFA.I,
  QFA.y=SHM$QFA.y,QFA.x=SHM$QFA.x,QFA.NoORF=SHM$QFA.NoORF,
  QFA.NoTIME=SHM$QFA.NoTIME,PRIORS=PRIORS,TUNING=TUNING)

ask_plot_simple()
filename="DEMO"
#pdf(paste("SHM_plot_",filename,".pdf",sep=""),useDingbats=F)
plot_SHM_simple(SHM_output,SHM)
#dev.off()
