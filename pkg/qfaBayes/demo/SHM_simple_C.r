data("URA3_Raw_extratrim_15")
data("priors_SHM")
PRIORS=as.double((priors_SHM)[[1]])[1:18]
qfa.variables(a)

Screen<-unique(a$Screen.Name)
SHM<-SHM_postpro(a=a,Treat=27,Screen=Screen,MPlate=15)
SHM_output<-SHM_main(burn=1,iters=1,thin=1,CAPL=4294,QFA.I=SHM$QFA.I,QFA.y=SHM$QFA.y,QFA.x=SHM$QFA.x,QFA.NoORF=SHM$QFA.NoORF,QFA.NoTIME=SHM$QFA.NoTIME,PRIORS=PRIORS)

