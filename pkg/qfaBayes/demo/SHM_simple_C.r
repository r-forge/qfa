data("URA3_Raw_extratrim")
data("priors_SHM")
PRIORS=as.double((priors_SHM)[[1]])[1:18]
qfa.variables(a)

Screen<-unique(a$Screen.Name)
SHM<-SHM_postpro(a=a,Treat=27,Screen=Screen,MPlate=15)
SHM_output<-SHM_main(burn=1,iters=1,thin=1,CAPL=4294,SHM$QFA.I,SHM$QFA.y,SHM$QFA.x,SHM$QFA.NoORF,SHM$QFA.NoTIME,PRIORS)

