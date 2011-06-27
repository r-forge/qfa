# Hierarchical.R
# Single Model

Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")[1]
DescripControl<-"ExptDescriptionCDC13.txt"
upd=2000000
iter=1000000
thin=100
data("AdamFull")

qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))[1]
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)

CustomModel="CustomModel10"
ControlFit<-qfa.Hierachical(a,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=FALSE,work="ModelHExample",CustomModel=CustomModel)
save(ControlFit,file=paste(CustomModel,"R","3",sep="."))

### Plots ###
qfaplots.H(ControlFit,CustomModel,LinearGaussian=TRUE)

# eof


