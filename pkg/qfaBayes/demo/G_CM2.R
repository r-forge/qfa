# Hierarchical.R
# Single Model
CustomModel="C_CM2"

Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")
DescripControl<-"ExptDescriptionCDC13.txt"
source("iter.R")
#upd=8
#iter=2
#thin=2
work=paste(upd,iter,thin,sep="_")
data("AdamFull")
qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))[1]
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)

inits=list('K'=log(0.089),'r'=log(1.82),'PO_L'=log(0.00224),'K_i_tau_L'=1.25,'r_i_tau_L'=3.429)

ControlFit<-qfa.Hierachical(a,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=FALSE,work="ModelHExample",CustomModel=CustomModel,inits=inits)
save(ControlFit,file=paste(CustomModel,work,"R",sep="."))

### Plots ###
qfaplots.H(ControlFit,CustomModel,LinearGaussian=TRUE)

# eof


