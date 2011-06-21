# Hierarchical.R
# Single Model

Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")[1]
DescripControl<-"ExptDescriptionCDC13.txt"
upd=200000
iter=100000
thin=100
data("Adam_cdc-1_SDLV2_REP1")

qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))[1]
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)

CustomModel="CustomModel7"
ControlFit<-qfa.Hierachical(a,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=FALSE,work="ModelHExample",CustomModel=CustomModel)

### Plots ###
qfaplots.H(CustomModel,ControlFit,LinearGaussian=TRUE)

# eof


