# Joint.R
# Joint Model
upd=20
iter=10
thin=1
work="Joint"
PlotOutput=FALSE
data("Adam_cdc-1_SDLV2_REP1")
qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))[1]
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)


data("cdc13-1_rad9D_SDLv2_Rpt1")

qfa.variables(b)
Screen<-as.character(unique(b$Screen.Name))[1]
b<-funcREMOVE(b,Screen,Treat,MPlate)
CustomModel="2JCMNNormals2"
JointFit<-qfa.Joint(a,b,Scaling=TRUE,iter,upd,thin,PlotOutput=FALSE,work,CustomModel=FALSE)
save(ControlFit,file=paste(CustomModel,"BASIC","R",sep="."))

qfaplots.J(JointFit,work)

# eof


