# Joint.R
# Joint Model
Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")
Work="testjoint"
DescripControl<-"ExptDescriptionCDC13.txt"
upd=2000
iter=1000
thin=1
work="Joint"
PlotOutput=FALSE
a=rod.read(file=Control,inoctimes=DescripControl)
qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))[1]
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)


Query<-c("cdc13-1_rad9D_SDLv2_Rpt1.txt","cdc13-1_rad9D_SDLv2_Rpt2.txt","cdc13-1_rad9D_SDLv2_Rpt3.txt","cdc13-1_rad9D_SDLv2_Rpt4.txt")
Work="testjoint"
DescripQuery<-"ExptDescriptionCDC13RAD9.txt"
upd=2
iter=5
thin=1
wrk=Work
PlotOutput=FALSE
b=rod.read(file=Query,inoctimes=DescripQuery)
qfa.variables(b)
Screen<-as.character(unique(b$Screen.Name))[1]
b<-funcREMOVE(b,Screen,Treat,MPlate)

JointFit<-qfa.Joint(a,b,Scaling=TRUE,iter,upd,thin,PlotOutput=FALSE,work,CustomModel=FALSE)
QFA.J.Plots(work,JointFit)

# eof

