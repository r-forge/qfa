# Hierarchical.R
# Single Model

require(qfaBayes)
Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")[1]
Query<-c("cdc13-1_rad9D_SDLv2_Rpt1.txt","cdc13-1_rad9D_SDLv2_Rpt2.txt","cdc13-1_rad9D_SDLv2_Rpt3.txt","cdc13-1_rad9D_SDLv2_Rpt4.txt")[1]
Work="testjoint"
DescripControl<-"ExptDescriptionCDC13.txt"
DescripQuery<-"ExptDescriptionCDC13RAD9.txt"
upd=200000
iter=100000
thin=100
wrk=Work
Cont=Control
Que=Query
PlotOutput=FALSE
DescripCont=DescripControl
DescripQue=DescripQuery

data("Adam_cdc-1_SDLV2_REP1")
qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))[1]
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)

CustomModel="CustomModel3"
BASIC<-qfa.Hierachical(a,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=PlotOutput,work="ModelHExample",CustomModel=CustomModel)
save(CustomModel,file=paste(CustomModel,"R",sep="."))#####
QFA.H.Plots(CustomModel,BASIC)

# eof

