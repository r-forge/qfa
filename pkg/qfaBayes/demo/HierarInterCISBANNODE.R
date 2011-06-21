# HierarInter.R

library(qfa, lib=”~/R”)
library(qfaBayes, lib=”~/R”)

### Control ###

Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")
Work="testjoint"
DescripControl<-"ExptDescriptionCDC13.txt"
upd=2000
iter=1000
thin=1
wrk=Work
PlotOutput=FALSE
a=rod.read(file=Control,inoctimes=DescripControl)
qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)

CustomModel="CustomModel7"

ControlFit<-qfa.Hierachical(a,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=FALSE,work="ModelHExample",CustomModel=CustomModel)

### Query ###

Query<-c("cdc13-1_rad9D_SDLv2_Rpt1.txt","cdc13-1_rad9D_SDLv2_Rpt2.txt","cdc13-1_rad9D_SDLv2_Rpt3.txt","cdc13-1_rad9D_SDLv2_Rpt4.txt")
Work="testjoint"
DescripQuery<-"ExptDescriptionCDC13RAD9.txt"
upd=2000
iter=1000
thin=1
wrk=Work
PlotOutput=FALSE
b=rod.read(file=Query,inoctimes=DescripQuery)

qfa.variables(b)
Screen<-as.character(unique(b$Screen.Name))
b<-funcREMOVE(b,Screen,Treat,MPlate)

CustomModel="CustomModel7"
QueryFit<-qfa.Hierachical(b,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=PlotOutput,work="ModelHExample",CustomModel=CustomModel)


### Hierarchical model plots for Control and Query ###
qfaplots.H("Control",ControlFit,LinearGaussian=TRUE)
qfaplots.H("Query",QueryFit,LinearGaussian=TRUE)

### Fit Interaction Model ###
InteractionFit<-qfa.Interaction(ControlFit,QueryFit,iter=300000,upd=200000,thin=200,PlotOutput=FALSE,work="work",CustomModel=FALSE,Priors=FALSE)

### Interaction model plots ###
qfaplots.I("Final",InteractionFit)

# eof


