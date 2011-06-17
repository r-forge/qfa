# Hierarchical.R
# Single Model

require(qfaBayes)
Control<-c("Adam_cdc13-1_SDLV2_REP1.txt","Adam_cdc13-1_SDLV2_REP2.txt","Adam_cdc13-1_SDLV2_REP3.txt","Adam_cdc13-1_SDLV2_REP4.txt")
Query<-c("cdc13-1_rad9D_SDLv2_Rpt1.txt","cdc13-1_rad9D_SDLv2_Rpt2.txt","cdc13-1_rad9D_SDLv2_Rpt3.txt","cdc13-1_rad9D_SDLv2_Rpt4.txt")
Work="testjoint"
DescripControl<-"ExptDescriptionCDC13.txt"
DescripQuery<-"ExptDescriptionCDC13RAD9.txt"
upd=10
iter=1
thin=1
wrk=Work
Cont=Control
Que=Query
PlotOutput=FALSE
DescripCont=DescripControl
DescripQue=DescripQuery

a=rod.read(file=Cont,inoctimes=DescripCont)
qfa.variables(a)
Screen<-as.character(unique(a$Screen.Name))
Treat<-as.character(unique(a$Treatment))[2]
MPlate<-unique(a$MasterPlate.Number)[15]
a<-funcREMOVE(a,Screen,Treat,MPlate)

CustomModel="CustomModel7"
BASICa<-qfa.Hierachical(a,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=PlotOutput,work="ModelHExample",CustomModel=CustomModel)
save(BASICa,file=paste("BASICa","R",sep="."))#####
###################
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

b=rod.read(file=Que,inoctimes=DescripQue)

qfa.variables(b)
Screen<-as.character(unique(b$Screen.Name))
b<-funcREMOVE(b,Screen,Treat,MPlate)

CustomModel="CustomModel7"
BASICb<-qfa.Hierachical(b,Scaling=TRUE,iter=iter,upd=upd,thin=thin,PlotOutput=PlotOutput,work="ModelHExample",CustomModel=CustomModel)
save(BASICb,file=paste("BASICb","R",sep="."))#####

#######################################

BASIC_I<-qfa.Interaction(BASICa,BASICb,iter=10000,upd=10000,thin=100,PlotOutput=FALSE,work="work",CustomModel=FALSE,Priors=FALSE)

QFA.I.Plots("Final",BASIC_I)


# eof


