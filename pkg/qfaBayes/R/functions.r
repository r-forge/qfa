### Filter by Screen Name, Temprature and Master Plate Number ###
funcREMOVE<-function(data,Screen,Treat,MPlate){
  data=data[data$Screen.Name%in%Screen,]
  data=data[data$Treatments%in%Treat,]
  data=data[data$MasterPlate.Number%in%MPlate,]
  data
}

### Orders dataset ###
funcIDORDER<-function(data){
  data$ID<-paste(data$Barcode,data$MasterPlate.Number,
    formatC(data$Row,digits=2),formatC(data$Col,digits=2),sep="")
  data<-data[order(paste(data$ORF,data$ID),data$Timeseries.order), ]
}

### Gives gene names ###
funcGENE<-function(x,data){
  data$Gene[data$ORF%in%x][1]
}

### Gives number of repeats for each ORF ###
funcNoORF<-function(x,data){
  length(unique((data$ID[data$ORF==x])))
}

### Gives number of time points for a repeat ###
funcNoTime<-function(x,data){
  length((data$ID[data$ID==x]))
}

### Gives running total of number of number of repeats for each ORF ###
funcNoSum<-function(x,NoORF_vec){
  sum(NoORF_vec[1:x])
}

### Adds NA values at the end of a repeat to give consistent row length for an array###
funcRowRep<-function(x,NoTime_vec,data_vec,dimr,dimc){
  c(data_vec[sum(1,NoTime_vec[1:x]):sum(NoTime_vec[1:(x+1)])],
    rep(NA,
    dimc-length(data_vec[sum(1,NoTime_vec[1:x]):sum(NoTime_vec[1:(x+1)])])))
}

### Adds NA values at the end of an ORF to give consistent Col length for an array ###
funcColORF<-function(x,NoSum_vec,data_vec,dimr,dimc){
  c(data_vec[(dimc*NoSum_vec[x]+1):(dimc*NoSum_vec[x+1])],
    rep(NA,
	dimr*dimc-length(data_vec[(dimc*NoSum_vec[x]+1):(dimc*NoSum_vec[x+1])])))
}

### Creates and transposes an Array ###
funcARRAYTRANS<-function(data_vec,dim){
  vec<-array(c(data_vec),dim=dim)
  vec<-aperm(vec, c(2,1,3))
  vec
}

### Sorts data into array with correct dimensions ###
funcXY<-function(data,M,N,NoTime_vec,NoSum_vec,dimr,dimc){
  XY<-unlist(lapply(1:M,funcRowRep,NoTime_vec=NoTime_vec,data_vec=data,dimr,
    dimc))
  XY<-unlist(lapply(1:N,funcColORF,NoSum_vec=NoSum_vec,data_vec=XY,dimr,dimc))
  dim<-c(dimc,dimr,N)
  XY<-funcARRAYTRANS(XY,dim)
  XY
}

### Creates and transposes an Array ###
funcARRAYTRANS_J<-function(data_vec,dim){
  vec<-array(c(data_vec),dim=dim)
  vec<-aperm(vec, c(2,1,3,4))
  vec
}

### Sorts data into array with correct dimensions (Joint Model Specific) ###
funcXY_J<-function(data,data_b,Ma,Mb,N,NoTime_vec,NoSum_vec,NoTime_vec_b,NoSum_vec_b,dimr,dimc){
  XY<-unlist(lapply(1:Ma,funcRowRep,NoTime_vec=NoTime_vec,data_vec=data,
    dimr,dimc))
  XY<-unlist(lapply(1:N,funcColORF,NoSum_vec=NoSum_vec,data_vec=XY,dimr,dimc))
  XY_b<-unlist(lapply(1:Mb,funcRowRep,NoTime_vec=NoTime_vec_b,data_vec=data_b,
    dimr,dimc))
  XY_b<-unlist(lapply(1:N,funcColORF,NoSum_vec=NoSum_vec_b,data_vec=XY_b,dimr,
    dimc))
  dim<-c(dimc,dimr,N,2)
  XY<-funcARRAYTRANS_J(c(XY,XY_b),dim)
  XY
}

### Scales data by maximum theroetical value for IOD ###
funcSCALING<-function(data,vec){
  lim<-max(data$Tile.Dimensions.Y)*max(data$Tile.Dimensions.X)*255
  vec<-vec/lim
  vec
}

###  Gives experiment variables from Colonyzer output###
qfa.variables<-function(data){
  Screen<-as.character(unique(data$Screen.Name))
  Treat<-as.character(unique(data$Treatments))
  MPlate<-unique(data$MasterPlate.Number)
  list(Screen=Screen,Treat=Treat,MPlate=MPlate)
}

#############

