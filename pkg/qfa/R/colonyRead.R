#### Read in Colonyzer output and return a rod.read-like object ####
# This function essentially replaces the database functionalities of ROD and is therefore quite slow.
# I recommend that you run it once, and write the returned table to file for future use.
# It is *extremely* fast and flexible compared to using ROD itself however...
####

colonyzer.read<-function(path=".",files=c(),experiment="ExptDescription.txt",ORF2gene="",libraries="LibraryDescriptions.csv",screenID=""){
	# Are we reading all files in a folder?
	if (length(strsplit(path,".")[[1]])>1){pathT=TRUE}else{pathT=FALSE}
	# If no Colonyzer files specified, use all .dat in working directory
	if (length(files)==0){fs<-list.files(path=path,pattern=".dat")} else {fs<-files}
	# If we're using a path, then add path to filenames
	if (pathT==TRUE) fs<-paste(path,fs,sep="/")
	print("List of data files to read:")
	print(fs)

	# Open the library descriptions
	libs=read.delim(libraries,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	libs$ORF=toupper(libs$ORF)

	# What library are we using? e.g. library="AndrewsOE_v2" library="SDL_v2"
	liblist=unique(libs$Library)
	NLIB=length(liblist)
	libdict=1:NLIB
	names(libdict)=liblist
	#libs=libs[libs$Library==library,]
	# Make an array for storing info about spots
	NROW=max(libs$Row,na.rm=TRUE)
	NCOL=max(libs$Column,na.rm=TRUE)
	NPLATE=max(libs$Plate,na.rm=TRUE)

	SPOTARRAY=array("missing",dim=c(NLIB,NPLATE,NROW,NCOL))
	# Fill the spot array object
	for (x in 1:length(libs$ORF)) SPOTARRAY[libdict[[libs[x,"Library"]]],libs[x,"Plate"],libs[x,"Row"],libs[x,"Column"]]=libs[x,"ORF"]
	getORF<-function(lib,plate,row,col) SPOTARRAY[libdict[[lib]],plate,row,col]

	# Open the experimental description
	expt=read.delim(experiment,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	# Watch the mismatch between library names (e.g. BooneSDLV2 and SDL_v2)
	# barcode-> plate dictionaries for all treatments and repeats
	barcStart=expt$Start.Time
	names(barcStart)=expt$Barcode

	barcTreat=expt$Treatment
	names(barcTreat)=expt$Barcode

	barcMed=expt$Medium
	names(barcMed)=expt$Barcode

	barcScreen=expt$Screen
	names(barcScreen)=expt$Barcode

	if(!"RepQuad"%in%colnames(expt)) stop("Woah there!  We need a RepQuad column in expt. description file")
	barcQuad=expt$RepQuad
	names(barcQuad)=expt$Barcode

	barcPlate=expt$Plate
	names(barcPlate)=expt$Barcode

	barcLib=expt$Library
	names(barcLib)=expt$Barcode
	
	if("Client"%in%colnames(expt)){
		barcClient=expt$Client
		names(barcClient)=expt$Barcode
	}
	
	if("ExptDate"%in%colnames(expt)){
		barcExptDate=expt$ExptDate
		names(barcExptDate)=expt$Barcode
	}
	
	if("User"%in%colnames(expt)){	
		barcUser=expt$User
		names(barcUser)=expt$Barcode
	}
	
	if("PI"%in%colnames(expt)){
		barcPI=expt$PI
		names(barcPI)=expt$Barcode
	}
	
	if("Condition"%in%colnames(expt)){
		barcCond=expt$Condition
		names(barcCond)=expt$Barcode
	}
	
	if("Inoc"%in%colnames(expt)){
		barcInoc=expt$Inoc
		names(barcInoc)=expt$Barcode
	}
		
	# Open the ORF2GENE file
	orf2gene=read.delim(ORF2gene,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	colnames(orf2gene)=c("orf","gene")
	orf2gene$orf=toupper(orf2gene$orf)
	# Add a "missing" row
	orf2gene=rbind(orf2gene,c("missing","missing"))
	orf2gene=rbind(orf2gene,c("MISSING","MISSING"))
	# Create an ORF2Gene dictionary
	getGene=orf2gene$gene
	names(getGene)=orf2gene$orf
	#colnames(iman)=c("FILENAME","ROW","COLUMN","TOPLEFTX","TOPLEFTY","WHITEAREA","TRIMMED","THRESHOLD","INTENSITY","EDGEPIXELS","COLR","COLG","COLB","BKR","BKG","BKB","EDGELEN","XDIM","YDIM")
	# ROD-like columns
	colNames=c("Image.Name","Row","Col","X.Offset","Y.Offset","Area","Trimmed","Threshold","Intensity","Edge.Pixels","Colony.Color.R","Colony.Color.G","Colony.Color.B","Background.Color.R","Background.Color.G","Background.Color.B","Edge.length","Tile.Dimensions.X","Tile.Dimensions.Y")
	colClasses=c("character",rep("numeric",length(colNames)-1))

	# Read in the image analysis output
	# Have to define colClasses here since from R3.1-10, automatic conversion for string representation of numbers
	if("data.table" %in% rownames(installed.packages())){
		library(data.table) # Faster version of read.delim...
		iman=do.call(rbind, lapply(fs, data.table::fread,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses=colClasses))
		iman=data.frame(iman)
	}else{
		iman=do.call(rbind, lapply(fs, read.delim,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses=colClasses))
	}
	# Sometimes users include multiple copies of the same image analysis files (e.g. in concatenated collection & separately)
	#iman=unique(iman) # This operation is quite time-consuming for large data.frames...
	# Colonyzer columns
	colnames(iman)=colNames

	# Create extra columns
	iman$Growth=iman$Trimmed/(255*iman$Tile.Dimensions.X*iman$Tile.Dimensions.Y)
	if(nchar(iman$Image.Name[1])==31){
		iman$Barcode=substr(iman$Image.Name,1,11)
		iman$Date.Time=substr(iman$Image.Name,13,31)
	}else{
		iman$Barcode=substr(iman$Image.Name,1,15)
		iman$Date.Time=substr(iman$Image.Name,17,35)
	}

	# Dump any images which are not in the experimental description file
	iman=iman[iman$Barcode%in%expt$Barcode,]
	smalliman=iman[(iman$Row==1)&(iman$Col==1),]

	# Create a dictionary for filename->photo number
	getPhotoNum<-function(filename){
		# Get plate name from filename
		if(nchar(iman$Image.Name[1])==31){
			platename=substr(filename,1,11)
		}else{
			platename=substr(filename,1,15)
		}
		# Filter iman data frame by filename
		#tmp=na.omit(smalliman[(smalliman$Barcode==platename),])
		tmp=smalliman[(smalliman$Barcode==platename),]
		tmp=tmp[order(tmp$Image.Name),]
		tmp$PhotoNum=1:length(tmp$Image.Name)
		return(as.numeric(tmp$PhotoNum[tmp$Image.Name==filename]))
	}

	fnames=unique(iman$Image.Name)
	photoNum=sapply(fnames,getPhotoNum)
	names(photoNum)=fnames

	iman$Inoc.Time=barcStart[iman$Barcode]
	iman$Treatments=barcTreat[iman$Barcode]
	iman$Medium=barcMed[iman$Barcode]
	iman$Screen.Name=barcScreen[iman$Barcode]
	iman$RepQuad=barcQuad[iman$Barcode]
	iman$MasterPlate.Number=barcPlate[iman$Barcode]
	iman$Timeseries.order=as.numeric(photoNum[iman$Image.Name])
	iman$Library.Name=barcLib[iman$Barcode]
	iman$ORF=mapply(getORF, iman$Library.Name, iman$MasterPlate.Number, iman$Row, iman$Col)
	iman$Gene=getGene[toupper(iman$ORF)]
	iman$ScreenID=rep(screenID,length(iman$Image.Name))
	if("Client"%in%colnames(expt)){iman$Client=barcClient[iman$Barcode]}
	if("ExptDate"%in%colnames(expt)){iman$ExptDate=barcExptDate[iman$Barcode]}
	if("User"%in%colnames(expt)){iman$User=barcUser[iman$Barcode]}
	if("PI"%in%colnames(expt)){iman$PI=barcPI[iman$Barcode]}
	if("Condition"%in%colnames(expt)){iman$Condition=barcCond[iman$Barcode]}
	if("Inoc"%in%colnames(expt)){iman$Inoc=barcInoc[iman$Barcode]}

	fmt="%Y-%m-%d_%H-%M-%S"
	t0<-as.POSIXlt(as.character(iman$Inoc.Time),format=fmt)
	t1<-as.POSIXlt(as.character(iman$Date.Time),format=fmt)
	iman$Expt.Time=as.numeric(difftime(t1,t0,units="days"))	

	# Print checks so people are sure they're using the correct data #
	print(paste("Number of Barcodes :",length(unique(iman$Barcode))))
	print(paste("Number of Plate Photos :",length(unique(iman$Image.Name))))
	print(paste("Number of ORFs :",length(unique(iman$ORF))))
	print(paste("Number of Culture Images :",length(iman$Image.Name)))
	print("Treatments :")
	print(unique(iman$Treatments))
	print("Media:")
	print(unique(iman$Medium))
	print("Screens:")
	print(unique(iman$Screen.Name))
	if("Client"%in%colnames(iman)){print("Client :"); print(unique(iman$Client))}
	if("ExptDate"%in%colnames(iman)){print("Experiment date :"); print(unique(iman$ExptDate))}
	if("User"%in%colnames(iman)){print("User :"); print(unique(iman$User))}
	if("PI"%in%colnames(iman)){print("PI :"); print(unique(iman$PI))}
	if("Condition"%in%colnames(iman)){print("Condition :"); print(unique(iman$Condition))}
	if("Inoc"%in%colnames(iman)){print("Inoculation type :"); print(unique(iman$Inoc))}	
	platesize<-max(as.numeric(iman$Row))*max(as.numeric(iman$Col))
	if (length(iman$Date.Time)%%platesize!=0){
		warning("Number of cultures not multiple of plate size")}
	print("Inoculation DateTimes:")
	print(unique(iman$Inoc.Time))

	return(iman)
}

# For testing older code, we can write a synthetic ROD object to file using this function.
# Hopefully this will not be required.

rod.write<-function(iman,outf){
# This is how we can construct a ROD file from the returned object...
ROD=data.frame(
	Barcode=iman$Barcode,
	Area=iman$Area,
	SpotRow=iman$Row,
	TrimmedArea=iman$Trimmed,
	SpotColumn=iman$Col,
	Intensity=iman$Intensity,
	EdgePixels=iman$Edge.Pixels,
	Threshold=iman$Threshold,
	XOffset=iman$X.Offset,
	YOffset=iman$Y.Offset,
	Treatments=iman$Treatments,
	Medium=iman$Medium,
	ImageName=iman$Image.Name,
	ORFName=iman$ORF,
	DateOfImage=iman$Date.Time,
	ScreenName=iman$Screen,
	LibraryName=iman$Library,
	MasterPlateNumber=iman$MasterPlate.Number,
	TimeseriesOrder=iman$Timeseries.order,
	ColonyColorR=iman$Colony.Color.R,
	ColonyColorG=iman$Colony.Color.G,
	ColonyColorB=iman$Colony.Color.B,
	BackgroundColorR=iman$Background.Color.R,
	BackgroundColorG=iman$Background.Color.G,
	BackgroundColorB=iman$Background.Color.B,
	Edgelength=iman$Edge.length,
	TileDimensionsX=iman$Tile.Dimensions.X,
	TileDimensionsY=iman$Tile.Dimensions.Y
)

# Write table as rod output
write(file=outf,"ROD Export Version 5")
write(file=outf,append=TRUE,"Barcode\tArea\tSpot Row\tTrimmed Area\tSpot Column\tIntensity\tEdge Pixels\tThreshold\tX Offset\tY Offset\tTreatments\tMedium\tImage Name\tORF Name\tDate of Image\tScreen Name\tLibrary Name\tMasterPlate Number\tTimeseries order\tColony Color R\tColony Color G\tColony Color B\tBackground Color R\tBackground Color G\tBackground Color B\tEdge length\tTile Dimensions X\tTile Dimensions Y")
write.table(ROD,file=outf,append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}
