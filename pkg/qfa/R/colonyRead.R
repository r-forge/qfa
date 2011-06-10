#### Read in Colonyzer output and return a rod.read-like object ####
# This function essentially replaces the database functionalities of ROD and is therefore quite slow.
# I recommend that you run it once, and write the returned table to file for future use.
# It is *extremely* fast compared to using ROD itself however...
####

colonyzer.read<-function(path=".",files=c(),experiment="ExptDescription.txt",ORF2gene="",libraries="LibraryDescriptions.csv",background=""){
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
	NROW=max(libs$Row)
	NCOL=max(libs$Column)
	NPLATE=max(libs$Plate)
	SPOTARRAY=array("missing",dim=c(NLIB,NPLATE,NROW,NCOL))
	# Fill the spot array object
	for (x in 1:length(libs$Plate)) SPOTARRAY[libdict[[libs[x,"Library"]]],libs[x,"Plate"],libs[x,"Row"],libs[x,"Column"]]=libs[x,"ORF"]
	getORF<-function(lib,plate,row,col) SPOTARRAY[libdict[[lib]],plate,row,col]

	# Open the experimental description
	expt=read.delim(experiment,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	# Watch the mismatch between library names (e.g. BooneSDLV2 and SDL_v2)
	# barcode-> plate dictionaries for all treatments and repeats
	starts=expt$Start.Time
	names(starts)=expt$Barcode
	barcStart<-function(barc) starts[[barc]]

	treats=expt$Treatment
	names(treats)=expt$Barcode
	barcTreat<-function(barc) treats[[barc]]

	meds=expt$Medium
	names(meds)=expt$Barcode
	barcMed<-function(barc) meds[[barc]]

	screens=expt$Screen
	names(screens)=expt$Barcode
	barcScreen<-function(barc) screens[[barc]]

	if(!"RepQuad"%in%colnames(expt)) stop("Woah there!  We need a RepQuad column in expt. description file")
	rquad=expt$RepQuad
	names(rquad)=expt$Barcode
	barcQuad<-function(barc) rquad[[barc]]

	plates=expt$Plate
	names(plates)=expt$Barcode
	barcPlate<-function(barc) plates[[barc]]

	libs=expt$Library
	names(libs)=expt$Barcode
	barcLib<-function(barc) libs[[barc]]

	# Open the ORF2GENE file
	orf2gene=read.delim(ORF2gene,sep="\t",header=FALSE,stringsAsFactors=FALSE)
	colnames(orf2gene)=c("orf","gene")
	orf2gene$orf=toupper(orf2gene$orf)
	# Create an ORF2Gene dictionary
	orfdict=orf2gene$gene
	names(orfdict)=orf2gene$orf
	getGene<-function(orf) orfdict[[orf]]

	bigd<-data.frame()
	for (f in fs){
		# Read in the image analysis output
		iman=read.delim(f,header=FALSE,sep="\t",stringsAsFactors=FALSE)
		# Colonyzer columns
		#colnames(iman)=c("FILENAME","ROW","COLUMN","TOPLEFTX","TOPLEFTY","WHITEAREA","TRIMMED","THRESHOLD","INTENSITY","EDGEPIXELS","COLR","COLG","COLB","BKR","BKG","BKB","EDGELEN","XDIM","YDIM")
		# ROD-like columns
		colnames(iman)=c("Image.Name","Row","Col","X.Offset","Y.Offset","Area","Trimmed","Threshold","Intensity","Edge.Pixels","Colony.Color.R","Colony.Color.G","Colony.Color.B","Background.Color.R","Background.Color.G","Background.Color.B","Edge.length","Tile.Dimensions.X","Tile.Dimensions.Y")

		# Create extra columns
		iman$Growth=iman$Trimmed
		iman$Barcode=substr(iman$Image.Name,1,11)
		iman$Date.Time=substr(iman$Image.Name,13,31)

		# Dump any images which are not in the experimental description file
		iman=iman[iman$Barcode%in%expt$Barcode,]

		# Create a dictionary for filename->photo number
		getPhotoNum<-function(filename){
			# Get plate name from filename
			platename=substr(filename,1,11)
			# Filter iman data frame by filename
			tmp=na.omit(iman[(iman$Barcode==platename)&(iman$Row==1)&(iman$Col==1),])
			tmp=tmp[order(tmp$Image.Name),]
			tmp$PhotoNum=1:length(tmp$Image.Name)
			return(as.numeric(tmp$PhotoNum[tmp$Image.Name==filename]))
		}

		fnames=unique(iman$Image.Name)
		pnum=sapply(fnames,getPhotoNum)
		names(pnum)=fnames
		photoNum<-function(fname) as.numeric(pnum[[fname]])

		iman$Inoc.Time=sapply(iman$Barcode,barcStart)
		iman$Treatments=sapply(iman$Barcode,barcTreat)
		iman$Medium=sapply(iman$Barcode,barcMed)
		iman$Screen.Name=sapply(iman$Barcode,barcScreen)
		iman$RepQuad=sapply(iman$Barcode,barcQuad)
		iman$MasterPlate.Number=sapply(iman$Barcode,barcPlate)
		iman$Timeseries.order=as.numeric(sapply(iman$Image.Name,photoNum))
		iman$Library.Name=sapply(iman$Barcode,barcLib)
		iman$ORF=mapply(getORF, iman$Library, iman$MasterPlate.Number, iman$Row, iman$Col)
		iman$Gene=sapply(iman$ORF,getGene)
		iman$Background=rep(background,length(iman$Image.Name))

		fmt="%Y-%m-%d_%H-%M-%S"
		t0<-as.POSIXlt(as.character(iman$Inoc.Time),format=fmt)
		t1<-as.POSIXlt(as.character(iman$Date.Time),format=fmt)
		iman$Expt.Time=as.numeric(difftime(t1,t0,units="days"))	

		# Bind to existing data
		bigd=rbind(bigd,iman)
	}

	# Print checks so people are sure they're using the correct data #
	print(paste("Number of Barcodes :",length(unique(bigd$Barcode))))
	print(paste("Number of Plate Photos :",length(unique(bigd$Image.Name))))
	print(paste("Number of ORFs :",length(unique(bigd$ORF))))
	print(paste("Number of Culture Images :",length(bigd$Image.Name)))
	print("Treatments :")
	print(unique(bigd$Treatments))
	print("Media:")
	print(unique(bigd$Medium))
	print("Screens:")
	print(unique(bigd$Screen.Name))
	platesize<-max(as.numeric(bigd$Row))*max(as.numeric(bigd$Col))
	if (length(bigd$Date.Time)%%platesize!=0){
		warning("Number of cultures not multiple of plate size")}
	print("Inoculation DateTimes:")
	print(unique(bigd$Inoc.Time))

	return(bigd)

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