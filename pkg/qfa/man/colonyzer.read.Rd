\name{colonyzer.read}
\alias{colonyzer.read}
\title{Reads raw cell density timecourse data from Colonyzer output files}
\description{
Reads in and binds together all of the Colonyzer output files in a directory, so they are checked and ready for bayesian or likelihood inference.
Colonyzer is an open source image analysis tool for quantifying cell densities on agar plates: http://research.ncl.ac.uk/colonyzer/
}
\usage{
colonyzer.read(path=".",files=c(),experiment="ExptDescription.txt",
ORF2gene="",libraries="LibraryDescriptions.csv",screenID="")
}
\arguments{
  \item{path}{The path to the folder containing the Colonyzer .dat files to be read: working directory by default.}
  \item{files}{Character vector giving locations of Colonyzer .dat files to be read - overrides path}
  \item{experiment}{Name of text file describing the inoculation times, library and plate number for unique plates.  Taken relative to path if specified.  File must be a tab-delimited text file with no header containing the following columns:
  \itemize{

\item Image.Name - Full name at image capture (includes barcode and date-time) of image from which data are derived
\item Row - Row number (counting from top of image) of culture in rectangular gridded array
\item Col - Column number (counting from left of image) of culture in rectangular gridded array
\item X.Offset - x-coordinate of top left corner of rectangular tile bounding culture (number of pixels from left of image)
\item Y.Offset - y-coordinate of top left corner of rectangular tile bounding culture (number of pixels from top of image)
\item Area - Culture area (pixels)
\item Trimmed - Integrated Optical Density, sum of pixel intensities within culture area
\item Threshold - Global pixel intensity threshold used for image segmentation (after lighting correction)
\item Intensity - Total pixel intensity for square tile containing culture
\item Edge Pixels - Number of pixels classified as culture on edge of square tile
\item Colony.Color.R - Culture red channel intensity
\item Colony.Color.G - Culture green channel intensity
\item Colony.Color.B - Culture blue channel intensity
\item Background.Color.R - Background red channel intensity (for current tile)
\item Background.Color.G - Background green channel intensity (for current tile)
\item Background.Color.B - Background blue channel intensity (for current tile)
\item Edge.length - Number of culture pixels classified as being microcolony edge pixels (useful for classifying contaminants in cultures grown from dilute inoculum)
\item Tile.Dimensions.X - Culture tile width (pixels)
\item Tile.Dimensions.Y - Culture tile height (pixels)
\item Growth - Default measure of cell density (direct copy of one of Trimmed, Threshold or Intensity)
\item Barcode - Unique plate identifier
\item Date.Time - Timestamp of image capture (extracted from image filename)
\item Inoc.Time - User specified date and time of inoculation (specified in ExptDescription.txt file)
\item Treatments - Conditions applied externally to plates (e.g. temperature(s) at which cultures were grown, UV irradiation applied, etc.)
\item Medium - Nutrients/drugs in plate agar
\item Screen.Name - Name of screen (identifies biological repeats, and experiment)
\item RepQuad - Integer identifying which of the quadrants of a 1536 plate were used to inoculate the current 384 plate (set equal to 1 for all cultures for 1536 format for example)
\item MasterPlate Number - Library plate identifier
\item Timeseries order - Sequential photograph number
\item Library.Name	- Name of library, specifying particular culture location
\item ORF - Systematic, unique identifier for genotype in this position in arrayed library
\item Gene - Standard gene name for genotype in this position in arrayed library.  Note that this can be set equal to ORF for example
\item ScreenID - Unique identifier for this QFA screen
\item Client - Client for whom screen was carried out
\item ExptDate - A representative/approximate date for the experiment (note that genome-wide QFA screens typically take weeks to complete)
\item User - Person who actually carried out screen
\item PI - Principal investigator leading project that screen is part of
\item Condition - The most important defining characteristic of screen, as specified by user (e.g. the temperature screen was carried out at if screen is part of multi-temperature set of screens, or the query mutation if part of a set of screens comparing query mutations, or the drugs present in the medium if part of a set of drug screens)
\item Inoc - Qualitative identifier of inoculation type (e.g. "DIL" for dilute inoculum, "CONC" for concentrated).  Used to distinguish between experiments carried out with different methods of inoculation.
\item Expt.Time - Time (days) since user-specified inoculation date (Inoc.Time) that current image was captured

}
}
  \item{ORF2gene}{Path to a tab-delimited text file containing two columns (with no headers) associating unique, systematic strain identifiers (e.g. yeast ORF Y-numbers) with human readable gene names (e.g. standard names from SGD).}
  \item{libraries}{Tab-delimited text file describing the array strains present in each row-column coordinate of each plate in a series of rectangular arrayed libraries.  Header row format is: "Library	ORF 	Plate	Row	Column	Notes".  Columns are:
  \itemize{
  \item Library - Library identifier (e.g. SDLV1)
  \item ORF - Systematic strain identifier
  \item Plate - Plate number
  \item Row - Row number
  \item Column - Column number
  \item Notes - Optional strain notes (e.g. is strain especially sick or missing?)
  }
  }
  \item{screenID}{Unique experiment identifier (e.g. QFA00001)}
}
\value{
An R data.frame where each row corresponds to a single observation on a single colony, with the value of the growth measurement in 'Growth', and the date and time of the measurement in 'Date.Time'. Other information about the observation is stored in the other columns.  Several columns returned are direct copies of Colonyzer output and mapped as follows:
\itemize{
\item Image.Name - Image Name
\item Row - Spot Row
\item Col - Spot Column
\item X.Offset - X Offset
\item Y.Offset - Y Offset
\item Area - Area
\item Trimmed - Trimmed Area
\item Threshold - Threshold
\item Intensity - Intensity
\item Edge.Pixels - Edge Pixels
\item Colony.Color.R - Colony Color R
\item Colony.Color.G - Colony Color G
\item Colony.Color.B - Colony Color B
\item Background.Color.R - Background Color R
\item Background.Color.G - Background Color G
\item Background.Color.B - Background Color B
\item Edge.length - Edge length
\item Tile.Dimensions.X - Tile Dimensions X
\item Tile.Dimensions.Y - Tile Dimensions Y
}

Extra columns are automatically added as follows.  Some of this information is derived from auxiliary files passed to the function such as the experimental description file, the orf-gene dictionary and the library description file:
\itemize{
\item Growth - A cell density surrogate built from trimmed Area normalised by tile area and maximum achievable pixel intensity: Trimmed/(Tile.Dimensions.X*Tile.Dimensions.Y*255)
\item Barcode - Plate identifier, essentially image name with date time and file extension stripped
\item Date.Time - Date time of image capture in YYYY-MM-DD_hh-mm-ss format
\item Inoc.Time - Date time that plate was inoculated.  If plate is grown at a high temperature, date time at which plate was moved into high temperature incubator.  The assumption in this case being that negligible growth occurred before plate temperature was shifted the the target temperature.
\item Treatments - Treatments applied to plate (e.g. temperature)
\item Medium - Medium contained in agar (e.g. nutrients or drugs added to agar)
\item Screen.Name - Unique identifier for experiment (usually identifies repeat number also if multiple repeats carried out).
\item RepQuad - Identifier for experiments scaling down from 1536 format plates to 384, indicating which quadrant on the original 1536 source plate the current 384 format plate belongs to.
\item MasterPlate.Number - Identifies which plate in the source library (as described in the library description file) corresponds to the current plate
\item Timeseries.order - Ordinal describing which photograph captured
\item Library.Name - Identifies which of the libraries identified in the library description file was used to construct this plate
\item ORF - Unique systematic identifier for the genotype of the strain at this location (e.g. yeast Y-number), as defined by library description file
\item Gene - Standard, human readable genotype identifier for the strain at this location, as defined by the ORF-Gene dictionary
\item Background - Tag identifying experiment, typically used to construct file names and axes titles in plots
\item Expt.Time - Number of days passed between inoculation (start of experiment) and current time
}

Finally, as well as returning the object above, this function prints a small report to screen, summarising the data returned.  This includes number of unique barcodes read, number of photos read, number of genotypes in experiment, number of unique culture observations made, a list of treatments applied, a list of media used, a list of unique screen names (e.g. replicates carried out), the plate dimensions (e.g. 1536, 384 or 96 format) and a list of unique inoculation dates.
}
\keyword{qfa}