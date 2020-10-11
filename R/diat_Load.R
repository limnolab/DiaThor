#' The function 'loadFile' loads the CSV or dataframe file, sets the Output folder for the package, and conducts both an exact and an heuristic search of the species' names
#' The function 'findAcronym' conducts both an exact and an heuristic search of the species' names and tries to convert it to its acronym in the internal database
#' @param species_df The data frame with your species data. Species as rows, Sites as columns. If empty, a dialog will prompt for a CSV file
#' @param isRelAb Boolean. If set to 'TRUE' it means that your species' data is the relative abundance of each species per site. If FALSE, it means that it the data corresponds to absolute densities. Default = FALSE
#' @param maxDistTaxa Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2
#' @description
#' The input file for the package is a dataframe or an external CSV file. Species should be listed as rows, with species' names in column 1.
#' The following columns have to be the abundance of each species (relative or absolute) in each sample (column).
#' The first row of the file has to contain the headers with the sample names.
#' If a dataframe is not specified as a parameter (species_df), the package will show a dialog box to search for the CSV file
#' A second dialog box will help set up an Output folder, where all outputs from the package will be exported to (dataframes, CSV files, plots in PDF)
#' The package downloads and installs  a wrapper for the Diat.Barcode project. Besides citing this package, the Diat.Barcode project should also be cited if the package is used, as follows:
#' \itemize{
#' \item Rimet, Frederic; Gusev, Evgenuy; Kahlert, Maria; Kelly, Martyn; Kulikovskiy, Maxim; Maltsev, Yevhen; Mann, David; Pfannkuchen, Martin; Trobajo, Rosa; Vasselon, Valentin; Zimmermann, Jonas; Bouchez, Agnès, 2018, "Diat.barcode, an open-access barcode library for diatoms", https://doi.org/10.15454/TOMBYZ, Portail Data Inra, V1
#' }
#' Sample data is taken from:
#' \itemize{
#' \item Frederic Rimet; Philippe Chaumeil; Francois Keck; Lenaïg Kermarrec; Valentin Vasselon; Maria Kahlert; Alain Franc; Agnes Bouchez. 2016. R-Syst::diatom: an open-access and curated barcode database for diatoms and freshwater monitoring. Database, 2016, 1-21.http://database.oxfordjournals.org/content/2016/baw016.full?keytype=ref&ijkey=H324uA95JzzEomz
#' }
#' \itemize{
#' \item Rimet F., Chaumeil P., Keck F., Kermarrec L., Vasselon V., Kahlert M., Franc A., Bouchez A., 2015. R-Syst::diatom: a barcode database for diatoms and freshwater biomonitoring - data sources and curation procedure. INRA Report, 14 pages http://dx.doi.org/10.5281/zenodo.31137
#' }
#' @examples
#' # EXAMPLE 1:
#' data("environmental_data") #ESTO HAY QUE HACER ALGUN EJEMPLO EMPAQUETADO DENTRO DEL PAQUETE
#' data("species_data")#ESTO HAY QUE HACER ALGUN EJEMPLO EMPAQUETADO DENTRO DEL PAQUETE
#' # EXAMPLE 2: Loads sample data where species are in absolute densities
#' data("environmental_data_example2")#ESTO HAY QUE HACER ALGUN EJEMPLO EMPAQUETADO DENTRO DEL PAQUETE
#' data("species_data_example2")#ESTO HAY QUE HACER ALGUN EJEMPLO EMPAQUETADO DENTRO DEL PAQUETE
#' # Calculates everything
#' thispackage::calculate_everything(species_df)
#' @concepts ecology, diatom, bioindicator, biotic indices
#' @export diat_findAcronyms
#' @export diat_loadData


########-------- FUNCTION TO FIND ACRONYMS FROM DATAFRAME  --------------#########
### INPUT: CSV file with the list of species in the "species" column, or dataframe with the same column
### OUTPUT: dataframe with an additional "acronym" column. NA = acronym not found

diat_findAcronyms <- function(species_df, maxDistTaxa=2){
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(species_df)) {
    print("Select CSV matrices with your sample data")
    Filters <- matrix(c("Comma Separated Values (CSV)", "*.csv"),
                      1, 2, byrow = TRUE)
    species_df <- as.data.frame(read.csv(file.choose())) #select the data matrix with a dialog box
    #handles cancel button
    if (missing(species_df)){
      stop("Calculations cancelled")
    }
  }

  acronymDB <- as.data.frame(DiaThor::acronyms)
  #acronymDB <- read.csv("../Indices/acronyms.csv") #uses the external csv file
  acronymDB$species <- as.character(acronymDB$species)
  acronymDB$acronym <- as.character(acronymDB$acronym)
  acronymDB$species <- trimws(acronymDB$species)
  acronymDB$acronym <- trimws(acronymDB$acronym)
  species_df$acronym <- NA
  species_df$acronym <- as.character(species_df$acronym)
  species_df$species <- as.character(species_df$species)

  print ("Searching for acronyms in exact match, please wait...")
  #exact match by species name (fast)
  species_df$acronym <- acronymDB$acronym[match(rownames(species_df), acronymDB$species)]
  #with the unmatched exact, try heuristic
  print ("Searching for acronyms in heuristic match, please wait...")
  #PROGRESS BAR
  ntot <- sum(is.na(species_df$acronym))
  pb <- txtProgressBar(min = 1, max = ntot, style = 3)
  for(i in 1:ntot) {
    if (is.na(species_df$acronym[i])){
      species_df$acronym[i] <- acronymDB$acronym[stringdist::amatch(species_df[i,"species"], trimws(acronymDB$species), maxDist=maxDistTaxa, matchNA = F)]
      setTxtProgressBar(pb, i)
    }
  }
  close(pb)

  print(paste("Acronyms were found for ", length(species_df$species) - sum(is.na(species_df$acronym)) , " species out of ", length(species_df$acronym), sep=""))
  return(species_df)
}

########-------- FUNCTION TO LOAD FILE  --------------#########
#### IN THIS SECTION WE READ THE MATRIX
### INPUT: CSV file with the samples * abundances or relative abundances
### OUTPUTS: a dataframe with the species as matched against the database, with species in RA; Taxaincluded and Taxaexcluded in CSV in the Output
### folder, detailing which taxa were recognized and which were not

diat_loadData <- function(species_df, isRelAb=FALSE, maxDistTaxa=2){
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(species_df)) {
    print("Select CSV matrices with your sample data")
    Filters <- matrix(c("Comma Separated Values (CSV)", "*.csv"),
                      1, 2, byrow = TRUE)
    species_df <- as.data.frame(read.csv(file.choose())) #select the data matrix with a dialog box
    #handles cancel button
    if (is.na(species_df)){
      stop("Calculations cancelled, no input data")
    }
  }
  print("Select Results folder")
  resultsPath <- choose.dir(default = "", caption = "Select folder for your Results")
  if (is.na(resultsPath)){
    stop("Calculations cancelled, no folder selected")
  }

  print("Result folder selected")
  if("species" %in% colnames(species_df)) {
    row.names(species_df)<-species_df[,"species"] #converts the species names to rownames
    species_df<- species_df[ , !(names(species_df) == "species")] #removes the species names column
  } else {
    #check if species names are in the first column
    stop("Calculations cancelled. File should contain a column named 'species' with the species list")
  }
  species_df <- species_df[order(row.names(species_df)),] #reorder alphabetically
  #if there is an acronym column, it removes it and stores it for later
  if("acronym" %in% colnames(species_df)) {
    print("acronyms found")
  } else {
    #IF THERE IS NO ACRONYM COLUMN, IT SHOULD BUILD IT, and search for the acronyms
    species_df$acronym <- NA
    acronymDB <- as.data.frame(DiaThor::acronyms)
    #acronymDB <- read.csv("../Indices/acronyms.csv") #uses the external csv file
    acronymDB$species <- as.character(acronymDB$species)
    acronymDB$species <- trimws(acronymDB$species)
    acronymDB$acronym <- as.character(acronymDB$acronym)
    acronymDB$acronym <- trimws(acronymDB$acronym)
    #species_df$species <- as.character(species_df$species) #species are now rownames
    print ("Searching for acronyms in exact match, please wait...")
    #exact match by species name (fast)
    species_df$acronym <- acronymDB$acronym[match(rownames(species_df), acronymDB$species)]
    #with the unmatched exact, try heuristic
    print ("Searching for acronyms in heuristic match, please wait...")
    #PROGRESS BAR FOR HEURISTIC SEARCH OF ACRONYMS
    ntot <- sum(is.na(species_df$acronym))
    pb <- txtProgressBar(min = 1, max = ntot, style = 3)
    for(i in 1:ntot) {
      if (is.na(species_df$acronym[i])){
        species_df$acronym[i] <- acronymDB$acronym[stringdist::amatch(rownames(species_df[i,]), trimws(acronymDB$species), maxDist=maxDistTaxa, matchNA = F)]
        setTxtProgressBar(pb, i)
      }
    }
    close(pb)
    print(paste("Acronyms were found for ", length(nrow(species_df)) - sum(is.na(species_df$acronym)) , " species out of ", length(species_df$acronym), sep=""))
  }
  #taxa unrecognized by acronym
  print("Exporting list of species with the acronyms found. NA= not found")
  write.csv(cbind(rownames(species_df), species_df$acronym), paste(resultsPath,"\\Recognized_acronyms.csv", sep=""))

  #removes NA
  species_df[is.na(species_df)] <- 0

  ########## LINK WITH DIAT.BARCODE DATABASE
  #internal Diat.Barcode number
  intversion <- "8.1"
  #get version number of latest Diat.Barcode
  dic <- read.csv("http://www.francoiskeck.fr/work/diatbarcode/dic_version.csv", header = TRUE, stringsAsFactors = FALSE)
  if (exists("dic")){
    #is able to check the version
    version <- dic$Version[which.max(as.numeric(as.POSIXlt(dic$Date, format = "%d-%m-%Y")))]
    #compare both. If updates are needed, attempt them
    if (version == intversion){
      #updates are not needed
      print ("No updates needed for the Diat.barcode database. Proceeding")
      #load("data/dbc_offline.RData") ##takes the internal database
      #dbc <- dbc_offline
      dbc <- DiaThor:::dbc_offline

    } else {
      #updates are needed
      ########--------  Diat.Barcode download attempt. If it fails, tries to use internal database
      print("Attempting to download diat.barcode from website")
      dbc <- diatbarcode::get_diatbarcode(version = "last") #loads the latest version of diat.barcode
      if (exists("dbc")){ #it if was able to download the new version, proceed
        print("Latest version of Diat.barcode succesfully downloaded. Remember to credit accordingly!")
      } else { #it if was unable to download the new version, use the internal database
        print("Latest version of Diat.barcode cannot be downloaded")
        print("Using internal database, Diat.barcode v.8.1 published on 10-06-2020. It might need to be updated")
        #load("data/dbc_offline.RData") ##takes the internal database
        #dbc <- dbc_offline
        dbc <- DiaThor:::dbc_offline
      }

    }
  } else {
    print("Latest version of Diat.barcode unknown")
    print("Using internal database, Diat.barcode v.8.1 published on 10-06-2020. It might need to be updated")
  }
  ### Double checks that database got loaded correctly or cancels alltogether
  if (exists("dbc")){
    print("Proceeding with calculations")
  } else { #if everything fails, cancels
    stop("Database could not be downloaded or retrieved from internal package. Cancelling calculations")
  }

  ########## END LINK WITH DIAT.BARCODE DATABASE

  #Remove duplicate by field "species"
  dbc2 <- as.data.frame(dbc[!duplicated(dbc[,"species"]),]) #transforms dbc to a dataframe

  ecodata <- dbc2[which(colnames(dbc2)=="species" ):ncol(dbc2)] #keeps only from the "species" column onwards, to keep the ecological data

  ########-------- MATCHES AND BINDS DATA SETS
  # NOTE: if several species match with the same name in the diat.barcode database (multiple strains), it keeps the
  # first occurence. Ecological indices usually match in all strains of the same species

  #fuzzy matching
  taxaInSp <- species_df[stringdist::ain(row.names(species_df), ecodata[,"species"], maxDist=maxDistTaxa, matchNA = F),]

  numberOfSamples <- ncol(taxaInSp) #Saves the number of columns that correspond to samples
  sampleNames <- colnames(taxaInSp) #Saves the names of the samples

  taxaInEco <- ecodata[stringdist::ain(ecodata[,"species"], row.names(species_df), maxDist=maxDistTaxa),] #FUZZY matching species in table of indices
  taxaInEco <- taxaInEco[!duplicated(taxaInEco[,"species"]),] #we have to make taxaInEco same length as taxaInSp removing duplicate species
  taxaInEco <- taxaInEco[order(taxaInEco$species),] #reorder alphabetically

  taxaIn <- cbind(taxaInSp, taxaInEco) #and we bind taxaInSp and taxaInEco
  #remove the acronym col and move it to the back
  acronym_col <- taxaIn$acronym
  acronym_col <- replace(acronym_col, acronym_col=="0", NA)
  taxaIn<- taxaIn[ , !(names(taxaIn) == "acronym")] #removes the acronym names column
  taxaIn$acronym <- acronym_col

  write.csv(taxaIn, paste(resultsPath,"\\Taxa included.csv", sep=""))
  print(paste("Number of taxa included:", nrow(taxaIn), "-- Detailed list in 'Taxa included.csv'"))

  #also makes a matrix for all the taxa left out, for the user to review
  '%out%' <- function(x,y)!('%in%'(x,y))
  taxaOut <- species_df[row.names(species_df)%out%row.names(taxaInSp),]
  print(paste("Number of taxa excluded:", nrow(taxaOut), "-- Detailed list in 'Taxa excluded.csv'" ))
  write.csv(taxaOut, paste(resultsPath,"\\Taxa excluded.csv", sep=""))

  #gets the column named "species" in the diat.barcode database, everything before that column should be a sample with abundance data
  lastcol = which(colnames(taxaIn)=="species")

  #gets sample names
  sampleNames <- colnames(taxaIn[1:(lastcol-1)])
  #CONVERT TO RELATIVE ABUNDANCE
  #Convert taxaIn sample data to Relative Abundance data
  if(isRelAb==FALSE){ #pass parameter when in function
    taxaInRA <- taxaIn
    print("Converting species' densities to relative abundance")
    for (i in 1:nrow(taxaIn)){
      for (j in 1:(lastcol-1)){
        if (is.na(taxaIn[i,j])){
          taxaInRA[i,j] <- 0
        } else {
          taxaInRA[i,j] <- (taxaIn[i,j]*100)/sum(taxaIn[,j])
        }

      }
    }
  }
  resultList <- list(as.data.frame(taxaInRA), as.data.frame(taxaIn), sampleNames, resultsPath)
  names(resultList) <- c("taxaInRA", "taxaIn", "sampleNames", "resultsPath")
  return(resultList)

}

