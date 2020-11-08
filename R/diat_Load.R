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

diat_findAcronyms <- function(species_df, maxDistTaxa=2, resultsPath){
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

  #Results fomder
  if(missing(resultsPath)) {
    print("Select Results folder")
    resultsPath <- choose.dir(default = "", caption = "Select folder for your Results")
  }
  if (is.na(resultsPath)){stop("Calculations cancelled, no folder selected")}

  acronymDB <- as.data.frame(DiaThor::acronyms) #loads taxonomy
  #acronymDB <- read.csv("Indices/acronyms.csv") #uses the external csv file


  #trims acronyms and converts them to characters, so it doesnt recognize them as factors (darn R)
  acronymDB$species <- as.character(acronymDB$species)
  acronymDB$species <- trimws(acronymDB$species)
  acronymDB$acronym <- as.character(acronymDB$acronym)
  acronymDB$acronym <- trimws(acronymDB$acronym)
  acronymDB$new_acronym <- as.character(acronymDB$new_acronym)
  acronymDB$new_acronym <- trimws(acronymDB$new_acronym)


  #FINDS ACRONYMS FROM SPECIES OR PROCEEDS TO ACRONYM UPDATE
  if("acronym" %in% colnames(species_df)) {
    #there is an acronym, updates
    print("Acronyms column found in file. Updating acronyms")
    species_df$acronym <- as.character(species_df$acronym)

  }else{
    #there is no acronym column, creates one
    species_df$acronym <- NA

    print ("Searching for acronyms in exact match, please wait...")
    #exact match by species name (fast)
    species_df$acronym <- acronymDB$acronym[match(trimws(rownames(species_df)), acronymDB$species)]
    #species_df_pre <- species_df
    #species_df <- species_df_pre
    #with the unmatched exact, try heuristic
    print ("Searching for acronyms in heuristic match, please wait...")
    #PROGRESS BAR FOR HEURISTIC SEARCH OF ACRONYMS
    ntot <- nrow(species_df)
    pb <- txtProgressBar(min = 1, max = ntot, style = 3)
    for(i in 1:ntot) {
      if (is.na(species_df$acronym[i])){
        species_df$acronym[i] <- acronymDB$acronym[stringdist::amatch(trimws(rownames(species_df[i,])), trimws(acronymDB$species), maxDist=maxDistTaxa, matchNA = F)]
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  #print(paste("Updated acronyms: ", updatedacronyms))
  print(paste("Acronyms were found for ", nrow(species_df) - sum(is.na(species_df$acronym)) , " species out of ", length(species_df$acronym), sep=""))

  #original species and acronym list
  speciesList <- trimws(rownames(species_df))
  acronymList <- trimws(species_df$acronym)
  #builds column with new name for species
    species_df$new_species <- NA

    #UPDATE ACRONYMS
  updatedacronyms <- 0
  species_df$acronym <- as.character(species_df$acronym)

  #keepupdating <- T
  #while (keepupdating == T){
    updatedacronyms <- 0
    for(i in 1:nrow(species_df)) {
      if (!is.na(species_df$acronym[i])){ #there is an acronym for this species, try updating
        updatedacronym <- acronymDB$new_acronym[match(species_df$acronym[i], acronymDB$acronym)] #get updated acronym
        updatedspecies <- acronymDB$new_species[match(species_df$acronym[i], acronymDB$acronym)] #get updated species
        if (is.na(updatedacronym) | updatedacronym == ""){
          #print(paste("i=", i, " already exists"))

        } else{
          #species_df$new_acronym[i] <- updatedacronym #update the acronym
          species_df$new_species[i] <- updatedspecies #update the species
          updatedacronyms <- updatedacronyms + 1
        }
      }
    }
    if (updatedacronyms > 0){
      print(paste(updatedacronyms, "updated acronyms were found and exported"))

    }

    #if (updatedacronyms==0){keepupdating <- F}
  #}



  #taxa recognized by acronym
  print("Exporting list of species with the acronyms found by heuristic search. NA= not found")
  outacronyms <- cbind(speciesList, species_df$new_species, acronymList, species_df$acronym)
  # outacronyms <- cbind(speciesList, acronymList, species_df$acronym)

  colnames(outacronyms) <- c("Original species","Updated species" ,"Original Acronym", "Updated acronym")
  write.csv(outacronyms, paste(resultsPath,"\\Recognized_acronyms.csv", sep=""))
  return(species_df)

}

########-------- FUNCTION TO LOAD FILE  --------------#########
#### IN THIS SECTION WE READ THE MATRIX
### INPUT: CSV file with the samples * abundances or relative abundances
### OUTPUTS: a dataframe with the species as matched against the database, with species in RA; Taxaincluded and Taxaexcluded in CSV in the Output
### folder, detailing which taxa were recognized and which were not

diat_loadData <- function(species_df, isRelAb=FALSE, maxDistTaxa=2, resultsPath){
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(species_df)) {
    print("Select CSV matrices with your sample data")
    Filters <- matrix(c("Comma Separated Values (CSV)", "*.csv"),
                      1, 2, byrow = TRUE)
    species_df <- as.data.frame(read.csv(file.choose())) #select the data matrix with a dialog box
  }

  #Results folder
  if(missing(resultsPath)) {
    print("Select Results folder")
    resultsPath <- choose.dir(default = "", caption = "Select folder for your Results")
  }
  if (is.na(resultsPath)){stop("Calculations cancelled, no folder selected")}
  print("Result folder selected")

  #check for duplicate species
  if("species" %in% colnames(species_df)) {
    print("Checking for duplicates")
    testdupl <- species_df[,"species"]
    if (length(testdupl[duplicated(testdupl)])> 0){ #check for duplicates
      print("Duplicate species found in input data")
      testdupl[duplicated(testdupl)] #show duplicates
      stop("Cancelling.Duplicate species found in input data") #abort script
    }
    speciesrow <- species_df[,"species"]
    row.names(species_df) <- speciesrow  #converts the species names to rownames
    species_df<- species_df[ , !(names(species_df) == "species")] #removes the species names column
  } else {
    #check if species names are in the first column
    stop("Calculations cancelled. File should contain a column named 'species' with the species list")
  }
  #Reorder alphabetically
  #species_df <- species_df[order(row.names(species_df)),] #reorder alphabetically


  # Find acronyms or update then
  print("Finding and updating acronyms")
  species_df <- diat_findAcronyms(species_df, maxDistTaxa, resultsPath)

  #removes NA in abundance data
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
    dbc <- DiaThor:::dbc_offline
  }
  ### Double checks that database got loaded correctly or cancels alltogether
  if (exists("dbc")){
  } else { #if everything fails, cancels
    stop("Database could not be downloaded or retrieved from internal package. Cancelling calculations")
  }
  ########## END LINK WITH DIAT.BARCODE DATABASE

  #Remove duplicate by field "species" in diat.barcode
  dbc2 <- as.data.frame(dbc[!duplicated(dbc[,"species"]),]) #transforms dbc to a dataframe
  ecodata <- dbc2[which(colnames(dbc2)=="species" ):ncol(dbc2)] #keeps only from the "species" column onwards, to keep the ecological data

  ########-------- MATCHES AND BINDS DATA SETS
  # NOTE: if several species match with the same name in the diat.barcode database (multiple strains), it keeps the
  # first occurence. Ecological indices usually match in all strains of the same species
  #fuzzy matching

  #NEW SEARCH
  rownames(species_df) <- trimws(rownames(species_df))
  taxaInSp <- as.data.frame(matrix(nrow=nrow(species_df), ncol = (ncol(species_df)+ncol(ecodata))))
  #species_df section
  taxaInSp[1:nrow(species_df),1:ncol(species_df)] <- species_df[1:nrow(species_df),1:ncol(species_df)] #copies species_df into taxaInSp
  colnames(taxaInSp)[1:ncol(species_df)] <- colnames(species_df) #copies column names
  rownames(taxaInSp) <- rownames(species_df) #copies row names
  taxaInSp$recognizedSp <- NA #creates a new column with the recognized species

  #ecodata_section
  lastcolspecies_df <-  which(colnames(taxaInSp)=="new_species") #gets last column of taxaInSp with species_df data
  colnames(taxaInSp)[(lastcolspecies_df+1):(ncol(taxaInSp)-1)] <- colnames(ecodata) #copies column names
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = nrow(taxaInSp), style = 3)
  for (i in 1:nrow(taxaInSp)){
    searchvectr <- ecodata[stringdist::ain(ecodata[,"species"],row.names(species_df)[i], maxDist=maxDistTaxa, matchNA = F),] #seaches species by species
    if (nrow(searchvectr)==1){ #if it finds only one species, add that
      taxaInSp[i,(lastcolspecies_df + 1):(ncol(taxaInSp)-1)] <- searchvectr
      taxaInSp[i,"recognizedSp"] <- searchvectr$species
    } else if (nrow(searchvectr)>1){ #if it finds multiple species, keeps the one with the lower distance
      searchvectr <- searchvectr[which(stringdist::stringdist(searchvectr$species, row.names(species_df)[i]) == min(stringdist::stringdist(searchvectr$species, row.names(species_df)[i]))),]
      if (nrow(searchvectr) > 1) { #still finds more than one with the same lower distance, keeps the first
        searchvectr <- searchvectr[1,]
      }
      taxaInSp[i,(lastcolspecies_df + 1):(ncol(taxaInSp)-1)] <- searchvectr
      taxaInSp[i,"recognizedSp"] <- searchvectr$species

    } else if (nrow(searchvectr)==0){ #species not found at all
      taxaInSp[i,"recognizedSp"] <- "Not found"
    }
    #update progressbar
    setTxtProgressBar(pb, i)
  }
  #close progressbar
  close(pb)
  #END NEW SEARCH
  taxaInEco <- taxaInSp
  taxaIn <- species_df #dataframe to be exported for indices
  taxaIncluded <- as.data.frame(rownames(taxaInEco)[which(taxaInEco$recognizedSp != "Not found")])
  taxaExcluded <- as.data.frame(rownames(taxaInEco)[which(taxaInEco$recognizedSp == "Not found")])

  #remove the acronym col and move it to the back
  acronym_col <- taxaIn$acronym
  acronym_col2 <- taxaInEco$acronym
  acronym_col <- replace(acronym_col, acronym_col=="0", NA)
  acronym_col2 <- replace(acronym_col2, acronym_col2=="0", NA)
  taxaIn<- taxaIn[ , !(names(taxaIn) == "acronym")] #removes the acronym names column
  taxaInEco<- taxaInEco[ , !(names(taxaInEco) == "acronym")] #removes the acronym names column
  taxaInEco$acronym <- acronym_col2
  taxaIn$acronym <- acronym_col

  #remove the updated species col and move it to the back in both taxaIn
  newspecies_col <- taxaIn$new_species
  newspecies_col <- replace(newspecies_col, newspecies_col=="0", NA)
  taxaIn<- taxaIn[ , !(names(taxaIn) == "new_species")] #removes the newspecies_col
  taxaIn$new_species <- newspecies_col

  newspecies_col2 <- taxaInEco$new_species
  newspecies_col2 <- replace(newspecies_col2, newspecies_col2=="0", NA)
  taxaInEco<- taxaInEco[ , !(names(taxaInEco) == "new_species")] #removes the newspecies_col
  taxaInEco$new_species <- newspecies_col2

  sampleNames <- colnames(taxaIn) #Saves the names of the samples
  removeelem <- c("acronym","new_species","species") #Removes columns not samples
  sampleNames <- sampleNames[!(sampleNames %in% removeelem)]


  if (nrow(taxaIncluded) == 0) {
    print("No taxa were recognized for morphology analysis. Check taxonomy instructions for the package")
  }
  #Exports included taxa in morphology analyses
  colnames(taxaIncluded) <- "Eco/Morpho"
  write.csv(taxaIncluded, paste(resultsPath,"\\Taxa included.csv", sep=""))
  print(paste("Number of taxa recognized for morphology:", nrow(taxaIncluded), "-- Detailed list in 'Taxa included.csv'"))

  #also makes a matrix for all the taxa left out, for the user to review
  colnames(taxaExcluded) <- "Eco/Morpho"
  write.csv(taxaExcluded, paste(resultsPath,"\\Taxa excluded.csv", sep=""))
  print(paste("Number of taxa excluded for morphology:", nrow(taxaExcluded), "-- Detailed list in 'Taxa excluded.csv'" ))

  #Has to clean the TaxaInEco for those species that were not found
  taxaInEco <- taxaInEco[which(taxaInEco$recognizedSp != "Not found"),]

  #creates a blank precision matrix for the indices
  precisionmatrix <- as.data.frame(sampleNames)
  names(precisionmatrix)[names(precisionmatrix)=="sampleNames"] <- "Sample"

  write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))

  #gets the column named "acronym" everything before that column should be a sample with abundance data
  lastcol = which(colnames(taxaIn)=="acronym")
  #
  # #gets sample names
  # sampleNames <- colnames(taxaIn[1:(lastcol-1)])

  #CREATES A RELATIVE ABUNDANCE MATRIX AS WELL FOR THOSE INDICES THAT USE IT
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
  } else {
    taxaInRA <- taxaIn
  }

  #CREATES THE EXPORT PRODUCTS
  resultList <- list(as.data.frame(taxaInRA), as.data.frame(taxaIn), sampleNames, resultsPath, taxaInEco)
  names(resultList) <- c("taxaInRA", "taxaIn", "sampleNames", "resultsPath", "taxaInEco")
  print("Data loaded")
  return(resultList)

}
