#' These functions calculate morphological and ecological parameters from the input data
#' @param resultLoad The result dataframe obtained from the loadFiles function
#' @param vandamReports Boolean. If set to 'TRUE' the detailed reports for the Van Dam classifications will be reported in the Output. Default = TRUE
#' @description
#' The input for these functions is the resulting dataframe obtained from the loadFiles function to calculate morphological parameters, diversity values (using the vegan package), ecological guilds and the Van Dam ecological preferences
#' The morphological data (size classes, chlorophlasts) is obtained from the Diat.Barcode project. Besides citing this package, the Diat.Barcode project should also be cited if the package is used, as follows:
#' \itemize{
#' \item Rimet, Frederic; Gusev, Evgenuy; Kahlert, Maria; Kelly, Martyn; Kulikovskiy, Maxim; Maltsev, Yevhen; Mann, David; Pfannkuchen, Martin; Trobajo, Rosa; Vasselon, Valentin; Zimmermann, Jonas; Bouchez, Agnès, 2018, "Diat.barcode, an open-access barcode library for diatoms", https://doi.org/10.15454/TOMBYZ, Portail Data Inra, V1
#' }
#' Guild classification is obtained from:
#' \itemize{
#' \item CITA GUILDS
#' }
#' Van Dam classification is obtained form:
#' \itemize{
#' \item CITA VAN DAM
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
#' @export diat_morpho
#' @export diat_size
#' @export diat_diversity
#' @export diat_guild
#' @export diat_vandam
#'


###### ---------- FUNCTION FOR MORPHOLOGICAL DATA: DONE ---------- ########
#### IN THIS SECTION WE CALCULATE MORPHOLOGICAL VARIABLES
### INPUT: result created in loadData(). If its in RA, biovolume cannot be calculated
### OUTPUTS: #number of chloroplasts #shape of chloroplasts #biovolume
diat_morpho <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]
  taxaInRA <- resultLoad[[1]]

  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaInRA)=="species")
  #gets sample names
  sampleNames <- colnames(taxaInRA[1:(lastcol-1)])


  ##### Number of Chloroplasts
  #creates an empty dataframe
  numcloroplastos <- c(unique(taxaInRA[,"chloroplast_number"]))
  numcloroplastos.result <- data.frame(matrix(ncol = length(numcloroplastos), nrow = (lastcol-1)))
  colnames(numcloroplastos.result) <- paste("chloroplasts - ", (numcloroplastos), sep ="")

  #PROGRESS BAR
  print("Calculating chloroplast quantity")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix

    numcloroplastos.ab <- NULL
    #sum abundances per category
    for (i in numcloroplastos){
      abc <- sum(taxaInRA[which(taxaInRA$chloroplast_number == i),sampleNumber])
      #round numbers and remove negatives
      if (abc<0){abc <- 0}
      abc <- round(abc, digits=3)
      numcloroplastos.ab <- c(numcloroplastos.ab, abc)
    }

    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    numcloroplastos.result[sampleNumber, ] <- numcloroplastos.ab
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  ##### END Number of Chloroplasts

  ##### Shape of Chloroplasts
  #creates results dataframe
  shpcloroplastos <- c(unique(taxaInRA[,"chloroplast_shape"]))
  shpcloroplastos.result <- data.frame(matrix(ncol = length(shpcloroplastos), nrow = (lastcol-1)))
  colnames(shpcloroplastos.result) <- paste("shape chloroplasts -", (shpcloroplastos), sep ="")

  #PROGRESS BAR
  print("Calculating chloroplast shapes")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    shpcloroplastos <- c(unique(taxaInRA[,"chloroplast_shape"]))
    shpcloroplastos.ab <- NULL
    #sum abundances per category
    for (i in shpcloroplastos){
      abc <- sum(taxaInRA[which(taxaInRA$chloroplast_shape == i),sampleNumber])
      #round numbers and remove negatives
      if (abc<0){abc <- 0}
      abc <- round(abc, digits=3)
      shpcloroplastos.ab <- c(shpcloroplastos.ab, abc)
    }

    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    shpcloroplastos.result[sampleNumber, ] <- shpcloroplastos.ab
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  ##### END Shape of Chloroplasts

  ##### Biovolume
  #creates result dataframe
  biovol.val.result <- data.frame(matrix(ncol = 1, nrow = (lastcol-1)))
  colnames(biovol.val.result) <- "Total Biovolume"
  #PROGRESS BAR
  print("Calculating biovolume")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)

  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    biovol <- c(unique(taxaIn[,"biovolume"]))
    biovol.val <- NULL
    #sum abundances*biovolume
    biovol.val <- sum((taxaIn[which(biovol != 0),sampleNumber])*taxaIn[which(biovol != 0),"biovolume"])
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    biovol.val.result[sampleNumber, ] <- biovol.val
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)

  ##### END Biovolume
  morphoresultTable <- NULL
  morphoresultTable <- data.frame(c( if(exists("numcloroplastos.result")){numcloroplastos.result}, if(exists("shpcloroplastos.result")){shpcloroplastos.result}, if(exists("biovol.val.result")){biovol.val.result}))
  rownames(morphoresultTable) <- sampleNames

  morphoresultTable <- list(as.data.frame(numcloroplastos.result), as.data.frame(shpcloroplastos.result), as.data.frame(if(exists("biovol.val.result")){biovol.val.result}))
  names(morphoresultTable) <- c("numcloroplastos.result", "shpcloroplastos.result", "biovol.val.result")
  return(morphoresultTable)
}



###### ---------- FUNCTION FOR SIZE CLASSES: DONE  ---------- ########
#### IN THIS SECTION WE CALCULATE SIZE CLASSES ACCORDING TO....
### INPUT: resultLoad created in loadData()
### OUTPUTS: dataframe with % of size classes

diat_size <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]
  taxaInRA <- resultLoad[[1]]

  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaInRA)=="species")
  #gets sample names
  sampleNames <- colnames(taxaInRA[1:(lastcol-1)])

  ##### SIZE CLASSES
  #creates results dataframe
  size_labels <- c("Size class 1", "Size class 2", "Size class 3", "Size class 4", "Size class 5", "Size class Indet", "Size Taxa used")
  size.results <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(size.results) <- size_labels

  #finds the size class column
  size_v <-  taxaInRA[,startsWith(colnames(taxaInRA), "classe_de_taille")]
  #remove the NA
  size_v[is.na(size_v)] = 0

  #PROGRESS BAR
  print("Calculating size classes")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix

    #sum abundances per category
    size_cat1 <- sum(taxaInRA[which(size_v == 1),sampleNumber])
    size_cat2 <- sum(taxaInRA[which(size_v == 2),sampleNumber])
    size_cat3 <- sum(taxaInRA[which(size_v == 3),sampleNumber])
    size_cat4 <- sum(taxaInRA[which(size_v == 4),sampleNumber])
    size_cat5 <- sum(taxaInRA[which(size_v == 5),sampleNumber])
    size_indet <- 100 - sum(size_cat1, size_cat2, size_cat3, size_cat4, size_cat5) #calculates indetermined
    #round numbers and remove negatives
    if (size_indet<0){size_indet <- 0}
    size_cat1 <- round(size_cat1, digits=3)
    size_cat2 <- round(size_cat2, digits=3)
    size_cat3 <- round(size_cat3, digits=3)
    size_cat4 <- round(size_cat4, digits=3)
    size_cat5 <- round(size_cat5, digits=3)
    size_indet <- round(size_indet, digits=3)

    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    sizetaxaused <- length(which(size_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    sizetaxaused_taxa <- taxaInRA[which(size_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    size_values <- c(size_cat1, size_cat2, size_cat3, size_cat4, size_cat5, size_indet, sizetaxaused)
    size.results[sampleNumber, ] <- size_values
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)

  }
  #close progressbar
  close(pb)
  ##### END SIZE CLASSES
  return(size.results)
}


###### ---------- FUNCTION FOR ECOLOGICAL INDICES: DONE   ---------- ########
#### IN THIS SECTION WE CALCULATE BASIC ECOLOGICAL INDICES WITH VEGAN
### INPUT: resultLoad Cannot be in Relative Abundance
### OUTPUTS: dataframe with diversity indices per sample

diat_diversity <- function(resultLoad){


  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]
  taxaInRA <- resultLoad[[1]]

  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaInRA)=="species")

  #remove NA to 0
  taxaIn[is.na(taxaIn)] <- 0

  sampleNames <- colnames(taxaIn[1:(lastcol-1)])
  diversityIndices <- data.frame(matrix(ncol = (lastcol-1), nrow = 3))

  #PROGRESS BAR
  print("Calculating diversity indices")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  #samples are columns

  for (i in 1:(lastcol-1)){

    richness <-vegan::specnumber(taxaIn[,i])
    shannon <- vegan::diversity(taxaIn[,i])
    evenness <- shannon/log(vegan::specnumber(taxaIn[,i]))
    diversityIndices[1,i] <- richness
    diversityIndices[2,i] <- shannon
    diversityIndices[3,i] <- evenness
    #update progressbar
    setTxtProgressBar(pb, i)
  }
  #close progressbar
  close(pb)
  #RESULTS
  diversity.results<-as.data.frame(t(diversityIndices)) #transposes the diversity matrix to plot
  rownames(diversity.results) <- colnames(taxaIn[1:(lastcol-1)])
  colnames(diversity.results) <- c("Richness", "Shannon H", "Evenness")

  return(diversity.results)
}


###### ---------- FUNCTION FOR ECOLOGICAL GUILDS: DONE   ---------- ########
#### IN THIS SECTION WE CALCULATE ECOLOGICAL GUILDS ACCORDING TO PASSY (2017)
### INPUT: taxaInRA created in loadData()
### OUTPUTS: dataframe with ecological guilds  per sample

diat_guild <- function(resultLoad){
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]
  taxaInRA <- resultLoad[[1]]

  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaInRA)=="species")

  #creates results dataframe
  guild_labels <- c("Guild: HP", "Guild: LP", "Guild: Mot", "Guild: Plank", "Guild: Indet", "Guilds Taxa used")
  guilds.results <- data.frame(matrix(ncol = 6, nrow = (lastcol-1)))
  colnames(guilds.results) <- guild_labels

  #PROGRESS BAR
  print("Calculating ecological guilds")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix

    guild_HP <- taxaInRA[,startsWith(colnames(taxaInRA), "high_profile_guild")]
    guild_LP <- taxaInRA[,startsWith(colnames(taxaInRA), "low_profile_guild")]
    guild_Mot <- taxaInRA[,startsWith(colnames(taxaInRA), "motile_guild")]
    guild_Plank <- taxaInRA[,startsWith(colnames(taxaInRA), "euplanctonic_guild")]
    #remove the NA
    guild_HP[is.na(guild_HP)] = 0
    guild_LP[is.na(guild_LP)] = 0
    guild_Mot[is.na(guild_Mot)] = 0
    guild_Plank[is.na(guild_Plank)] = 0

    #sum abundances per guild
    guild_HP_ab <- sum(taxaInRA[which(guild_HP == 1),sampleNumber])
    guild_LP_ab <- sum(taxaInRA[which(guild_LP == 1),sampleNumber])
    guild_Mot_ab <- sum(taxaInRA[which(guild_Mot == 1),sampleNumber])
    guild_Plank_ab <- sum(taxaInRA[which(guild_Plank == 1),sampleNumber])
    guild_indet <- 100 - sum(guild_HP_ab, guild_LP_ab, guild_Mot_ab, guild_Plank_ab) #calculates indetermined


    #round numbers and remove negatives
    if (guild_indet<0){guild_indet <- 0}
    guild_HP_ab <- round(guild_HP_ab, digits=3)
    guild_LP_ab <- round(guild_LP_ab, digits=3)
    guild_Mot_ab <- round(guild_Mot_ab, digits=3)
    guild_Plank_ab <- round(guild_Plank_ab, digits=3)
    guild_indet <- round(guild_indet, digits=3)

    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    guildtaxaused <- length(which(guild_HP == 1 & taxaInRA[,sampleNumber] > 0))
    guildtaxaused <- guildtaxaused + length(which(guild_LP == 1 & taxaInRA[,sampleNumber] > 0))
    guildtaxaused <- guildtaxaused + length(which(guild_Mot == 1 & taxaInRA[,sampleNumber] > 0))
    guildtaxaused <- guildtaxaused + length(which(guild_Plank == 1 & taxaInRA[,sampleNumber] > 0))

    #which taxa were used? to export
    guildtaxaused_taxa <- taxaInRA[which(guild_HP == 1 & taxaInRA[,sampleNumber] > 0),"species"]
    guildtaxaused_taxa <- c(guildtaxaused_taxa, taxaInRA[which(guild_LP == 1 & taxaInRA[,sampleNumber] > 0),"species"])
    guildtaxaused_taxa <- c(guildtaxaused_taxa, taxaInRA[which(guild_Mot == 1 & taxaInRA[,sampleNumber] > 0),"species"])
    guildtaxaused_taxa <- c(guildtaxaused_taxa, taxaInRA[which(guild_Plank == 1 & taxaInRA[,sampleNumber] > 0),"species"])

    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    guild_values <- c(guild_HP_ab, guild_LP_ab, guild_Mot_ab, guild_Plank_ab, guild_indet, guildtaxaused)
    guilds.results[sampleNumber, ] <- guild_values
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)

  ##### END GUILDS
  return(guilds.results)

}


###### ---------- FUNCTION FOR VANDAM ECOLOGICAL CHARACTERIZATION: DONE   ---------- ########
#### IN THIS SECTION WE CALCULATE ECOLOGICAL PREFERENCES ACCORDING TO VAN DAM (REF)
### INPUT: resultLoad created in loadData()
### OUTPUTS: dataframe with VanDam's ecological values
### PARAMETERS: vandamReports (Boolean), if a detailed report of the taxa used for each sub index is exported or not
#vdams = salinity
#vdamnh = N-Heterotrophie
#vdamo2 = Oxygen
#vdamsap = saprobity
#vdamtrop = trophic
#vdamaero = aerophilie

diat_vandam <- function(resultLoad, vandamReports=T){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]
  taxaInRA <- resultLoad[[1]]

  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaInRA)=="species")

  resultsPath <- resultLoad[[4]]
  # CHECKS IF RESULTSPATH EXISTS, OR ASKS FOR IT, FOR THE DETAILED REPORTS
  if (vandamReports==T & is.na(resultsPath)==T){
    print("Select Results folder for the detailed VanDam reports")
    resultsPath <- choose.dir(default = "", caption = "Select folder for your Results")
  }

  #EXPORT REPORTS FOR EACH SAMPLE OF TAXA USED
  if (vandamReports==T & is.na(resultsPath)==F){
    print("Exporting detailed reports for VanDam ecological preferences")
    print(resultsPath)
    vandamtaxafile = paste("VanDam Taxa used.txt", sep="")
    write("TAXA USED FOR EACH ECOLOGICAL VARIABLE USING VANDAM's CLASSIFICATION", paste(resultsPath, "\\", vandamtaxafile, sep=""))
    write("These taxa were included because: ",paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)
    write("a) they had a reliable classification in VanDam's classification system",paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)
    write("b) they had a relative abundance in the sample > 0",paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)
  }

  ##### VANDAM SALINITY
  #creates results dataframe
  vdams_labels <- c("VD Salinity 1", "VD Salinity 2", "VD Salinity 3", "VD Salinity 4", "VD Salinity Indet", "VD Salinity Taxa used")
  vdamSalinity <- data.frame(matrix(ncol = 6, nrow = (lastcol-1)))
  colnames(vdamSalinity) <- vdams_labels

  #finds the column
  vdams_v <- (taxaInRA[,"vdams"])
  #remove the NA
  vdams_v[is.na(vdams_v)] = 0
  print("Calculating Van Dam salinity")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdams_cat1 <- sum(taxaInRA[which(vdams_v == 1),sampleNumber])
    vdams_cat2 <- sum(taxaInRA[which(vdams_v == 2),sampleNumber])
    vdams_cat3 <- sum(taxaInRA[which(vdams_v == 3),sampleNumber])
    vdams_cat4 <- sum(taxaInRA[which(vdams_v == 4),sampleNumber])
    vdams_indet <- 100 - sum(vdams_cat1, vdams_cat2, vdams_cat3, vdams_cat4) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    Vdamstaxaused <- length(which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamstaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    vdams_values <- c(vdams_cat1, vdams_cat2, vdams_cat3, vdams_cat4, vdams_indet, Vdamstaxaused)
    vdamSalinity[sampleNumber, ] <- vdams_values
    if (vandamReports==T & exists("resultsPath")) {write(paste("SALINITY- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
    if (vandamReports==T & exists("resultsPath")) {write(Vdamstaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
  }

  ##### END VANDAM SALINITY

  ##### VANDAM N-HETEROTROPHIE
  #creates results dataframe
  vdamnh_labels <- c("VD N-Het 1", "VD N-Het 2", "VD N-Het 3", "VD N-Het 4", "VD N-Het Indet", "VD N-Het Taxa used")
  vdamNHeterotrophy <- data.frame(matrix(ncol = 6, nrow = (lastcol-1)))
  colnames(vdamNHeterotrophy) <- vdamnh_labels
  #finds the column
  vdamnh_v <- (taxaInRA[,"vdamnh"])
  #remove the NA
  vdamnh_v[is.na(vdamnh_v)] = 0
  print("Calculating Van Dam N-heterotrophy")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamnh_cat1 <- sum(taxaInRA[which(vdamnh_v == 1),sampleNumber])
    vdamnh_cat2 <- sum(taxaInRA[which(vdamnh_v == 2),sampleNumber])
    vdamnh_cat3 <- sum(taxaInRA[which(vdamnh_v == 3),sampleNumber])
    vdamnh_cat4 <- sum(taxaInRA[which(vdamnh_v == 4),sampleNumber])
    vdamnh_indet <- 100 - sum(vdamnh_cat1, vdamnh_cat2, vdamnh_cat3, vdamnh_cat4) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    Vdamnhtaxaused <- length(which(vdamnh_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamnhtaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    vdamnh_values <- c(vdamnh_cat1, vdamnh_cat2, vdamnh_cat3, vdamnh_cat4, vdamnh_indet, Vdamnhtaxaused)
    vdamNHeterotrophy[sampleNumber, ] <- vdamnh_values
    if (vandamReports==T & exists("resultsPath")) {write(paste("N-HETEROTROPHY- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
    if (vandamReports==T & exists("resultsPath")) {write(Vdamnhtaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
  }
  ##### END VANDAM N-HETEROTROPHIE

  ##### VANDAM OXYGEN REQUIREMENTS
  #creates results dataframe
  vdamo2_labels <- c("VD Oxygen 1", "VD Oxygen 2", "VD Oxygen 3", "VD Oxygen 4","VD Oxygen 5", "VD Oxygen Indet", "VD Oxygen Taxa used")
  vdamOxygen <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(vdamOxygen) <- vdamo2_labels
  #finds the column
  vdamo2_v <- (taxaInRA[,"vdamo2"])
  #remove the NA
  vdamo2_v[is.na(vdamo2_v)] = 0
  print("Calculating Van Dam oxygen requirements")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamo2_cat1 <- sum(taxaInRA[which(vdamo2_v == 1),sampleNumber])
    vdamo2_cat2 <- sum(taxaInRA[which(vdamo2_v == 2),sampleNumber])
    vdamo2_cat3 <- sum(taxaInRA[which(vdamo2_v == 3),sampleNumber])
    vdamo2_cat4 <- sum(taxaInRA[which(vdamo2_v == 4),sampleNumber])
    vdamo2_cat5 <- sum(taxaInRA[which(vdamo2_v == 5),sampleNumber])
    vdamo2_indet <- 100 - sum(vdamo2_cat1, vdamo2_cat2, vdamo2_cat3, vdamo2_cat4, vdamo2_cat5) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    Vdamo2taxaused <- length(which(vdamo2_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamo2taxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    vdamo2_values <- c(vdamo2_cat1, vdamo2_cat2, vdamo2_cat3, vdamo2_cat4, vdamo2_cat5, vdamo2_indet, Vdamo2taxaused)
    vdamOxygen[sampleNumber, ] <- vdamo2_values
    if (vandamReports==T & exists("resultsPath")) {write(paste("OXYGEN REQUIREMENTS- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
    if (vandamReports==T & exists("resultsPath")) {write(Vdamo2taxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
  }
  ##### END VANDAM OXYGEN REQUIREMENTS

  ##### VANDAM SAPROBITY
  #creates results dataframe
  vdamsap_labels <- c("VD Saprobity 1", "VD Saprobity 2", "VD Saprobity 3", "VD Saprobity 4","VD Saprobity 5", "VD Saprobity Indet", "VD Saprobity Taxa used")
  vdamSaprobity <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(vdamSaprobity) <- vdamsap_labels
  #finds the column
  vdamsap_v <- (taxaInRA[,"vdamsap"])
  #remove the NA
  vdamsap_v[is.na(vdamsap_v)] = 0
  print("Calculating Van Dam saprobity")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamsap_cat1 <- sum(taxaInRA[which(vdamsap_v == 1),sampleNumber])
    vdamsap_cat2 <- sum(taxaInRA[which(vdamsap_v == 2),sampleNumber])
    vdamsap_cat3 <- sum(taxaInRA[which(vdamsap_v == 3),sampleNumber])
    vdamsap_cat4 <- sum(taxaInRA[which(vdamsap_v == 4),sampleNumber])
    vdamsap_cat5 <- sum(taxaInRA[which(vdamsap_v == 5),sampleNumber])
    vdamsap_indet <- 100 - sum(vdamsap_cat1, vdamsap_cat2, vdamsap_cat3, vdamsap_cat4, vdamsap_cat5) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    vdamsaptaxaused <- length(which(vdamsap_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamsaptaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    vdamsap_values <- c(vdamsap_cat1, vdamsap_cat2, vdamsap_cat3, vdamsap_cat4, vdamsap_cat5, vdamsap_indet, vdamsaptaxaused)
    vdamSaprobity[sampleNumber, ] <- vdamsap_values
    if (vandamReports==T & exists("resultsPath")) {write(paste("SAPROBITY- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
    if (vandamReports==T & exists("resultsPath")) {write(Vdamsaptaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
  }
  ##### END VANDAM SAPROBITY

  ##### VANDAM MOISTURE (AERO) STATE
  #creates results dataframe
  vdamaero_labels <- c("VD Aero 1", "VD Aero 2", "VD Aero 3", "VD Aero 4","VD Aero 5", "VD Aero Indet", "VD Aero Taxa used")
  vdamAero <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(vdamAero) <- vdamaero_labels
  #finds the column
  vdamaero_v <- (taxaInRA[,"vdamaero"])
  #remove the NA
  vdamaero_v[is.na(vdamaero_v)] = 0
  print("Calculating Van Dam aero")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamaero_cat1 <- sum(taxaInRA[which(vdamaero_v == 1),sampleNumber])
    vdamaero_cat2 <- sum(taxaInRA[which(vdamaero_v == 2),sampleNumber])
    vdamaero_cat3 <- sum(taxaInRA[which(vdamaero_v == 3),sampleNumber])
    vdamaero_cat4 <- sum(taxaInRA[which(vdamaero_v == 4),sampleNumber])
    vdamaero_cat5 <- sum(taxaInRA[which(vdamaero_v == 5),sampleNumber])
    vdamaero_indet <- 100 - sum(vdamaero_cat1, vdamaero_cat2, vdamaero_cat3, vdamaero_cat4, vdamaero_cat5) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    vdamaerotaxaused <- length(which(vdamaero_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamaerotaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    vdamaero_values <- c(vdamaero_cat1, vdamaero_cat2, vdamaero_cat3, vdamaero_cat4, vdamaero_cat5, vdamaero_indet, vdamaerotaxaused)
    vdamAero[sampleNumber, ] <- vdamaero_values
    if (vandamReports==T & exists("resultsPath")) {write(paste("MOISTURE- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
    if (vandamReports==T & exists("resultsPath")) {write(Vdamaerotaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
  }
  ##### END VANDAM MOISTURE (AERO) STATE

  ##### VANDAM TROPHIC STATE
  #creates results dataframe
  vdamtrop_labels <- c("VD Trophic 1", "VD Trophic 2", "VD Trophic 3", "VD Trophic 4","VD Trophic 5", "VD Trophic Indet", "VD Trophic Taxa used")
  vdamTrophic <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(vdamTrophic) <- vdamtrop_labels
  #finds the column
  vdamtrop_v <- (taxaInRA[,"vdamtrop"])
  #remove the NA
  vdamtrop_v[is.na(vdamtrop_v)] = 0
  print("Calculating Van Dam trophic state")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamtrop_cat1 <- sum(taxaInRA[which(vdamtrop_v == 1),sampleNumber])
    vdamtrop_cat2 <- sum(taxaInRA[which(vdamtrop_v == 2),sampleNumber])
    vdamtrop_cat3 <- sum(taxaInRA[which(vdamtrop_v == 3),sampleNumber])
    vdamtrop_cat4 <- sum(taxaInRA[which(vdamtrop_v == 4),sampleNumber])
    vdamtrop_cat5 <- sum(taxaInRA[which(vdamtrop_v == 5),sampleNumber])
    vdamtrop_indet <- 100 - sum(vdamtrop_cat1, vdamtrop_cat2, vdamtrop_cat3, vdamtrop_cat4, vdamtrop_cat5) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    vdamtroptaxaused <- length(which(vdamtrop_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamtroptaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    vdamtrop_values <- c(vdamtrop_cat1, vdamtrop_cat2, vdamtrop_cat3, vdamtrop_cat4, vdamtrop_cat5, vdamtrop_indet, vdamtroptaxaused)
    vdamTrophic[sampleNumber, ] <- vdamtrop_values
    if (vandamReports==T & exists("resultsPath")) {write(paste("TROPHIC STATE- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
    if (vandamReports==T & exists("resultsPath")) {write(Vdamtroptaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=T)}
  }
  ##### END VANDAM TROPHIC STATE

  ### VANDAM RESULTS SUMMARY TABLES
  vandam.results <- data.frame(c(vdamSalinity, vdamNHeterotrophy, vdamOxygen, vdamSaprobity, vdamAero, vdamTrophic))

  if (vandamReports==T & exists("resultsPath")) {print(paste("VanDam detailed reports - File exported as ", paste(resultsPath, "\\", vandamtaxafile, sep="")))}
  ### END VANDAM RESULTS

  return(vandam.results)
}
