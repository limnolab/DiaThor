#' Calculates the Diatom Index for Soda Pans (DISP)
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Stenger-Kovács, C., Körmendi, K., Lengyel, E., Abonyi, A., Hajnal, É., Szabó, B., Buczkó, K. & Padisák, J. (2018). Expanding the trait-based concept of benthic diatoms: Development of trait-and species-based indices for conductivity as the master variable of ecological status in continental saline lakes. Ecological Indicators, 95, 63-74.
#' }
#'
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi.org/10.1016/j.ecolind.2019.105951
#' }
#' @examples
#' \donttest{
#' # Example using sample data included in the package (sampleData):
#' data("diat_sampleData")
#' # First, the diat_loadData() function has to be called to read the data
#' # The data will be stored into a list (loadedData)
#' # And an output folder will be selected through a dialog box if resultsPath is empty
#' # In the example, a temporary directory will be used in resultsPath
#' df <- diat_loadData(diat_sampleData, resultsPath = tempdir())
#' dispResults <- diat_disp(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_disp


###### ---------- FUNCTION FOR DISP INDEX (Leclercq & Maq. 1988)---------- ########
### INPUT: resultLoad Data
### OUTPUTS: dataframe with DISP index per sample
diat_disp <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[1]] #1 = Relative Abundance, 2 = abundance

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  dispDB <- diathor::disp

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)

  # #the ones still not found (NA), try against fullspecies
  taxaIn$disp_v <- NA
  taxaIn$disp_s <- NA
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$disp_s[i]) | is.na(taxaIn$disp_v[i])){
      taxaIn$disp_v[i] <- dispDB$disp_v[match(trimws(rownames(taxaIn[i,])), trimws(dispDB$fullspecies))]
      taxaIn$disp_s[i] <- dispDB$disp_s[match(trimws(rownames(taxaIn[i,])), trimws(dispDB$fullspecies))]
    }
  }

  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #######--------DISP INDEX START --------#############
  print("Calculating DISP index")
  #creates results dataframe
  disp.results <- data.frame(matrix(ncol = 2, nrow = (lastcol-1)))
  colnames(disp.results) <- c("DISP", "Precision")
  #finds the column
  disp_s <- (taxaIn[,"disp_s"])
  disp_v <- (taxaIn[,"disp_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    DISPtaxaused <- (length(which(disp_s * as.double(taxaIn[,sampleNumber]) > 0))*100 / length(disp_s))
    #remove the NA
    disp_s[is.na(disp_s)] = 0
    disp_v[is.na(disp_v)] = 0
    DISP <- sum((taxaIn[,sampleNumber]*as.double(disp_s)*as.double(disp_v)))/sum(taxaIn[,sampleNumber]*as.double(disp_v)) #raw value
    disp.results[sampleNumber, ] <- c(DISP, DISPtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------DISP INDEX: END--------############
  #PRECISION
  resultsPath <- resultLoad[[4]]
  precisionmatrix <- read.csv(file.path(resultsPath, "Precision.csv"))
  precisionmatrix <- cbind(precisionmatrix, disp.results$Precision)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="disp.results$Precision"] <- "DISP"
  write.csv(precisionmatrix, file.path(resultsPath, "Precision.csv"))
  #END PRECISION

  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(taxaIn$disp_s > 0)]
  inclusionmatrix <- read.csv(file.path(resultsPath, "Taxa included.csv"))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "DISP")
  #creates a new data matrix to append the existing Taxa Included file
  newinclusionmatrix <- as.data.frame(matrix(nrow=max(length(taxaIncluded), nrow(inclusionmatrix)), ncol=ncol(inclusionmatrix)+1))
  for (i in 1:ncol(inclusionmatrix)){
    newinclusionmatrix[1:nrow(inclusionmatrix),i] <- as.character(inclusionmatrix[1:nrow(inclusionmatrix),i])
  }
  if (nrow(newinclusionmatrix) > length(taxaIncluded)){
    newinclusionmatrix[1:length(taxaIncluded), ncol(newinclusionmatrix)] <- taxaIncluded
  } else {
    newinclusionmatrix[1:nrow(newinclusionmatrix), ncol(newinclusionmatrix)] <- taxaIncluded
  }
  inclusionmatrix <- newinclusionmatrix
  colnames(inclusionmatrix) <- colnamesInclusionMatrix
  inclusionmatrix <- inclusionmatrix[-(1:which(colnames(inclusionmatrix)=="Eco.Morpho")-1)]
  write.csv(inclusionmatrix, file.path(resultsPath,"Taxa included.csv"))
  #END TAXA INCLUSION

  #EXCLUDED TAXA
  taxaExcluded <- taxaIn[!('%in%'(taxaIn$species,taxaIncluded)),"species"]
  exclusionmatrix <- read.csv(file.path(resultsPath, "Taxa excluded.csv"))
  #creates a new data matrix to append the existing Taxa Included file
  newexclusionmatrix <- as.data.frame(matrix(nrow=max(length(taxaExcluded), nrow(exclusionmatrix)), ncol=ncol(exclusionmatrix)+1))
  for (i in 1:ncol(exclusionmatrix)){
    newexclusionmatrix[1:nrow(exclusionmatrix),i] <- as.character(exclusionmatrix[1:nrow(exclusionmatrix),i])
  }
  if (nrow(newexclusionmatrix) > length(taxaExcluded)){
    newexclusionmatrix[1:length(taxaExcluded), ncol(newexclusionmatrix)] <- taxaExcluded
  } else {
    newexclusionmatrix[1:nrow(newexclusionmatrix), ncol(newexclusionmatrix)] <- taxaExcluded
  }
  exclusionmatrix <- newexclusionmatrix
  colnames(exclusionmatrix) <- colnamesInclusionMatrix
  exclusionmatrix <- exclusionmatrix[-(1:which(colnames(exclusionmatrix)=="Eco.Morpho")-1)]
  #write.csv(exclusionmatrix, paste(resultsPath,"\\Taxa excluded.csv", sep=""))
  write.csv(exclusionmatrix, file.path(resultsPath,"Taxa excluded.csv"))
  #END EXCLUDED TAXA

  rownames(disp.results) <- resultLoad[[3]]
  return(disp.results)
}

