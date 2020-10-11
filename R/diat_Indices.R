
###### ---------- FUNCTION FOR IPS INDEX: NEW TEST DONE   ---------- ########
#### IN THIS SECTION WE CALCULATE IPS INDEX (REF)
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with IPS index per sample
#' @export
diat_ips <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaIn)=="species")

  #######--------IPS INDEX START (indice poluto sensible)--------#############
  print("Calculating IPS index")
  #creates results dataframe
  ips.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(ips.results) <- c("IPS", "IPS20", "IPS Taxa used")
  #finds the column
  ips_s <- (taxaIn[,"ipss"])
  ips_vv <- (taxaIn[,"ipsvv"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    IPStaxaused <- (length(ips_s) - sum(is.na(taxaIn$ips_s )))*100 / length(ips_s)
    #remove the NA
    ips_s[is.na(ips_s)] = 0
    ips_vv[is.na(ips_vv)] = 0
    IPS <- sum((taxaIn[,sampleNumber]*as.double(ips_s)*as.double(ips_vv)))/sum(taxaIn[,sampleNumber]*as.double(ips_vv)) #raw value
    IPS20 <- (IPS*4.75)-3.75 #STANDARDIZED VALUE TO 20
    ips.results[sampleNumber, ] <- c(IPS, IPS20,IPStaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------IPS INDEX: END--------############
  return(ips.results)
}


###### ---------- FUNCTION FOR TDI INDEX: NEW TEST DONE  (Trophic Index - Kelly)  ---------- ########
#### IN THIS SECTION WE CALCULATE TDI INDEX (Trophic Index - Kelly
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with TDI index per sample
#### WARNING: The formula used is [TDI= sum(a.s.v)/ sum(a.v)], NOT the formula
#### written in the OMNIDIA help file, which is [TDI= sum(a.s.v)/ sum(s.v)].
#### This corrected formula matches the OMNIDIA results.
#' @export
diat_tdi <- function(resultLoad){


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

  #######--------TDI INDEX START --------#############
  print("Calculating TDI index")
  #creates results dataframe
  tdi.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(tdi.results) <- c("TDI20", "TDI100", "TDI Taxa used")
  #finds the column
  tdi_s <- (taxaIn[,"kellys"])
  tdi_vv <- (taxaIn[,"kellyv"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    TDItaxaused <- (length(tdi_s) - sum(is.na(taxaIn$tdi_s )))*100 / length(tdi_s)
    #remove the NA
    tdi_s[is.na(tdi_s)] = 0
    tdi_vv[is.na(tdi_vv)] = 0
    TDI <- sum((taxaIn[,sampleNumber]*as.double(tdi_s)*as.double(tdi_vv)))/sum(taxaIn[,sampleNumber]*as.double(tdi_vv)) #raw value
    TDI20 <- (-4.75*TDI)+24.75
    TDI100 <- (TDI*25)-25
    tdi.results[sampleNumber, ] <- c(TDI20, TDI100,TDItaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------TDI INDEX: END--------############
  return(tdi.results)
}


###### ---------- FUNCTION FOR IDP INDEX: TEST IN OMNIDIA (Pampean Index - Gomez & Licursi)  ---------- ########
#### IN THIS SECTION WE CALCULATE IDP INDEX (Pampean Index - Gomez & Licursi)
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with IDP index per sample
#' @export
diat_idp <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #idpDB <- read.csv("../Indices/idp.csv") #uses the external csv file
  idpDB <- DiaThor::idp
  idpDB2 <- idpDB
  idpDB <- idpDB

  #exact matches species in input data to acronym from index
  taxaIn$idp_v <- idpDB$idp_v[match(taxaIn$acronym, trimws(idpDB$acronym))]


  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$idp_v[i])){
      taxaIn$idp_v[i] <- idpDB$idp_v[match(taxaIn$species[i], trimws(idpDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS


  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaIn)=="species")

  #######--------IDP INDEX START --------#############
  print("Calculating IDP index")
  idp.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(idp.results) <- c("IDP", "IDP20", "IDP Taxa used")
  #finds the column
  idp_v <- (taxaIn[,"idp_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    IDPtaxaused <- (length(idp_v) - sum(is.na(taxaIn$idp_v )))*100 / length(idp_v)
    #remove the NA
    idp_v[is.na(idp_v)] = 0
    IDP <- sum((taxaIn[,sampleNumber]*as.double(idp_v)))/sum(taxaIn[which(idp_v > 0),sampleNumber]) #raw value
    IDP20 <- 20-(4.75*IDP)
    idp.results[sampleNumber, ] <- c(IDP, IDP20,IDPtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #outputs: IDP, IDP20
  #close progressbar
  close(pb)
  #######--------IDP INDEX: END--------############
  return(idp.results)
}


###### ---------- FUNCTION FOR CEE INDEX: INCOMPLETE: CANT CALCULATE?  (Descy & Coste 1991)---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with CEE index per sample
#' @export
diat_cee <- function(resultLoad){

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

  #######--------CEE INDEX START --------#############
  print("Calculating CEE index")
  cee.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(cee.results) <- c("CEE", "CEE20", "CEE Taxa used")
  #finds the column
  cee_v <- (taxaIn[,"cee_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    CEEtaxaused <- (length(cee_v) - sum(is.na(cee_v)))
    #remove the NA
    cee_v[is.na(cee_v)] = 0
    CEE <- sum((taxaIn[,sampleNumber]*as.double(cee_v)))/sum(taxaIn[which(cee_v > 0),sampleNumber]) #raw value
    CEE20 <- 1+(1.91*CEE)
    cee.results[sampleNumber, ] <- c(CEE, CEE20,CEEtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #outputs: IDP, IDP20
  #close progressbar
  close(pb)
  #######--------CEE INDEX: END--------############
  return(cee.results)
}

###### ---------- FUNCTION FOR DES INDEX: NEW TEST DONE (Descy 1979) ---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with DES index per sample
#' @export
diat_des <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]


  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #desDB <- read.csv("../Indices/des.csv") #uses the external csv file
  desDB <- DiaThor::des
  #exact matches species in input data to acronym from index
  taxaIn$des_v <- desDB$des_v[match(taxaIn$acronym, trimws(desDB$acronym))]
  taxaIn$des_s <- desDB$des_s[match(taxaIn$acronym, trimws(desDB$acronym))]


  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$des_s[i]) | is.na(taxaIn$des_v[i])){
      taxaIn$des_v[i] <- desDB$des_v[match(taxaIn$species[i], trimws(desDB$fullspecies))]
      taxaIn$des_s[i] <- desDB$des_s[match(taxaIn$species[i], trimws(desDB$fullspecies))]
      #taxaIn$des_v[i] <- desDB$des_v[match(row.names(taxaIn[i,]), trimws(desDB$fullspecies))]
      #taxaIn$des_s[i] <- desDB$des_s[match(row.names(taxaIn[i,]), trimws(desDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS

  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="species")

  #######--------DES INDEX START --------#############
  print("Calculating DES index")
  #creates results dataframe
  des.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(des.results) <- c("DES", "DES20", "DES Taxa used")
  #finds the column
  des_s <- (taxaIn[,"des_s"])
  des_v <- (taxaIn[,"des_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    DEStaxaused <- (length(des_s) - sum(is.na(taxaIn$des_s )))*100 / length(des_s)

    #remove the NA
    des_s[is.na(des_s)] = 0
    des_v[is.na(des_v)] = 0
    DES <- sum((taxaIn[,sampleNumber]*as.double(des_s)*as.double(des_v)))/sum(taxaIn[,sampleNumber]*as.double(des_v)) #raw value
    DES20 <- (4.75*DES)-3.75
    des.results[sampleNumber, ] <- c(DES, DES20,DEStaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------DES INDEX: END--------############
  return(des.results)
}

###### ---------- FUNCTION FOR EPID INDEX: NEW TEST DONE (Dell'Uomo) ---------- ########
### INPUT: resultLoad Data. Data needs to be in RA for this index, so if it isn't, the function converts it
### OUTPUTS: dataframe with EPID index per sample
#' @export
diat_epid <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]] #input data
  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #epidDB <- read.csv("../Indices/epid.csv") #uses the external csv file
  epidDB <- DiaThor::epid
  #exact matches species in input data to acronym from index
  taxaIn$epid_v <- epidDB$epid_v[match(taxaIn$acronym, trimws(epidDB$acronym))]
  taxaIn$epid_s <- epidDB$epid_s[match(taxaIn$acronym, trimws(epidDB$acronym))]


  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$epid_s[i]) | is.na(taxaIn$epid_v[i])){
      taxaIn$epid_v[i] <- epidDB$epid_v[match(taxaIn$species[i], trimws(epidDB$fullspecies))]
      taxaIn$epid_s[i] <- epidDB$epid_s[match(taxaIn$species[i], trimws(epidDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaIn)=="species")

  #######--------EPID INDEX START --------#############
  print("Calculating EPID index")
  #creates results dataframe
  epid.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(epid.results) <- c("EPID", "EPID20", "EPID Taxa used")
  #finds the column
  epid_s <- (taxaIn[,"epid_s"])
  epid_v <- (taxaIn[,"epid_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    EPIDtaxaused <- (length(epid_s) - sum(is.na(taxaIn$epid_s )))*100 / length(epid_s)
    #remove the NA
    epid_s[is.na(epid_s)] = 0
    epid_v[is.na(epid_v)] = 0
    EPID <- sum((taxaIn[,sampleNumber]*as.double(epid_s)*as.double(epid_v)))/sum(taxaIn[,sampleNumber]*as.double(epid_v)) #raw value
    EPID20 <- (-4.75*EPID)+20
    epid.results[sampleNumber, ] <- c(EPID, EPID20,EPIDtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------EPID INDEX: END--------############
  return(epid.results)
}

###### ---------- FUNCTION FOR IDAP INDEX: NEW TEST DONE (Indice Diatomique Artois-Picardie) ---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with IDAP index per sample
#' @export
diat_idap <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #idapDB <- read.csv("../Indices/idap.csv") #uses the external csv file
  idapDB <- DiaThor::idap
  #exact matches species in input data to acronym from index
  taxaIn$idap_v <- idapDB$idap_v[match(taxaIn$acronym, trimws(idapDB$acronym))]
  taxaIn$idap_s <- idapDB$idap_s[match(taxaIn$acronym, trimws(idapDB$acronym))]


  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$idap_s[i]) | is.na(taxaIn$idap_v[i])){
      taxaIn$idap_v[i] <- idapDB$idap_v[match(taxaIn$species[i], trimws(idapDB$fullspecies))]
      taxaIn$idap_s[i] <- idapDB$idap_s[match(taxaIn$species[i], trimws(idapDB$fullspecies))]
    }
  }


  ### FINISH NEW CORRECTIONS


  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaIn)=="species")

  #######--------IDAP INDEX START --------#############
  print("Calculating IDAP index")
  #creates results dataframe
  idap.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(idap.results) <- c("IDAP", "IDAP20", "IDAP Taxa used")
  #finds the column
  idap_s <- (taxaIn[,"idap_s"])
  idap_v <- (taxaIn[,"idap_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    IDAPtaxaused <- (length(idap_s) - sum(is.na(taxaIn$idap_s )))*100 / length(idap_s)

    #remove the NA
    idap_s[is.na(idap_s)] = 0
    idap_v[is.na(idap_v)] = 0
    IDAP <- sum((taxaIn[,sampleNumber]*as.double(idap_s)*as.double(idap_v)))/sum(taxaIn[,sampleNumber]*as.double(idap_v)) #raw value
    IDAP20 <- (4.75*IDAP)-3.75
    idap.results[sampleNumber, ] <- c(IDAP, IDAP20,IDAPtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------IDAP INDEX: END--------############
  return(idap.results)
}



###### ---------- FUNCTION FOR IDCH INDEX: NEED ACRONYMS IN CVS FILE (Swiss Diatom Index)  ---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with IDCH index per sample
#' @export
diat_idch <- function(resultLoad){
  print("starting idch")
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #idchDB <- read.csv("../Indices/idch.csv") #uses the external csv file
  idchDB <- DiaThor::idch
  #exact matches species in input data to acronym from index
  taxaIn$idch_v <- idchDB$idch_v[match(taxaIn$acronym, trimws(idchDB$acronym))]
  taxaIn$idch_s <- idchDB$idch_s[match(taxaIn$acronym, trimws(idchDB$acronym))]
  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$idch_s[i]) | is.na(taxaIn$idch_v[i])){
      taxaIn$idch_v[i] <- idchDB$idch_v[match(taxaIn$species[i], trimws(idchDB$fullspecies))]
      taxaIn$idch_s[i] <- idchDB$idch_s[match(taxaIn$species[i], trimws(idchDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS
  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaIn)=="species")

  #######--------IDCH INDEX START --------#############
  print("Calculating IDCH index")
  #creates results dataframe
  idch.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(idch.results) <- c("IDCH", "IDCH20", "IDCH Taxa used")
  #finds the column
  idch_s <- (taxaIn[,"idch_s"])
  idch_v <- (taxaIn[,"idch_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    IDCHtaxaused <- (length(idch_s) - sum(is.na(taxaIn$idch_s )))*100 / length(idch_s)
    #remove the NA
    idch_s[is.na(idch_s)] = 0
    idch_v[is.na(idch_v)] = 0
    IDCH <- sum((taxaIn[,sampleNumber]*as.double(idch_s)*as.double(idch_v)))/sum(taxaIn[,sampleNumber]*as.double(idch_v)) #raw value
    IDCH20 <- 22.714-(2.714*IDCH)
    idch.results[sampleNumber, ] <- c(IDCH, IDCH20,IDCHtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------IDCH INDEX: END--------############
  print("starting end")
  return(idch.results)

}

###### ---------- FUNCTION FOR ILM INDEX: NEW TEST DONE (Leclercq & Maq. 1988)---------- ########
### INPUT: resultLoad Data
### OUTPUTS: dataframe with ILM index per sample
#' @export
diat_ilm <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #ilmDB <- read.csv("../Indices/ilm.csv") #uses the external csv file
  ilmDB <- DiaThor::ilm
  #exact matches species in input data to acronym from index
  taxaIn$ilm_v <- ilmDB$ilm_v[match(taxaIn$acronym, trimws(ilmDB$acronym))]
  taxaIn$ilm_s <- ilmDB$ilm_s[match(taxaIn$acronym, trimws(ilmDB$acronym))]


  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$ilm_s[i]) | is.na(taxaIn$ilm_v[i])){
      taxaIn$ilm_v[i] <- ilmDB$ilm_v[match(taxaIn$species[i], trimws(ilmDB$fullspecies))]
      taxaIn$ilm_s[i] <- ilmDB$ilm_s[match(taxaIn$species[i], trimws(ilmDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaIn)=="species")

  #######--------ILM INDEX START --------#############
  print("Calculating ILM index")
  #creates results dataframe
  ilm.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(ilm.results) <- c("ILM", "ILM20", "ILM Taxa used")
  #finds the column
  ilm_s <- (taxaIn[,"ilm_s"])
  ilm_v <- (taxaIn[,"ilm_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    ILMtaxaused <- (length(ilm_s) - sum(is.na(taxaIn$ilm_s )))*100 / length(ilm_s)

    #remove the NA
    ilm_s[is.na(ilm_s)] = 0
    ilm_v[is.na(ilm_v)] = 0
    ILM <- sum((taxaIn[,sampleNumber]*as.double(ilm_s)*as.double(ilm_v)))/sum(taxaIn[,sampleNumber]*as.double(ilm_v)) #raw value
    ILM20 <- (4.75*ILM)-3.75
    ilm.results[sampleNumber, ] <- c(ILM, ILM20,ILMtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------ILM INDEX: END--------############
  return(ilm.results)
}

###### ---------- FUNCTION FOR LOBO INDEX: NEW TEST DONE (Lobo et al. 2002)---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with LOBO index per sample
#' @export
diat_lobo <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]


  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #loboDB <- read.csv("../Indices/lobo.csv") #uses the external csv file
  loboDB <- DiaThor::lobo
  #exact matches species in input data to acronym from index
  taxaIn$lobo_v <- loboDB$lobo_v[match(taxaIn$acronym, trimws(loboDB$acronym))]
  taxaIn$lobo_s <- loboDB$lobo_s[match(taxaIn$acronym, trimws(loboDB$acronym))]


  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$lobo_s[i]) | is.na(taxaIn$lobo_v[i])){
      taxaIn$lobo_v[i] <- loboDB$lobo_v[match(taxaIn$species[i], trimws(loboDB$fullspecies))]
      taxaIn$lobo_s[i] <- loboDB$lobo_s[match(taxaIn$species[i], trimws(loboDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS

  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaIn)=="species")

  #######--------LOBO INDEX START --------#############
  print("Calculating LOBO index")
  #creates results dataframe
  lobo.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(lobo.results) <- c("LOBO", "LOBO20", "LOBO Taxa used")
  #finds the column
  lobo_s <- (taxaIn[,"lobo_s"])
  lobo_v <- (taxaIn[,"lobo_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    LOBOtaxaused <- (length(lobo_s) - sum(is.na(taxaIn$lobo_s )))*100 / length(lobo_s)
    #remove the NA
    lobo_s[is.na(lobo_s)] = 0
    lobo_v[is.na(lobo_v)] = 0
    LOBO <- sum((taxaIn[,sampleNumber]*as.double(lobo_s)*as.double(lobo_v)))/sum(taxaIn[,sampleNumber]*as.double(lobo_v)) #raw value
    LOBO20 <- (6.333*LOBO)-5.333
    lobo.results[sampleNumber, ] <- c(LOBO, LOBO20,LOBOtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------LOBO INDEX: END--------############
  return(lobo.results)
}

###### ---------- FUNCTION FOR SHE INDEX: INCOMPLETE: CANT CALCULATE?  (Schiefele & Schreiner)---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with SHE index per sample
#' @export
diat_she <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaInRA <- resultLoad[[1]]

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #sheDB <- read.csv("../Indices/she.csv") #uses the external csv file
  sheDB <- DiaThor::she
  #exact matches species in input data to acronym from index
  taxaInRA$she_v <- sheDB$she_v[match(taxaInRA$acronym, trimws(sheDB$acronym))]

  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaInRA)) {
    if (is.na(taxaInRA$she_v[i])){
      taxaInRA$she_v[i] <- sheDB$she_v[match(taxaInRA$species[i], trimws(sheDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS

  #removes NA from taxaInRA
  taxaInRA[is.na(taxaInRA)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaInRA)=="species")

  #######--------SHE INDEX START --------#############
  print("Calculating SHE index")
  #creates results dataframe
  she.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(she.results) <- c("SHE", "SHE20", "SHE Taxa used")
  #finds the column
  she_v <- (taxaInRA[,"she_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    SHEtaxaused <- (length(she_v) - sum(is.na(taxaIn$she_v )))*100 / length(she_v)

    #remove the NA
    she_v[is.na(she_v)] = 0
    SHE <- sum((taxaInRA[,sampleNumber]*as.double(she_v)))/sum(taxaInRA[,sampleNumber]*as.double(she_v)) #raw value
    SHE20 <- (-4.75*SHE)-3.75
    she.results[sampleNumber, ] <- c(SHE, SHE20,SHEtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------SHE INDEX: END--------############
  return(she.results)
}

###### ---------- FUNCTION FOR WAT INDEX: INCOMPLETE: CANT CALCULATE? (Watanabe et al. 1990) ---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with WAT index per sample
#' @export
diat_wat <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaInRA <- resultLoad[[1]]

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #watDB <- read.csv("../Indices/wat.csv") #uses the external csv file
  watDB <- DiaThor::wat
  #exact matches species in input data to acronym from index
  taxaInRA$wat_v <- watDB$wat_v[match(taxaInRA$acronym, trimws(watDB$acronym))]

  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaInRA)) {
    if (is.na(taxaInRA$wat_v[i])){
      taxaInRA$wat_v[i] <- watDB$wat_v[match(taxaInRA$species[i], trimws(watDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS


  #removes NA from taxaInRA
  taxaInRA[is.na(taxaInRA)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaInRA)=="species")

  #######--------WAT INDEX START --------#############
  print("Calculating WAT index")
  #creates results dataframe
  wat.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(wat.results) <- c("WAT", "WAT20", "WAT Taxa used")

  #finds the column
  wat_v <- as.data.frame(taxaInRA[,"wat_v"])
  colnames(wat_v) <- "wat_v"
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    WATtaxaused <- (length(she_v) - sum(is.na(taxaIn$wat_v )))*100 / length(wat_v)
    #remove the NA
    wat_v[is.na(wat_v)] = 0
    WAT <- 50+ 0.5*(sum(taxaInRA[which(wat_v$wat_v=="2"),sampleNumber])- sum(taxaInRA[which(wat_v$wat_v=="1"),sampleNumber]))
    #WAT <- sum((taxaInRA[,sampleNumber]*as.double(wat_v)))/sum(taxaInRA[,sampleNumber]*as.double(wat_v)) #raw value
    WAT20 <- (0.190*WAT)+1
    wat.results[sampleNumber, ] <- c(WAT, WAT20,WATtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------WAT INDEX: END--------############
  return(wat.results)
}

###### ---------- FUNCTION FOR SLA INDEX: NEW TEST DONE  (Sladecek 1986)---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with SLA index per sample
#' @export
diat_sla <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #slaDB <- read.csv("../Indices/sla.csv") #uses the external csv file
  slaDB <- DiaThor::sla
  #exact matches species in input data to acronym from index
  taxaIn$sla_v <- slaDB$sla_v[match(taxaIn$acronym, trimws(slaDB$acronym))]
  taxaIn$sla_s <- slaDB$sla_s[match(taxaIn$acronym, trimws(slaDB$acronym))]


  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$sla_s[i]) | is.na(taxaIn$sla_v[i])){
      taxaIn$sla_v[i] <- slaDB$sla_v[match(taxaIn$species[i], trimws(slaDB$fullspecies))]
      taxaIn$sla_s[i] <- slaDB$sla_s[match(taxaIn$species[i], trimws(slaDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS

  #removes NA from taxaInRA
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaIn)=="species")

  #######--------SLA INDEX START --------#############
  print("Calculating SLA index")
  #creates results dataframe
  sla.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(sla.results) <- c("SLA", "SLA20", "SLA Taxa used")
  #finds the column
  sla_s <- (taxaIn[,"sla_s"])
  sla_v <- (taxaIn[,"sla_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    SLAtaxaused <- (length(sla_s) - sum(is.na(taxaIn$sla_s )))*100 / length(sla_s)

    #remove the NA
    sla_s[is.na(sla_s)] = 0
    sla_v[is.na(sla_v)] = 0
    SLA <- sum((taxaIn[,sampleNumber]*as.double(sla_s)*as.double(sla_v)))/sum(taxaIn[,sampleNumber]*as.double(sla_v)) #raw value
    SLA20 <- 20-(4.75*SLA)
    sla.results[sampleNumber, ] <- c(SLA, SLA20,SLAtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------SLA INDEX: END--------############
  return(sla.results)
}

###### ---------- FUNCTION FOR SPEAR INDEX: (Wood et al. 2019) ---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with SPEAR index per sample
#' @export
diat_spear <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaInRA <- resultLoad[[1]]

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #spearDB <- read.csv("../Indices/spear.csv") #uses the external csv file
  spearDB <- DiaThor::spear
  #exact matches species in input data to acronym from index
  taxaInRA$spear_v <- spearDB$spear_v[match(taxaInRA$acronym, trimws(spearDB$acronym))]

  #the ones not found exact (NA), try against fullspecies
  for (i in 1:nrow(taxaInRA)) {
    if (is.na(taxaInRA$spear_v[i])){
      taxaInRA$spear_v[i] <- spearDB$spear_v[match(taxaInRA$species[i], trimws(spearDB$fullspecies))]
    }
  }
  ### FINISH NEW CORRECTIONS

  #removes NA from taxaInRA
  taxaInRA[is.na(taxaInRA)] <- 0

  #gets the column named "species", everything before that is a sample
  lastcol = which(colnames(taxaInRA)=="species")

  #######--------SPEAR INDEX START --------#############
  print("Calculating SPEAR index")
  #creates results dataframe
  spear.results <- data.frame(matrix(ncol = 2, nrow = (lastcol-1)))
  colnames(spear.results) <- c("SPEAR", "SPEAR Taxa used")
  #finds the column
  spear_v <- (taxaInRA[,"spear_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    SPEARtaxaused <- (length(spear_v) - sum(is.na(taxaInRA$spear_v )))*100 / length(spear_v)
    #remove the NA
    spear_v[is.na(spear_v)] = 0

    SPEAR <- sum((log(taxaInRA[,sampleNumber]+1)*as.double(spear_v)))/sum(log(taxaInRA[,sampleNumber]+1)) #raw value
    spear.results[sampleNumber, ] <- c(SPEAR, SPEARtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------SPEAR INDEX: END--------############
  return(spear.results)
}
