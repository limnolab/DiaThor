#' This function is the master function of the package. It calculates all possible outputs from the data, and places them in the Output folder
#' @param species_df The data frame with your species data. Species as rows, Sites as columns. If empty, a dialog will prompt for a CSV file
#' @param isRelAb Boolean. If set to 'TRUE' it means that your species' data is the relative abundance of each species per site. If FALSE, it means that it the data corresponds to absolute densities. Default = FALSE
#' @param maxDistTaxa Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2
#' @param calculateguilds Boolean. If set to 'TRUE' the percentage of abundance of each diatom guild will be calculated. Default = TRUE
#' @param vandamReports Boolean. If set to 'TRUE' the detailed reports for the Van Dam classifications will be reported in the Output. Default = TRUE
#' @param singleResult Boolean. If set to 'TRUE' all results will go into a single output file. If FALSE, separate output files will be generated. Default = TRUE
#' @param plotAll Boolean. If set to 'TRUE', plots will be generated for each Output in a PDF file. Default = TRUE
#' @param exportFormat Integer. If = 1: A CSV (external file) will be generated with the output matrices; 2: an internal R dataframe will be generated; 3: both a CSV and an internal R dataframe are generated. Default: 3
#' @param exportName String. Prefix for all export files. Default: "Diato_results"
#' @param color Color code (hex). Default color for bar charts and lolipop plots. Default "#0073C2"
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
#' @export diathorAll



###### ---------- MASTER FUNCTION, CALCULATES EVERYTHING WITH ALL POSSIBLE OUTPUTS BY DEFAULT  ---------- ########

diathorAll <- function(species_df, isRelAb = F, maxDistTaxa = 2, calculateguilds = T, vandam = T, vandamReports = T, singleResult = T, exportFormat = 3, exportName = "Diato_results", plotAll = T, color = "#0073C2"){
  resultmat <- diat_loadData(species_df)
  morpho.results <- diat_morpho(resultmat)
  numcloroplastos.result <- morpho.results[[1]]
  shpcloroplastos.result <- morpho.results[[2]]
  biovol.val.result <- morpho.results[[3]]
  size.results <- diat_size(resultmat)
  diversity.results <- diat_diversity(resultmat)
  guilds.results <- diat_guild(resultmat)
  vandam.results <- diat_vandam(resultmat)
  ips.results <- diat_ips(resultmat)
  tdi.results <- diat_tdi(resultmat)
  idp.results <- diat_idp(resultmat)
  #cee.results <- calculate_cee(resultmat)
  des.results <- diat_des(resultmat)
  epid.results <- diat_epid(resultmat)
  idap.results <- diat_idap(resultmat)
  #idch.results <- diat_idch(resultmat) -- NEEDS ALL THE ACRONYMS IN THE FILE
  lobo.results <- diat_lobo(resultmat)
  #she.results <- calculate_she(resultmat)
  #wat.results <- calculate_wat(resultmat)
  sla.results <- diat_sla(resultmat)
  spear.results <- diat_spear(resultmat)

  #sampledata
  sampleNames <- resultmat[[3]]
  resultsPath <- resultmat[[4]]

  ########### RESULTS TABLES ############
  if (singleResult == T){ #by default exports a single spreadsheet with all the results
    singleTable <- NULL
    singleTable <- data.frame(c(if(exists("diversity.results")){diversity.results},
                     if(exists("numcloroplastos.result")){numcloroplastos.result},
                     if(exists("shpcloroplastos.result")){shpcloroplastos.result},
                     if(exists("biovol.val.result")){biovol.val.result},
                     if(exists("size.results")){size.results},
                     if(exists("guilds.results")){guilds.results},
                     if(exists("vandam.results")){vandam.results},
                     if(exists("ips.results")){ips.results},
                     if(exists("tdi.results")){tdi.results},
                     if(exists("idp.results")){idp.results},
                     if(exists("des.results")){des.results},
                     if(exists("epid.results")){epid.results},
                     if(exists("idap.results")){idap.results},
                     if(exists("idch.results")){idch.results},
                     if(exists("lobo.results")){lobo.results},
                     if(exists("sla.results")){sla.results},
                     if(exists("spear.results")){spear.results}
                     ))

    rownames(singleTable) <- sampleNames
  } else { #separate files for each result
    listOfTables <- list(if(exists("diversity.results")){diversity.results},
                     if(exists("numcloroplastos.result")){numcloroplastos.result},
                     if(exists("shpcloroplastos.result")){shpcloroplastos.result},
                     if(exists("biovol.val.result")){biovol.val.result},
                     if(exists("size.results")){size.results},
                     if(exists("guilds.results")){guilds.results},
                     if(exists("vandam.results")){vandam.results},
                     if(exists("ips.results")){ips.results},
                     if(exists("tdi.results")){tdi.results},
                     if(exists("idp.results")){idp.results},
                     if(exists("des.results")){des.results},
                     if(exists("epid.results")){epid.results},
                     if(exists("idap.results")){idap.results},
                     if(exists("idch.results")){idch.results},
                     if(exists("lobo.results")){lobo.results},
                     if(exists("sla.results")){sla.results},
                     if(exists("spear.results")){spear.results}
    )

    names(listOfTables) <- c(if(exists("diversity.results")){"Diversity"},
                         if(exists("numcloroplastos.result")){"Chloroplast number"},
                         if(exists("shpcloroplastos.result")){"Chloroplast shape"},
                         if(exists("biovol.val.result")){"Biovolume"},
                         if(exists("size.results")){"Size classes"},
                         if(exists("guilds.results")){"Guilds"},
                         if(exists("vandam.results")){"VanDam"},
                         if(exists("ips.results")){"IPS index"},
                         if(exists("tdi.results")){"TDI index"},
                         if(exists("idp.results")){"IDP index"},
                         if(exists("des.results")){"DES index"},
                         if(exists("epid.results")){"EPID index"},
                         if(exists("idap.results")){"IDAP index"},
                         if(exists("idch.results")){"IDCH index"},
                         if(exists("lobo.results")){"LOBO index"},
                         if(exists("sla.results")){"SLA index"},
                         if(exists("spear.results")){"SPEAR index"}
                         )

  }

  #EXPORT AS CSV
  if (exportFormat == 1) {
    if (singleResult == T) {
      filename = paste(exportName, " - Results", ".csv")
      write.csv(singleTable, paste(resultsPath, "\\", filename, sep=""))
    } else {
      for (i in seq_along(listOfTables)) {
        filename = paste(exportName, " - ",names(listOfTables)[i], ".csv")
        write.csv(listOfTables[[i]], paste(resultsPath, "\\", filename, sep=""))
      }
    }
  }

  #EXPORT AS INTERNAL DATAFRAME
  if (exportFormat == 2) {
    if (singleResult == T) {
      .GlobalEnv$Results <- singleTable
    } else {
      .GlobalEnv$Diversity <- as.data.frame(listOfTables[[1]])
      .GlobalEnv$ChloroplastNumber <- as.data.frame(listOfTables[[2]])
      .GlobalEnv$ChloroplastShape <- as.data.frame(listOfTables[[3]])
      .GlobalEnv$Biovolume <- as.data.frame(listOfTables[[4]])
      .GlobalEnv$SizeClasses <- as.data.frame(listOfTables[[5]])
      .GlobalEnv$Guilds <- as.data.frame(listOfTables[[6]])
      .GlobalEnv$VanDam <- as.data.frame(listOfTables[[7]])
      .GlobalEnv$IPS <- as.data.frame(listOfTables[[8]])
      .GlobalEnv$TDI <- as.data.frame(listOfTables[[9]])
      .GlobalEnv$DES <- as.data.frame(listOfTables[[10]])
      .GlobalEnv$EPID <- as.data.frame(listOfTables[[11]])
      .GlobalEnv$IDAP <- as.data.frame(listOfTables[[12]])
      .GlobalEnv$IDCH <- as.data.frame(listOfTables[[13]])
      .GlobalEnv$LOBO <- as.data.frame(listOfTables[[14]])
      .GlobalEnv$SLA <- as.data.frame(listOfTables[[15]])
      .GlobalEnv$SPEAR <- as.data.frame(listOfTables[[16]])
    }
  }


  #EXPORT AS BOTH CSV AND INTERNAL DATAFRAME - Default
  if (exportFormat == 3) {
    if (singleResult == T) {
      filename = paste(exportName, " - Results", ".csv")
      write.csv(singleTable, paste(resultsPath, "\\", filename, sep=""))
      .GlobalEnv$Results <- singleTable
    } else {
      .GlobalEnv$Diversity <- as.data.frame(listOfTables[[1]])
      .GlobalEnv$ChloroplastNumber <- as.data.frame(listOfTables[[2]])
      .GlobalEnv$ChloroplastShape <- as.data.frame(listOfTables[[3]])
      .GlobalEnv$Biovolume <- as.data.frame(listOfTables[[4]])
      .GlobalEnv$SizeClasses <- as.data.frame(listOfTables[[5]])
      .GlobalEnv$Guilds <- as.data.frame(listOfTables[[6]])
      .GlobalEnv$VanDam <- as.data.frame(listOfTables[[7]])
      .GlobalEnv$IPS <- as.data.frame(listOfTables[[8]])
      .GlobalEnv$TDI <- as.data.frame(listOfTables[[9]])
      .GlobalEnv$IDP <- as.data.frame(listOfTables[[10]])
      .GlobalEnv$EPID <- as.data.frame(listOfTables[[11]])
      .GlobalEnv$IDAP <- as.data.frame(listOfTables[[12]])
      .GlobalEnv$IDCH <- as.data.frame(listOfTables[[13]])
      .GlobalEnv$LOBO <- as.data.frame(listOfTables[[14]])
      .GlobalEnv$SLA <- as.data.frame(listOfTables[[15]])
      .GlobalEnv$SPEAR <- as.data.frame(listOfTables[[16]])
      for (i in seq_along(listOfTables)) {
        filename = paste(exportName, " - ",names(listOfTables)[i], ".csv")
        write.csv(listOfTables[[i]], paste(resultsPath, "\\", filename, sep=""))
      }
    }
  }
  ########### END RESULTS TABLES ############

  ########### START PLOTS ############

  loli.plot <- function(result, ylabel, ylow, yhigh, color = "#0073C2"){
    x <- rownames(result)
    data <- result
    data[is.na(data)] = 0

    for(i in 1:ncol(result)) {
      y <- data[,i]

      print(ggplot2::ggplot(data, aes(x=x, y=y)) +
              geom_segment( aes(x=x, xend=x, y=0, yend=y), color="grey") +
              geom_point( color=color, size=4) +
              theme_light() +
              theme(
                panel.grid.major.x = element_blank(),
                panel.border = element_blank(),
                axis.ticks.x = element_blank()
              ) +
              ylim(ylow, yhigh) +
              xlab(rownames(result)) +
              ylab(ylabel)
      )
    }
  }

  percentbarchart.plot <- function(result, title){
    x <- colnames(result) #checks if taxa used colum exists and removes it
    x <- substr(x, nchar(x)-8, nchar(x)) == 'Taxa used' #checks if taxa used colum exists and removes it
    if (tail(x, n=1)==T) {
      result <- result[,-ncol(result)] #and removes it
    }
    sampleCol <- rep(sampleNames, ncol(result)) #gets sample names
    result <- tidyr::gather(result) #uses tidyr to rearrange the dataframe in a single column
    result$sampleCol <- sampleCol #adds another column with the sample names

    colors <- c("#CC1C00", "#5C88DA", "#84BD00", "#FFCD00", "#7C878E", "#E64B35", "#4DBBD5", "#01A087", "#3C5488", "#F39B7F")

    print(ggplot2::ggplot(result, aes(fill=key, y=value, x=sampleCol)) +
            geom_bar(position="fill", stat="identity") +
            scale_fill_manual(values=colors) +
            xlab("Samples") +
            ylab(title)) #graphs
  }

  result.plots.allToPDF <- function(color ="#0073C2"){
    print("Exporting all plots to PDF, please wait...")
    # Open a pdf file
    pdf(paste(resultsPath, "\\", "ResultPlots.pdf", sep=""))

    #Plots all resulting graphs (if exist)
    if(exists("diversity.results")){

      loli.plot(as.data.frame(diversity.results[,1]), "Species richess", 0, max(as.data.frame(diversity.results[,1])))
      loli.plot(as.data.frame(diversity.results[,2]), "Shannon's diversity", 0, max(as.data.frame(diversity.results[,2])))
      loli.plot(as.data.frame(diversity.results[,3]), "Evenness", 0, max(as.data.frame(diversity.results[,3])))
    }
    if(exists("biovol.val.result")){
      loli.plot(as.data.frame(biovol.val.result[,1]), "Biovolume", 0, max(as.data.frame(biovol.val.result[,1])))
    }
    if(exists("numcloroplastos.result")){
      percentbarchart.plot(numcloroplastos.result, "Number of chloroplasts") #default: piled bars
    }
    if(exists("shpcloroplastos.result")){
      percentbarchart.plot(shpcloroplastos.result, "Shape of chloroplasts") #default: piled bars
    }
    if(exists("size.results")){
      percentbarchart.plot(size.results, "Size classes") #default: piled bars
    }
    if(exists("guilds.results")){
      percentbarchart.plot(guilds.results, "Guilds")#default: piled bars
    }
    if(exists("vdamSalinity")){
      percentbarchart.plot(vdamSalinity, "Salinity")
    }
    if(exists("vdamNHeterotrophy")){
      percentbarchart.plot(vdamNHeterotrophy, "N-Heterotrophy")
    }
    if(exists("vdamOxygen")){
      percentbarchart.plot(vdamOxygen, "Oxygen preferences")
    }
    if(exists("vdamSaprobity")){
      percentbarchart.plot(vdamSaprobity, "Saprobity")
    }
    if(exists("vdamAero")){
      percentbarchart.plot(vdamAero, "Moisture")
    }
    if(exists("vdamTrophic")){
      percentbarchart.plot(vdamTrophic, "Trophic state")
    }

    if(exists("ips.results")){
      loli.plot(as.data.frame(ips.results[,1]), "IPS", 0, 10) #raw index
      loli.plot(as.data.frame(ips.results[,2]), "IPS - Standardized", 0, 20) #standard 20
    }
    if(exists("tdi.results")){
      loli.plot(as.data.frame(tdi.results[,1]), "TDI - Standardized", 0, 20) #standard 20
      loli.plot(as.data.frame(tdi.results[,2]), "TDI - %", 0, 100) #%
    }
    if(exists("idp.results")){
      loli.plot(as.data.frame(idp.results[,1]), "IDP", 0, 10) #raw index
      loli.plot(as.data.frame(idp.results[,2]), "IDP - Standardized", 0, 20) #standard 20
    }
    if(exists("epid.results")){
      loli.plot(as.data.frame(epid.results[,1]), "EPID", 0, 10) #raw index
      loli.plot(as.data.frame(epid.results[,2]), "EPID - Standardized", 0, 20) #standard 20
    }
    if(exists("idap.results")){
      loli.plot(as.data.frame(idap.results[,1]), "IDAP", 0, 10) #raw index
      loli.plot(as.data.frame(idap.results[,2]), "IDAP - Standardized", 0, 20) #standard 20
    }
    if(exists("idch.results")){
      loli.plot(as.data.frame(idch.results[,1]), "IDCH", 0, 10) #raw index
      loli.plot(as.data.frame(idch.results[,2]), "IDCH - Standardized", 0, 20) #standard 20
    }
    if(exists("lobo.results")){
      loli.plot(as.data.frame(lobo.results[,1]), "LOBO", 0, 10) #raw index
      loli.plot(as.data.frame(lobo.results[,2]), "LOBO - Standardized", 0, 20) #standard 20
    }
    if(exists("sla.results")){
      loli.plot(as.data.frame(sla.results[,1]), "SLA", 0, 10) #raw index
      loli.plot(as.data.frame(sla.results[,2]), "SLA - Standardized", 0, 20) #standard 20
    }
    if(exists("spear.results")){
      loli.plot(as.data.frame(spear.results[,1]), "SPEAR", 0, 10) #raw index
    }

    # Close the pdf file
    dev.off()
    print("Plots exported!")

  }

  if (plotAll == TRUE){
    result.plots.allToPDF()
  }

  ########### END PLOTS ############




}
