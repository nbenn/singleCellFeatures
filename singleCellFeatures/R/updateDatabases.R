#' Retrieve the settingsDatabase object 
#' 
#' In order to update the settingsDatabase object, this helper function fetches
#' the current settingsDatabase object for the user to modify and save using
#' updateDatabaseSettingsSet
#'
#' @return The current settingsDatabase object
#'
#' @examples
#' # load current settingsDatabase
#' current <- updateDatabaseSettingsGet()
#' # save updated version
#' updateDatabaseSettingsSet(current, list(newKey="newValue"))
#' # using the settingsDatabase object
#' data(settingsDatabase)
#' str(settings.database)
#' # just updating the settingsDatabase object
#' updateDatabaseSettingsSet(updateDatabaseSettingsGet())
#' 
#' @export
updateDatabaseSettingsGet <- function() {
  # load current settingsDatabase
  data(settingsDatabase, envir=environment())
  return(settings.database)
}

#' Save the settingsDatabase object
#' 
#' In order to update the settingsDatabase object, updateDatabaseSettingsGet
#' can be used to fetch the current settingsDatabase object, which then can be
#' modifyed and saved using this function
#'
#' @param old.list Required parameter holding the old list/a new version of
#'                 settingsDatabase with e.g. a key-value pair removed
#' @param new.list Optional parameter holding one or more key-value pairs as a
#'                 list which will be appended to old.list
#'
#' @return The updated settingsDatabase object which is also saved to the /data
#'         folder. For the new file to be available, the package has to be re-
#'         loaded
#'
#' @examples
#' # load current settingsDatabase
#' current <- updateDatabaseSettingsGet()
#' # save updated version
#' updateDatabaseSettingsSet(current, list(newKey="newValue"))
#' # using the settingsDatabase object
#' data(settingsDatabase)
#' str(settings.database)
#' # just updating the settingsDatabase object
#' updateDatabaseSettingsSet(updateDatabaseSettingsGet())
#' 
#' @export
updateDatabaseSettingsSet <- function(old.list, new.list=NULL) {
  # concatenante current settingsDatabase and new key-value list
  settings.database <- c(old.list, new.list)
  save(settings.database,
       file=paste0(settings.database$package, "/data/settingsDatabase.rda"),
       compression_level=1)
}

#' Generate the plateDatabase object
#' 
#' Generate the plateDatabase object (holding all neccessary information to
#' locate any plate on openBIS/in a local directory structure created by
#' BeeDataSetDownloader) using genome aggregate files.
#'
#' @return NULL, invisibly. The updated plateDatabase object is saved to the
#'         /data folder. For  the new file to be available, the package has to
#'         be reloaded.
#'
#' @examples
#' updateDatabasePlate()
#' # using the plateDatabase object
#' data(plateDatabase)
#' str(plate.database)
#' 
#' @export
updateDatabasePlate <- function() {
  
  processPathogen <- function(path) {
    
    processPlate <- function(barcode, data) {
      # extract all data corresponding to the current plate
      plate <- data[data$Barcode == barcode,]
      if(nrow(plate) != 384) {
        warning("incomplete plate (", barcode, "), expected 384 wells but ",
                "got ", nrow(plate), " wells.")
      }
      # get all values of interest and check if they are identical over the
      # whole plate
      space <- unique(plate$Space)
      if(length(space) != 1) {
        warning("different spaces within plate: ", barcode)
      }
      group <- unique(plate$Group)
      if(length(group) != 1) {
        warning("different groups within plate: ", barcode)
      }
      exper <- unique(plate$Experiment)
      if(length(exper) != 1) {
        warning("different experiments within plate: ", barcode)
      }
      return(c(barcode, space, group, exper))
    }
    
    message("processing ", tail(unlist(strsplit(path, "/")), n=1))
    # load the genome aggregate file of the current pathogen
    pathogen <- read.delim(path, as.is=TRUE)
    # find all plates of the current pathogen and
    plates <- unique(pathogen$Barcode)
    # process plates individually
    result <- sapply(plates, processPlate, pathogen)
    return(t(result))
  }
  
  # load path of genome aggregate files
  data(settingsDatabase, envir=environment())
  # find all genome aggregate files
  files <- list.files(path=settings.database$gen.aggr, pattern="\\.csv$",
                      full.names=TRUE)
  n.cores <- detectCores()
  registerDoMC(cores=n.cores)
  message("found ", length(files), " .csv files; using ", n.cores, " cores.")
  # process each of the files (they are per pathogen)
  plate.database <- foreach(i=1:length(files), .combine=rbind) %doparMC% {
    processPathogen(files[[i]])
  }
  message(paste(plate.database$out, collapse="\n"))
  plate.database <- plate.database$res
  colnames(plate.database) <- c("Barcode", "Space", "Group", "Experiment")
  plate.database <- as.data.frame(plate.database, stringsAsFactors=FALSE)
  rownames(plate.database) <- NULL
  save(plate.database,
       file=paste0(settings.database$package, "/data/plateDatabase.rda"),
       compression_level=1)
  invisible(NULL)
}

#' Generate the wellDatabase{Pathogen} objects
#' 
#' For each specified pathogen, generates the corresponding 
#' wellDatabase{Pathogen} object (holding well type, target gene ID and name for
#' each well) using genome aggregate files, as well as a kinome aggregate file
#' which contains additional information of control types which are not
#' currently available for genome plates.
#'
#' @param pathogens An optional list of pathogen names to have their
#'                  corresponding database updated. If not specified, all
#'                  pathogens are updated.
#'
#' @return NULL, invisibly. The updated plateDatabase object for each of the
#'         specified pathogens is saved to the /data folder. For the new files
#'         to be available, the package has to be reloaded.
#'
#' @examples
#' updateDatabaseWells(c("adeno", "brucella"))
#' # using the wellDatabase{Pathogen} objects
#' data(wellDatabaseAdeno)
#' str(well.database.adeno)
#' 
#' @export
updateDatabaseWells <- function(pathogens=NULL) {
  
  processPathogen <- function(path, supp.data) {
    message("processing ", tail(unlist(strsplit(path, "/")), n=1))
    # read genome aggregate file of the current pathogen
    pathogen.data.all <- read.delim(path, as.is=TRUE)
    # reduce the dataset a handful of cols (keep filesize/loading times down)
    pathogen.data.gen <- pathogen.data.all[c("Barcode", "WellRow", "WellColumn",
                                             "WellType", "ID", "Name")]
    # get the name of the current pathogen: "PATHOGEN_TEAM"
    pathogen.name <- unique(pathogen.data.all$Group)
    if(length(pathogen.name) != 1) stop("different group names within pathogen")
    # get the name of the current pathogen: strip "_TEAM"
    pathogen.name <- unlist(strsplit(pathogen.name, "_"))
    if(length(pathogen.name) != 2 | pathogen.name[2] != "TEAM") {
      stop("something is not right with the group name")
    }
    # get the name of the current pathogen: "pathogen"
    pathogen.name <- tolower(pathogen.name[1])
    # extract data for current pathogen from kinome aggregate
    pathogen.data.kin <- supp.data[grep(paste0("^", toupper(pathogen.name), 
                                               "-TEAM"), supp.data$Experiment),]
    # find all plates in the current genome dataset
    barcode.gen <- unique(pathogen.data.gen$Barcode)
    # find all plates in the current kinome dataset    
    barcode.kin <- unique(pathogen.data.kin$Barcode)
    # find all plates in both the current genome and kinome datasets      
    matches <- barcode.kin %in% barcode.gen
    barcode.upd <- barcode.kin[matches]
    message("augmenting ", length(barcode.upd), " of the ", length(barcode.gen),
            " found plates with kinome aggregate information.")
    # include all wells for each plate to be updated with kinome data
    update <- paste(rep(barcode.upd, each=384),
                    rep(rep(LETTERS[1:16], each=24), length(barcode.upd)),
                    rep(rep(1:24, 16), length(barcode.upd)), sep=":")
    # get indices of all wells to be updated
    ind.gen <- match(update, paste(pathogen.data.gen$Barcode,
                                   pathogen.data.gen$WellRow,
                                   pathogen.data.gen$WellColumn, sep=":"))
    # get indices of all wells containing the data for the update
    ind.kin <- match(update, paste(pathogen.data.kin$Barcode,
                                   pathogen.data.kin$WellRow,
                                   pathogen.data.kin$WellColumn, sep=":"))
    # update some of the cols
    pathogen.data.gen[ind.gen,]$Name <- pathogen.data.kin[ind.kin,]$GeneName
    pathogen.data.gen[ind.gen,]$ID   <- pathogen.data.kin[ind.kin,]$GeneID    
    # set the name the final object will have
    object.name   <- paste0("well.database.", pathogen.name)
    # set file name for the result    
    file.name     <- paste0(settings.database$package, "/data/wellDatabase", 
                            toupper(substring(pathogen.name, 1, 1)), 
                            substring(pathogen.name, 2),".rda")
    # assign to object name to the object (when the file is loaded later on, it
    # will have this name)
    assign(object.name, pathogen.data.gen)
    save(list=object.name, file=file.name, compression_level=1)
  }
  
  # load path of genome/kinome aggregate files
  data(settingsDatabase, envir=environment())
  # search for genome aggregate files  
  files <- list.files(path=settings.database$gen.aggr, pattern="\\.csv$",
                      full.names=TRUE)
  # if a list of pathogens ist specified, drop the other files
  if(!is.null(pathogens)) {
    matches <- unlist(lapply(pathogens, grep, files, ignore.case=TRUE))
    matches <- unique(matches)
    if(length(matches) < 1) {
      stop("no pathogens found matching your description.")      
    }
    files <- files[matches]
  }
  # search for kinome aggregate file
  aux <- list.files(path=settings.database$kin.aggr, pattern="\\.csv$",
                      full.names=TRUE)
  if(length(aux) != 1) stop("no supporting kinome data file found")
  # load kinome aggregate file  
  kin <- read.table(aux, header = TRUE, sep = ";", fill = TRUE,
                    stringsAsFactors = FALSE, comment.char = "")
  n.cores <- detectCores()
  registerDoMC(cores=n.cores)
  message("found 1 kinome .csv file:\n", aux, "\nand ", length(files), 
          " genome .csv files; using ", n.cores, " cores.")
  # process the data one pathogen at the time
  result <- foreach(i=1:length(files)) %doparMC% {
    processPathogen(files[[i]], kin)
  }
  message(paste(result$out, collapse="\n"))
  invisible(NULL)
}