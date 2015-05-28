#' Update the featureDatabase object
#' 
#' In order to update the featureDatabase object, all *.txt files in the folder 
#' /path/to/meta.dir/FeatureLists are read, processed and the feature
#' "Cells.Infection_IsInfected" is added. The filenames determine the list slot
#' names, so they should be formatted as PATHOGEN-*.txt.
#'
#' @return NULL (invisibly). The updated featureDatabase object saved to the
#'         /data folder. For the new file to be available, the package has to
#'         be reloaded.
#'
#' @examples
#' updateDatabaseFeatures()
#' 
#' @export
updateDatabaseFeatures <- function() {
  # load config (for data path)
  config <- configGet()
  dir <- paste0(config$dataStorage$metaDir, "/", "FeatureLists")
  files <- list.files(dir, pattern="\\.txt$", full.names=TRUE)
  if(length(files) != 10) warning("expecting 10 files, got ", length(files))
  filenames <- basename(files)
  filenames <- sapply(filenames, function(x) {
    tolower(unlist(strsplit(unlist(strsplit(x, "[.]"))[1], "-"))[1])
  })
  feature.database <- lapply(files, read.delim, header=FALSE, sep=" ",
                     stringsAsFactors=FALSE)
  names(feature.database) <- filenames
  feature.database <- lapply(feature.database, function(x) {
    feat <- x[,1]
    feat <- sapply(feat, function(f) {
      if(length(grep(".mat$", f)) != 1) {
        warning("feature ", f, "does not end in .mat")
      }
      else f <- unlist(strsplit(f, ".mat"))
    })
    names(feat) <- NULL
    feat <- c(feat, "Cells.Infection_IsInfected")
    return(feat)
  })

  lengths <- sapply(feature.database, function(x) length(x))
  message("successfully read ", length(files), " feature lists:\n", 
          paste("  using", lengths, "features for", filenames,
                collapse="\n", sep=" "))

  save(feature.database,
       file=paste0(config$singleCellFeatures$sourceDir, "/",
                   "data/featureDatabase.rda"),
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

  # load HCS_ANALYSIS_CELL_FEATURES_CC_MAT.tsv file
  config <- configGet()
  filename <- paste0(config$dataStorage$metaDir, "/",
    "HCS_ANALYSIS_CELL_FEATURES_CC_MAT.tsv")
  all <- read.delim(filename, stringsAsFactors = FALSE)
  # some datasets occur multiple times per plate
  all.barco <- unique(all$Sample)
  # for each unique plate barcode
  plate.database <- sapply(
    all.barco,
    function(x, dat) {
      hits <- dat[dat$Sample == x, ]
      # if multiple datasets present, select most recent one for DataID field
      if(nrow(hits) > 1) {
        hits <- hits[order(hits$Code, decreasing=TRUE),][1,]
      }
      barcode    <- as.character(x)
      group      <- as.character(hits["Project"])
      experiment <- as.character(hits["Experiment"])
      data.id    <- as.character(hits["PermID"])
      space      <- unlist(strsplit(as.character(hits["Experiment.Identifier"]),
                                    "/"))
      if(!all(space[1] == "" & space[3] == group & space[4] == experiment))
        stop("unexpected formating of Experiment.Identifier column")
      return(c(barcode, space[2], group, experiment, data.id))
    },
    all)
  plate.database <- t(plate.database)
  plate.database <- as.data.frame(plate.database, stringsAsFactors=FALSE)
  colnames(plate.database) <- c("Barcode", "Space", "Group", "Experiment",
    "DataID")
  plate.database <- plate.database[order(plate.database$Experiment),]
  rownames(plate.database) <- NULL

  save(plate.database,
       file=paste0(config$singleCellFeatures$sourceDir, "/",
                   "data/plateDatabase.rda"),
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

  processPathogen <- function(path) {
    # read genome aggregate file of the current pathogen
    pathogen.all <- read.delim(path, stringsAsFactors=FALSE)
    # reduce the dataset a handful of cols (keep filesize/loading times down)
    pathogen.red <- pathogen.all[c("Barcode", "WellRow", "WellColumn",
                                   "WellType", "ID_openBIS", "ID",
                                   "Name")]
    colnames(pathogen.red) <- c("barcode", "well.row", "well.col", "well.type",
                                "id.openBIS", "id.infx",  "name")
    # get the name of the current pathogen: "PATHOGEN_TEAM"
    pathogen.name <- unique(pathogen.all$Group)
    if(length(pathogen.name) != 1) stop("different group names within pathogen")
    if(pathogen.name != "MOCK") {
      # get the name of the current pathogen: strip "_TEAM"
      pathogen.name <- unlist(strsplit(pathogen.name, "_"))
      if(length(pathogen.name) != 2 | pathogen.name[2] != "TEAM") {
        stop("something is not right with the group name")
      }
    }
    # get the name of the current pathogen: "pathogen"
    pathogen.name <- tolower(pathogen.name[1])
    # set the name the final object will have
    object.name   <- paste0("well.database.", pathogen.name)
    # set file name for the result
    file.name     <- paste0(config$singleCellFeatures$sourceDir,
                            "/", "data/wellDatabase",
                            toupper(substring(pathogen.name, 1, 1)),
                            substring(pathogen.name, 2),".rda")
    # assign to object name to the object (when the file is loaded later on, it
    # will have this name)
    assign(object.name, pathogen.red)
    save(list=object.name, file=file.name, compression_level=1)
  }

  # load path of genome/kinome aggregate files
  config <- configGet()
  # search for genome aggregate files
  dir <- paste0(config$dataStorage$metaDir, "/", "Aggregates")
  files <- list.files(path=dir, pattern="\\.csv$", full.names=TRUE)
  # if a list of pathogens ist specified, drop the other files
  if(!is.null(pathogens)) {
    matches <- unlist(lapply(pathogens, grep, files, ignore.case=TRUE))
    matches <- unique(matches)
    if(length(matches) < 1) {
      stop("no pathogens found matching your description.")
    }
    files <- files[matches]
  }
  #registerDoParallel(cores=detectCores())
  #message("found ", length(files), " aggregate files; using ",
  #        getDoParWorkers(), " cores.")
  message("found ", length(files), " aggregate files.")
  # process the data one pathogen at the time
  l_ply(files, processPathogen,
        .progress=getOption("singleCellFeatures.progressBars"))
  invisible(NULL)
}

#' Assess the coverage provided by available metadata
#' 
#' The plate database is compared to all well databases in order to find all 
#' plates that are missing in well databases. The plate database is assumed
#' to be complete in the sense that it is generated through a data set search
#' in openBIS and should therefore include all plates for which single cell
#' data is available. The well databases however are extracted from aggregates
#' which only become available periodically.
#'
#' @param verbose    A verbosity argument, specifying whether to print all
#'                   plates that are missing of just a summary.
#' @param is.startup Only relevant when this function is used upon package
#'                   loading to output via packageStartupMessage.
#'
#' @return NULL, invisibly. All information is printed.
#'
#' @examples
#' wellDatabaseCoverage()
#' wellDatabaseCoverage(TRUE)
#' 
#' @export
wellDatabaseCoverage <- function(verbose=FALSE, is.startup=FALSE) {
  data(plateDatabase, envir=environment())
  groups      <- unique(plate.database$Group)
  experiments <- lapply(groups, function(group, dat) {
    hits <- dat[dat$Group == group,]
    return(unique(hits$Experiment))
  }, plate.database)
  names(experiments) <- groups
  
  well.data <- lapply(groups, function(group) {
    name <- tolower(unlist(strsplit(group, "_"))[1])
    name <- paste0(toupper(substring(name, 1, 1)),
      substring(name, 2))
    dataset.name <- paste0("wellDatabase", name)
    object.name  <- paste0("well.database.", tolower(name))
    suppressWarnings(
      well.db <- try({
        data(list=dataset.name, envir=environment())
        get(object.name)
      },
      silent = TRUE)
    )
    if(class(well.db) == "try-error") {
      if(is.startup) {
        packageStartupMessage("for ", name, ", no well database found")
      } else {
        warning("for ", name, ", no well database found", call.=FALSE)
      }
      return(NULL)
    } else {
      barcodes <- unique(well.db$barcode)
      n.duplic <- data.frame(table(well.db$barcode))
      duplic <- n.duplic[n.duplic$Freq != 384,]
      if(nrow(duplic) > 0) {
        if(is.startup) {
          packageStartupMessage("for ", name, ", some metadata might be ",
                                "incomplete for plates:\n",
                                paste("  ", duplic[,1], ": ", duplic[,2],
                                      " wells", sep="", collapse="\n"))
        } else {
          warning("for ", name, ", some metadata might be incomplete for ",
                  "plates:\n", paste("  ", duplic[,1], ": ", duplic[,2],
                                     " wells", sep="", collapse="\n"),
                  call.=FALSE)
        }
      }
      return(barcodes)
    }
  })
  names(well.data) <- groups

  missing <- lapply(groups, function(group, experiments, plate.db, well.db) {
    res <- lapply(experiments[[group]], function(experiment, plate, well) {
      res <- setdiff(plate[plate$Experiment == experiment, ]$Barcode, well)
      if(length(res) > 0) return(res)
      else return(NULL)
    }, plate.db[plate.db$Group == group,], well.db[[group]])
    names(res) <- experiments[[group]]
    res <- res[!sapply(res, is.null)]
    return(res)
  }, experiments, plate.database, well.data)

  names(missing) <- groups
  if(is.startup) {
    packageStartupMessage("missing metadata for ", length(unlist(missing)),
                          " plates.\ncoverage: ",
                          1 - length(unlist(missing)) / nrow(plate.database))
  } else {
    message("missing metadata for ", length(unlist(missing)), " plates.\n",
            "coverage: ", 1 - length(unlist(missing)) / nrow(plate.database))    
  }
  if(verbose) {
    drop <- lapply(groups, function(group, miss) {
      plates      <- miss[[group]]
      experiments <- names(plates)
      message("for group ", group, ": missing ", length(unlist(plates)),
              " plates.")
      drop <- lapply(experiments, function(experiment, barcodes) {
        bc <- plates[[experiment]]
        message("  for experiment ", experiment, ", missing ", length(bc),
                " plates.")
        bc <- paste("\"", bc, "\"", collapse="  ", sep="")
        colwidth <- max(nchar(bc))
        bc.pad <- stri_pad_right(bc, colwidth)
        str <- paste(bc.pad, collapse="  ", sep="")
        str <- stri_wrap(str, normalize=FALSE)
        message("    ", paste(str, collapse="\n    "))
      }, plates)
    }, missing)
  } else {
    if(is.startup) {
      packageStartupMessage("run wellDatabaseCoverage(TRUE) to show all ",
                            "missing plates.")
    } else {
      message("run wellDatabaseCoverage(TRUE) to show all missing plates.")
    }
  }

  invisible(NULL)
}