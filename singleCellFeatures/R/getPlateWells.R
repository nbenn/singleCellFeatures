#' Fetch plate wells specified by experiment, barcode, row and column
#' 
#' Given results from findPlateWellsFromGene, or a table containing columns for
#' experiment, barcode, well row and well column, fetch all needed plates (from
#' openBIS if needed, or from matlab files if present) and extract the specified
#' wells (result are saved such that plates do not have to be re-imported later 
#' on)
#'
#' @param coordinates    Table holding locations of wells (columns for
#'                       experiment, barcode, well row and well column)
#' @param force.download Force downloading from openBIS
#' @param force.import   Force importing from matlab files
#' @param keep.files     Keep matlab files after successful import and 
#'                       well data extraction
#'
#' @return A nested list with the top level hierarchy corresponding to a plate 
#'         and the level below to a well and for structure below, see 
#'         importCompletePlate
#'
#' @examples
#' # get gene locations
#' mtor.bruc.du <- findPlateWellsFromGene("brucella", "du-k", "MTOR")
#' scra.bruc.du <- findPlateWellsFromGene("brucella", "du-k", "SCRAMBLED")
#' # combine for faster fetching
#' targets <- rbind(mtor.bruc.du[1:6,], scra.bruc.du[c(6,18,78,90,150,162),])
#' 
#' data <- getPlateWells(targets)
#' 
#' @export

getPlateWells <- function(coordinates, force.download=FALSE, force.import=FALSE,
                          keep.files=FALSE) {
  # function to process several wells on a single plate
  processPlate <- function(plate, coords, force.download, force.import, 
                           keep.files) {
    # load file containing some needed paths
    data(paths, envir = environment())
    # find all wells on the current plate
    current.coords <- coords[coords$Barcode==plate,]
    experiment <- unique(current.coords$Experiment)
    barcode <- unique(current.coords$Barcode)
    
    if(length(experiment) != 1 | length(barcode) != 1) {
      stop("non unique barcode/experiment within single plate")
    }
    
    pathogen  <- toupper(unlist(strsplit(experiment, "-"))[1])
    # plate location in openBIS
    plateLoca <- paste0("/INFECTX_PUBLISHED/", toupper(pathogen), "_TEAM/", 
                        experiment, "/", barcode)
    # local plate path
    platePath <- paste0(baseDir, plateLoca)
    
    # search for previous rds files (one for each well)    
    foundRDS <- apply(current.coords, 1, function(coordinate) {
      prevRDS.name <- paste0(coordinate["Barcode"], "_", coordinate["WellRow"], 
                             coordinate["WellColumn"], ".rds")
      prevRDS.name <- gsub(" ", "", prevRDS.name, fixed = TRUE)
      prevRDS.file <- list.files(path=platePath, pattern=prevRDS.name, 
                                 full.names=TRUE, recursive=TRUE)
      if (length(prevRDS.file) == 1 & !force.download & !force.import) {
        return(TRUE)
      } else return(FALSE)
    })
    
    # if for each well, an rds file was found, load the all
    if(all(foundRDS)) {
      cat("\nfound previous results stored in .rds files.\n")
      result <- apply(current.coords, 1, function(coordinate) {
        # search for previous rds file  
        prevRDS.name <- paste0(coordinate["Barcode"], "_",
                               coordinate["WellRow"], coordinate["WellColumn"],
                               ".rds")
        prevRDS.name <- gsub(" ", "", prevRDS.name, fixed = TRUE)
        prevRDS.file <- list.files(path=platePath, pattern=prevRDS.name, 
                                   full.names=TRUE, recursive=TRUE)
        return(readRDS(prevRDS.file[1]))
      })
      names(result) <- paste0(current.coords$WellRow, current.coords$WellColumn)
    } else {
      # if not all needed wells are already storen in rds files, search for
      # matlab files
      matlab.files <- list.files(path=platePath, pattern="\\.mat$", 
                                 full.names=TRUE, recursive=TRUE)
      # if no matlab files are found, download the whole plate from openBIS
      if (length(matlab.files) == 0 | force.download) {
        cat("\ndownloading plate from openBIS.\n")
        # get password from file
        pwd.openBIS <- readRDS(filepaths$openBIS.pwd)
        fetchCommand <- paste0(
          "cd ", filepaths$bee.dir, "; export BEESOFTSRC='", 
          filepaths$bee.softsrc, "'; ./BeeDataSetDownloader.sh --user ",
          "'nbennett@student.ethz.ch' --password '", pwd.openBIS,
          "' --outputdir '", baseDir, "' --plateid '^", plateLoca,
          "' --files '.*.mat'", " --verbose '10'")
        system(fetchCommand)
      } else {
        cat("\nusing the previously dowloaded .mat files.\n")
      }
      # now the whole plate is present in .mat files -> import most of it
      data <- importCompletePlate(platePath, 
                                  features=c("^Cells.", "^Nuclei.", 
                                             "^PeriNuclei.", "^VoronoiCells."))
      # get linear indices of all wells 
      wells <- mapply(getLinearizedPlateIndex, current.coords$WellRow, 
                      current.coords$WellColumn)
      # extract wells of interest from whole plate
      result <- extractPartialPlate(data, wels=wells, imgs=NULL)
      # save rds files per well for faster reloading
      lapply(result$data, function(x) {
        name.RDS <- paste0(barcode, "_", x$meta$WellRow, x$meta$WellColumn,
                           ".rds")
        saveRDS(x, paste(platePath, name.RDS, sep="/"))
      })
    }
    # unless keep.files=TRUE, all matlab files of the current plate are deleted
    if(!keep.files) {
      matlab.files <- list.files(path=platePath, pattern="\\.mat$",
                                 full.names=TRUE, recursive=TRUE)
      if(length(matlab.files) > 0) {
        cat("\n\ndeleting downloaded .mat files.\n\n")    
        iterator <- 1
        repeat {
          target1 <- unlist(strsplit(matlab.files[iterator], 
                                     paste0("/", barcode, "/")))
          target2 <- unlist(strsplit(target1[2], "/"))
          if (length(target2 == 3)) break
          iterator <- iterator + 1  
        }
        target <- paste(target1[1], barcode, target2[1], sep="/")
        unlink(target, recursive=TRUE) 
      }
    }
    if(all(foundRDS)) return(result)
    else return(result$data)
  }
  
  # load file containing some needed paths
  data(paths, envir = environment())
  baseDir   <- filepaths$openBIS.data
  # find all plates that occur in input
  plates <- unique(coordinates$Barcode)
  # process plate by plate (interleaving plates is slow)
  result <- lapply(plates, processPlate, coordinates, force.download, 
                   force.import, keep.files)
  names(result) <- plates
  return(result)
}