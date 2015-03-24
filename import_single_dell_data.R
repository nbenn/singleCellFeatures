## Download data
# cd /Users/nbennett/ETH/MasterThesis/software/tginfectx/trunk/openBIS/Tools/
#   BeeDataSetDownloader
# export BEESOFTSRC="/Users/nbennett/ETH/MasterThesis/software/tginfectx/trunk/
#   screeningBee"
# ./BeeDataSetDownloader.sh 
#   --user "nbennett@student.ethz.ch"
#   --password "****"
#   --outputdir "/Users/nbennett/Polybox/MasterThesis/openBISDownload"
#   --plateid "^/INFECTX_PUBLISHED/ADENO_TEAM/ADENO-AU-K1/KB2-03-1I"
#   --files ".*(Intensity).*.mat"
#   --verbose "10"

## Functions
importPlate <- function(path, features=NULL) {
  # Read all available single cell feature data for a given plate
  #
  # Args:
  #   path:     Folder path to folder corresponding to a single plate, holding
  #             single cell feature data in .mat files somewhere below
  #   features: Vector holding strings corresponding to features to be imported
  #
  # Returns:
  #   A nested list with the top level hierarchy corresponding to a feature and
  #   the next lower level corresponding to the data in a single image (there 
  #   are 3456 = 9 images x 24 rows x 36 cols slots on this level)
  
  readFeatureFile <- function(filepath) {
    require("R.matlab")
    # nesting of list contains no information
    data <- readMat(filepath)[[1]][[1]][[1]][[1]]
    return(data)
  }
  
  getFeatureName <- function(filepath) {
    # extract the filename from the path
    name <- tail(unlist(strsplit(filepath, "[/]")), n=1)
    # split into substrings separated by "."
    name <- unlist(strsplit(name, "[.]"))
    # remove last substring, hopefully ".mat"
    name <- name[-length(name)]
    name <- paste(name, collapse=".")
    return(name)
  }
  
  # list all files ending in .mat below the input dir
  filenames <- list.files(path = path, pattern="\\.mat$", full.names=TRUE,
                          recursive=TRUE)
  
  # if features specified, drop other files
  if (!is.null(features)) {
    # for each feature entry, get all partial matches
    keep <- lapply(features, grep, x=filenames, ignore.case=TRUE)
    # remove duplicates from match index
    keep <- unique(unlist(keep))
    # drop unmatched features
    filenames <- filenames[keep]    
  }
  
  cat("\nfetching ", length(filenames), " files.\n", sep = "")
  
  # get the data
  import.data <- lapply(filenames, readFeatureFile)
  import.name <- lapply(filenames, getFeatureName)
  names(import.data) <- import.name
  
  cat("successfully read features:\n ", paste(import.name, collapse="\n  "))
  
  return(import.data)
}

importExperiment <- function(path, plates, features) {
  # Read all available single cell feature data for a given experiment,
  # spanning multiple plates
  #
  # Args:
  #   path:     Folder path to folder corresponding to an experiment
  #   plates:   Names of plates belonging to the given experiment that are to
  #             be imported
  #   features: Vector holding strings corresponding to features to be imported
  #
  # Returns:
  #   A nested list with the top level hierarchy corresponding to a plate and
  #   the level below to a feature
  
  # get all directories below the specified path
  all <- list.dirs(path=path, full.names=TRUE, recursive=TRUE)
  # array (a col for each plate and a row for each dir in all) of bool,
  # specifying the dirs that match each plate (one TRUE per col) 
  matches <- sapply(plates, function(pattern, paths) {
    sapply(paths, function(x) {
      grepl(pattern, tail(unlist(strsplit(x, "[/]")), n=1))
    })
  }, all)
  # collapse array to vector using OR operation
  matches <- apply(matches, 1, any)
  # drop unused dirs
  plate.paths <- all[matches]
  
  # import data
  import.data <- lapply(plate.paths, importPlate, features)
  import.name <- lapply(plate.paths, function(x) {
    name <- paste(tail(unlist(strsplit(x, "[/]")), n=2), collapse=".")
    return(gsub("-", "_", name))
  })
  names(import.data) <- import.name

  return(import.data)
}

## Loading data

# import a single plate
path     <- paste("/Users/nbennett/Polybox/MasterThesis/openBISDownload/",
                  "INFECTX_PUBLISHED/ADENO_TEAM/ADENO-AU-K1/KB2-02-1I", sep="")
features <- c("cells", "AreaShape")
dat.plat <- importPlate(path, features)

# import whole experiment
path     <- paste("/Users/nbennett/Polybox/MasterThesis/openBISDownload/",
                  "INFECTX_PUBLISHED/ADENO_TEAM", sep="")
plates   <- c("KB2-03-1I", "KB2-03-1J")
features <- c("Cells.AreaShape_Area", "Cells.AreaShape_Eccentricity")
dat.expe <- importExperiment(path, plates, features)

