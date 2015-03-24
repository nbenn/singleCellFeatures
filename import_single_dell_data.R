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
    name <- tail(unlist(strsplit(filepath, "[/]")),n=1)
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
  
  cat("fetching ", length(filenames), " files.\n", sep = "")
  
  # get the data
  import.data <- lapply(filenames, readFeatureFile)
  import.name <- lapply(filenames, getFeatureName)
  names(import.data) <- import.name
  
  cat("successfully read features:\n ", paste(import.name, collapse="\n  "))
  
  return(import.data)
}

## Loading data

data <- importPlate("/Users/nbennett/Polybox/MasterThesis/openBISDownload/INFECTX_PUBLISHED/ADENO_TEAM/ADENO-AU-K1/KB2-02-1I", c("cells", "AreaShape"))
names(data)
