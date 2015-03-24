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
  #
  # TODO(nbenn): Add metadata to lowest list hierarchy, flatten hierarchy by
  #              one level
  
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

extractImageData <- function(data, row, col, img=5) {
  # Given feature data from a whole plate, extract the data for a single image
  # in a specified well
  #
  # Args:
  #   data: Data for a whole plate, list of features
  #   row:  Row of the well, specified either by A:P, a:P, or 1:16
  #   col:  Column number of the well, specified by 1:24
  #   img:  Image number within the well, 1:9
  #
  # Returns:
  #   A list of features with all values for the given image
  
  # Input validation
  # 9 entries per well, 16x24 wells
  if(any(sapply(data, function(x) length(x) != 3456))) 
    stop("expecting data to be of length 3456")
  # 9 images per well
  img <- as.integer(img)
  if(img < 1 | img > 9) stop("img has to be in the range 1:9")
  # 24 columns per plate
  col <- as.integer(col)
  if(col < 1 | col > 24) stop("col has to be in the range 1:24")
  # 16 rows indicated by A/a:P/p
  if(!is.integer(row)) {
    if (row == "a" | row == "A") row <- 1
    else if (row == "b" | row == "B") row <- 2
    else if (row == "c" | row == "C") row <- 3
    else if (row == "d" | row == "D") row <- 4
    else if (row == "e" | row == "E") row <- 5
    else if (row == "f" | row == "F") row <- 6
    else if (row == "g" | row == "G") row <- 7
    else if (row == "h" | row == "H") row <- 8
    else if (row == "i" | row == "I") row <- 9
    else if (row == "j" | row == "J") row <- 10
    else if (row == "k" | row == "K") row <- 11
    else if (row == "l" | row == "L") row <- 12
    else if (row == "m" | row == "M") row <- 13  
    else if (row == "n" | row == "N") row <- 14 
    else if (row == "o" | row == "O") row <- 15 
    else if (row == "p" | row == "P") row <- 16
    else row <- as.integer(row)
  }
  # or 1:16
  if(row < 1 | row > 16) stop("row has to be in the range 1:16, a:p or A:P")
  
  # the data is in row major order, in groups of 9 images
  index <- ((row-1)*24+(col-1))*9+img

  # extract data of index for all features
  return(lapply(data, function(x) x[[index]][[1]]))
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

KB2_02_1I.f12.7 <- extractImageData(data=dat.plat, row="F", col=12, img=7)
