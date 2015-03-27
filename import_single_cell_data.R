## Download data
# cd /Users/nbennett/ETH/MasterThesis/software/tginfectx/trunk/openBIS/Tools/
#   BeeDataSetDownloader
# export BEESOFTSRC="/Users/nbennett/ETH/MasterThesis/software/tginfectx/trunk"
# ./BeeDataSetDownloader.sh 
#   --user "nbennett@student.ethz.ch"
#   --password "****"
#   --outputdir "/Users/nbennett/Polybox/MasterThesis/openBISDownload"
#   --plateid "^/INFECTX_PUBLISHED/ADENO_TEAM/ADENO-AU-K1/KB2-03-1I"
#   --files ".*(Intensity).*.mat"
#   --verbose "10"

## Functions
getLinearizedPlateIndex <- function(row, col, img=NULL) {
  # Given row, column and possible image coordinates on a plate, get the
  # corresponding index of a linearized representation
  #
  # Args:
  #   row:  Row of the well, specified either by A:P, a:P, or 1:16
  #   col:  Column number of the well, specified by 1:24
  #   img:  NULL/Image number within the well, 1:9
  #
  # Returns:
  #   The linearized index; if img==NULL, range 1:384, else range: 1:3456
  
  # Input validation
  col <- as.integer(col)
  # 24 columns per plate  
  if (col < 1 | col > 24) stop("col has to be in the range 1:24")
  # 16 rows indicated by A/a:P/p
  if (!is.integer(row)) {
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
  if (row < 1 | row > 16) stop("row has to be in the range 1:16, a:p or A:P")
  
  if (!is.null(img)) {
    # 9 images per well
    img <- as.integer(img)
    if(img < 1 | img > 9) stop("img has to be in the range 1:9")
    # the data is in row major order, in groups of 9 images
    index <- ((row-1)*24+(col-1))*9+img
  } else {
    # the data is in row major order
    index <- (row-1)*24+col
  }
  
  return(index)    
}

importCompletePlate <- function(path, features=NULL) {
  # Read all available/the specified single cell feature data for a plate given
  # by its path
  #
  # Args:
  #   path:     Folder path to folder corresponding to a single plate, holding
  #             single cell feature data in .mat files somewhere below
  #   features: Vector holding strings corresponding to features to be imported
  #
  # Returns:
  #   A nested list with the following structure:
  #   --> plate
  #     |--> metadata for plate (14 fields)
  #     |--> data (list of length 384; 24 cols by 16 rows)
  #        |--> A1
  #           |--> ...
  #       ...
  #        |--> P24
  #           |--> metadata for well (40 fields)
  #           |--> data (list of length 9; 9 images per well)
  #              |--> img_11
  #                 |--> ...
  #             ...
  #              |--> img_33
  #                 |--> Feature1
  #                ...
  #                 |--> FeatureN

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
  
  getMetadata <- function(plate.path) {
    # figure out plate name
    plate.name <- tail(unlist(strsplit(plate.path, "[/]")), n=1)
    previous.filename <- paste(plate.path, paste(plate.name, "metadata.rds", 
                                                 sep="_"), sep="/")
    # try to get metadata from previous run
    if (file.exists(previous.filename)) {
      cat("\nfetching metadata from previous run.\n", sep = "")              
      result <- readRDS(previous.filename)
    } else {
      # if no previous results are available, fetch from aggregate
      cat("\nfetching metadata from aggregate.\n", sep = "")        
      metad.path <- paste(plate.path, "../../../../../Aggregates", sep="/")
      # figure out experiment name
      exper.name <- tail(unlist(strsplit(plate.path, "[/]")), n=2)[1]
      # figure out pathogen name
      patho.name <- unlist(strsplit(exper.name, "[-]"))[1]
      metad.file <- list.files(path=metad.path, pattern=patho.name,
                               ignore.case=TRUE, full.names=TRUE)
      if (length(metad.file) != 1) {
        stop("found ", length(metad.file), " metadata files instead of 1.")
      }    
      # load metadata file
      metad.all  <- read.delim(metad.file, as.is=TRUE)
      # drop all except the plate of interest
      metad.plat <- metad.all[metad.all$Barcode == plate.name,]
      if (nrow(metad.plat) != 384) {
        stop("found information on ", nrow(metad.plat), 
             " wells instead of 384.")}
      # the first 14 cols are concerned with per-plate metadata
      plate <- metad.plat[,1:14]
      # per plate metadata is reduced to one entry per col
      plate <- apply(plate, 2, unique)
      if (length(plate) != 14) {
        print(plate)
        stop("non-identical per plate metadata")
      }
      # the rest of the metadata is per well
      well <- metad.plat[,-(1:14)]
      # add the linearized index for easier matching later on
      well$SerialIndex <- apply(well, 1, function(x) {
        row <- x["WellRow"]
        col <- x["WellColumn"]
        return(getLinearizedPlateIndex(row, col, NULL))
      })
      # save & return a list of per-plate and per-well metadata
      result <- list(plate=plate, well=well)
      saveRDS(result, previous.filename) 
    }
    return(result)
  }
  
  # list all files ending in .mat below the input dir
  filenames <- list.files(path = path, pattern="\\.mat$", full.names=TRUE,
                          recursive=TRUE)
  # get the metadata file
  metadata <- getMetadata(path)
  
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
  
  # reorder list nesting to match the specified target
  result <- lapply(1:384, function(index.plate, data, meta) {
    images <- lapply(1:9, function(index.img, index.plate, data) {
      lapply(data, function(x, i, j) {
        return(x[[((i-1)*9+j)]][[1]])
      }, index.plate, index.img)
    }, index.plate, data)
    names(images) <- c("img_11", "img_12", "img_13", 
                       "img_21", "img_22", "img_23", 
                       "img_31", "img_32", "img_33")
    return(list(meta=meta[meta$SerialIndex==index.plate,], data=images))
  }, import.data, metadata$well)
  well.names <- unlist(lapply(c("A", "B", "C", "D", "E", "F", "G", "H", "I", 
                                "J", "K", "L", "M", "N", "O", "P"), 
                              paste, 1:24, sep=""))
  names(result) <- well.names
  
  return(list(meta=metadata$plate, data=result))
}

importMultiplePlates <- function(path, plates=NULL, features=NULL) {
  # Read all available/the specified set of single cell feature data for all 
  # plates/the specified set below the dir path
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
  # if no plates specified, find all folders below path that correspond to a
  # plate
  if (is.null(plates)) {
    # get pathlengths of all found folders 
    folder.depth <- sapply(all, function(x) length(unlist(strsplit(x, "[/]"))))
    # exclude all but the folders with longest paths (the data folders)
    data.folders <- all[ifelse(folder.depth == max(folder.depth), TRUE, FALSE)]
    # drop the last two hierarchies to get from data folders to plate folders
    plate.paths <- sapply(data.folders, function(x) {
      paste(unlist(strsplit(x, "[/]"))[1:(max(folder.depth)-2)], collapse="/")
    })
  } else {
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
  }
  
  # import data
  import.data <- lapply(plate.paths, importCompletePlate, features)
  import.name <- lapply(plate.paths, function(x) {
    name <- tail(unlist(strsplit(x, "[/]")), n=1)
    return(gsub("-", "_", name))
  })
  names(import.data) <- import.name

  return(import.data)
}

extractPartialPlate <- function(data, rows=NULL, cols=NULL, wels=NULL, 
                                imgs=5) {
  # Given feature data from a whole plate, extract the data for specified rows/
  # cols/wells/images
  #
  # Args:
  #   data: Data for a whole plate, multiply nested list of features
  #   rows: Row(s) of the well(s), specified either by a vector consisting of 
  #         the characters A:P, a:P, or 1:16
  #   cols: Column(s) of the well(s), specified by a vector of integers in the
  #         range of 1:24
  #   wels: Well index(es), specified by a vector of integers in the range of
  #         1:384
  #   imgs: Image number(s) within the well, specified by a vector of integers 
  #         in the range of 1:9
  #
  # Returns:
  #   A multiply nested list of features with all values for the given plate
  #   coordinates (for structure, refer to importCompletePlate)
  
  # Input validation
  # 16x24 wells: check input data, as well as vector for wells
  if (length(data$data) != 384) stop("expecting $data to be of length 384")
  if (!(is.vector(wels, mode="numeric") | is.null(wels))) {
    stop("expecting NULL or an integer vector for wels")
  }
  if (!is.null(wels)) wels <- as.integer(wels)
  if (any(ifelse(wels < 1 | wels > 384, TRUE, FALSE))) {
    stop("expecting 1 <= wels <= 384")
  }
  # 9 images per well
  if (any(sapply(data$data, function(x) length(x$data) != 9))) {
    stop("expecting $data$plate$data to be of length 9")
  }
  if (!(is.vector(imgs, mode="numeric") | is.null(imgs))) {
    stop("expecting NULL or an integer vector for imgs")
  }
  if (!is.null(imgs)) imgs <- as.integer(imgs)
  if (any(ifelse(imgs < 1 | imgs > 9, TRUE, FALSE))) {    
    stop("expecting 1 <= imgs <= 9")
  }
  # type/range of cols/rows will be checked in getLinearizedPlateIndex
  if (!(is.vector(rows) | is.null(rows))) {
    stop("expecting NULL or an integer vector for rows")
  }
  if (!(is.vector(cols) | is.null(cols))) {
    stop("expecting NULL or an integer vector for cols")
  }
  
  # for the specified rows and cols, calculate the linear well indices
  if (is.null(rows) & !is.null(cols)) {
    wells.extra <- as.vector(sapply(cols, function(i) {
      sapply(1:16, function(j, i) getLinearizedPlateIndex(j, i, NULL), i)
    }))
    wels <- c(wels, wells.extra)
  } else if (!is.null(rows) & is.null(cols)) {
    wells.extra <- as.vector(sapply(rows, function(i) {
      sapply(1:24, function(j, i) getLinearizedPlateIndex(i, j, NULL), i)
    }))
    wels <- c(wels, wells.extra)
  } else if (!is.null(rows) & !is.null(cols)) {
    wells.extra <- as.vector(sapply(rows, function(i, cols) {
      sapply(cols, function(j, i) getLinearizedPlateIndex(i, j, NULL), i)
    }, cols))
    wels <- c(wels, wells.extra)
  } else if (is.null(rows) & is.null(cols) & is.null(wels)) {
    wels <- 1:384
  }
  # drop duplicates that can occur through specifying the same well twice
  wels <- sort(unique(wels))
  
  # extract the data from input
  if (is.null(imgs)) {
    well.data <- lapply(wels, function (i, data) data$data[[i]], data)
  } else {
    imgs <- sort(unique(imgs))
    well.data <- lapply(wels, function (i, j, data) {
      images <- lapply(j, function(j, data) data$data[[j]], data$data[[i]])
      names(images) <- names(data$data[[i]]$data)[j]
      return(list(meta=data$data[[i]]$meta, data=images))
    }, imgs, data) 
  }
  names(well.data) <- names(data$data)[wels]
    
  return(list(meta=data$meta, data=well.data))
}

## Loading data
# import a single plate
path     <- paste("/Users/nbennett/Polybox/MasterThesis/openBISDownload/",
                  "INFECTX_PUBLISHED/ADENO_TEAM/ADENO-AU-K2/KB2-03-1J", sep="")
features <- c("location")
dat.plat <- importCompletePlate(path, features)

sum(sapply(dat.plat$data$B10$data, function(x) length(x[[1]])))

subset.img <- extractPartialPlate(dat.plat)
subset.wel <- extractPartialPlate(dat.plat, cols=c(3,4), rows=c("A","B"), 
                                  wels=c(333, 21), imgs=NULL)

# import whole experiment
path     <- paste("/Users/nbennett/Polybox/MasterThesis/openBISDownload/",
                  "INFECTX_PUBLISHED/ADENO_TEAM", sep="")
plates   <- c("KB2-03-1I", "KB2-03-1J")
features <- c("Cells.AreaShape_Area", "Cells.AreaShape_Eccentricity")
dat.expe <- importMultiplePlates(path, NULL, features)

