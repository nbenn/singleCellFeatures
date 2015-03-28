#' Read specified single cell features for a single plate
#'
#' Read all available/the specified single cell feature data for a plate given
#' by its path
#'
#' @param path     Folder path to folder corresponding to a single plate, 
#'                 holding single cell feature data in .mat files somewhere 
#'                 below
#' @param features Vector holding strings corresponding to features to be 
#'                 imported
#'
#' @return A nested list; for structure see examples section
#'
#' @examples
#' path     <- paste("/Users/nbennett/Polybox/MasterThesis/openBISDownload/",
#'                   "INFECTX_PUBLISHED/ADENO_TEAM/ADENO-AU-K2/KB2-03-1J", 
#'                   sep="")
#' features <- c("location")
#' dat.plat <- importCompletePlate(path, features)
#'
#' \dontrun{
#'   --> structure of outplut list, representing a plate
#'     |--> metadata for plate (14 fields)
#'     |--> data (list of length 384; 24 cols by 16 rows)
#'        |--> A1
#'           |--> ...
#'       ...
#'        |--> P24
#'           |--> metadata for well (40 fields)
#'           |--> data (list of length 9; 9 images per well)
#'              |--> img_11
#'                 |--> ...
#'             ...
#'              |--> img_33
#'                 |--> Feature1
#'                ...
#'                 |--> FeatureN
#' }
#'
#' @export

importCompletePlate <- function(path, features=NULL) {

  readFeatureFile <- function(filepath) {
    require("R.matlab")
    # nesting of list contains no information
    data <- readMat(filepath)[[1]][[1]][[1]][[1]]
    if (length(data) == 3456) {
      return(data)
    } else {
      cat("expecting .mat files to have 3456 cols, got ", length(data), 
          " instead; skipping ", filepath, "\n", sep="")
      return(NULL)
    }
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