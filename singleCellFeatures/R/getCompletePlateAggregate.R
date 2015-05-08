#' Get the complete plate aggregate for a PlateLocation/WellLocation object
#'
#' Given a PlateLocation/WellLocation object, return the complete and
#' unprocessed genome and kinome plate aggregate information.
#'
#' @param a PlateLocation/WellLocation object
#' 
#' @return A list with slots for kinome (kin) and genome (gen) aggregate data
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' aggregate <- getCompletePlateAggregate(plate)
#' 
#' @export
getCompletePlateAggregate <- function(x) {
  UseMethod("getCompletePlateAggregate", x)
}

#' @export
getCompletePlateAggregate.DataLocation <- function(x) {
  # search for aggregate chache file
  if(file.exists(getMetadataCacheFilename(x))) {
    return(readRDS(getMetadataCacheFilename(x)))
  } else {
    # figure out pathogen name
    patho.name <- getPathogen(x)
    # get path info
    data(settingsDatabase, envir = environment())
    # get genome and kinome aggreagtes
    gen.file <- list.files(path=settings.database$gen.aggr, pattern=patho.name,
                           ignore.case=TRUE, full.names=TRUE)
    kin.file <- list.files(path=settings.database$kin.aggr, pattern="\\.csv$",
                           full.names=TRUE)
    if (length(gen.file) != 1) {
      stop("found ", length(gen.file), " genome aggregate files instead of 1.")
    }
    if (length(kin.file) != 1) {
      stop("found ", length(kin.file), " kinome aggregate files instead of 1.")
    }    
    # load metadata files
    gen.aggr  <- read.delim(gen.file, as.is=TRUE)
    kin.aggr  <- read.table(kin.file, header = TRUE, sep = ";", fill = TRUE,
                            stringsAsFactors = FALSE, comment.char = "")
    
    # drop all except the plate of interest
    gen.aggr <- gen.aggr[gen.aggr$Barcode == getBarcode(x),]
    kin.aggr <- kin.aggr[kin.aggr$Barcode == getBarcode(x),]
    if (nrow(gen.aggr) != 384) {
      warning("found genome aggregate information on ", nrow(gen.aggr), 
              " wells instead of 384.")
    }
    if (nrow(kin.aggr) != 384) {
      warning("found kinome aggregate information on ", nrow(kin.aggr), 
              " wells instead of 384.")
    }
    if (nrow(kin.aggr) == 0) kin.aggr <- NULL
    result <- list(gen=gen.aggr, kin=kin.aggr)
    saveRDS(result, getMetadataCacheFilename(x))
  }
  return(result)
}

#' @export
getCompletePlateAggregate.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}