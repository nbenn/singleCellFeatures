#' Constructor for WellLocation objects
#' 
#' Given a well specification consisting of plate barcode, row and column, an
#' S3 object of type WellLocation is created
#'
#' @param plate A plate barcode as a character string
#' @param row   Row specification (a:p, A:P, or 1:16)
#' @param col   Column specification (1:24)
#'
#' @return A WellLocation object: a list with slots plate, row, column, index
#'         (the linearized well index), space, group and experiment
#'
#' @examples
#' well <- WellLocation("J107-2C", "H", 6)
#' 
#' @export
WellLocation <- function(plate, row, col) {
  ind <- getWellIndex1D(row, col)
  row <- getWellIndex2D(ind)$row
  if (!is.character(plate)) stop("plate must be character string")
  data(plateDatabase, envir=environment())
  match <- plate.database[plate.database$Barcode==plate,]
  if(nrow(match) != 1) stop("could not find the specified plate")
  result <- list(plate = plate, row=row, column=col, index=ind, 
                 space=match$Space, group=match$Group, 
                 experiment=match$Experiment)
  return(structure(result, class = c("WellLocation", "DataLocation")))
}

#' Get the well name of a WellLocation object
#'
#' Given a WellLocation object, return the well name
#'
#' @param A WellLocation object
#' 
#' @return The well name of the WellLocation object (as a string)
#' 
#' @examples
#' plate <- WellLocation("J101-2C", "B", 15)
#' well.name <- getWellName(plate)
#' 
#' @export
getWellName <- function(x) {
  UseMethod("getWellName", x)
}
#' @export
getWellName.WellLocation <- function(x) {
  return(paste0(x$row, x$column))
}
#' @export
getWellName.default <- function(x) {
  stop("can only deal with WellLocation objects.")
}

#' Constructor for PlateLocation objects
#' 
#' Given a plate specification (consisting of a plate barcode) an S3 object of
#' type PlateLocation is created
#'
#' @param barcode A plate barcode as a character string
#'
#' @return A PlateLocation object: a list with slots plate, space, group and
#'         experiment
#'
#' @examples
#' plate <- PlateLocation("J107-2C")
#' 
#' @export
PlateLocation <- function(barcode) {
  if (!is.character(barcode)) stop("barcode must be character string")
  data(plateDatabase, envir=environment())
  match <- plate.database[plate.database$Barcode==barcode,]
  if(nrow(match) != 1) stop("could not find the specified plate")
  result <- list(plate = barcode, space=match$Space, group=match$Group, 
                 experiment=match$Experiment)
  return(structure(result, class = c("PlateLocation", "DataLocation")))
}

#' Get the barcode of a PlateLocation/WellLocation object
#'
#' Give a PlateLocation/WellLocation object, return the barcode
#'
#' @param a PlateLocation/WellLocation object
#' 
#' @return The barcode of the PlateLocation/WellLocation object (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' barcode <- getBarcode(plate)
#' 
#' @export
getBarcode <- function(x) {
  UseMethod("getBarcode", x)
}
#' @export
getBarcode.DataLocation <- function(x) {
  return(x$plate)
}
#' @export
getBarcode.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) or ",
       "Data (ImageData/WellData/PlateData) objects.")
}

#' Get the pathogen name of a PlateLocation/WellLocation object
#'
#' Give a PlateLocation/WellLocation object, return the pathogen name
#'
#' @param a PlateLocation/WellLocation object
#' 
#' @return The pathogen name of the PlateLocation/WellLocation object (as a
#'         string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' pathogen <- getPathogen(plate)
#' 
#' @export
getPathogen <- function(x) {
  UseMethod("getPathogen", x)
}
#' @export
getPathogen.DataLocation <- function(x) {
  name1 <- unlist(strsplit(x$group, "_"))[1]
  name2 <- unlist(strsplit(x$experiment, "-"))[1]
  if(name1 != name2) stop("unrecognized naming scheme")
  return(tolower(name1))
}
#' @export
getPathogen.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}

#' Get the local path of a PlateLocation/WellLocation object 
#'
#' Given a PlateLocation/WellLocation object, return the local path
#'
#' @param a PlateLocation/WellLocation object
#' 
#' @return The local path of the PlateLocation/WellLocation object (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- getLocalPath(plate)
#' 
#' @export
getLocalPath <- function(x) {
  UseMethod("getLocalPath", x)
}
#' @export
getLocalPath.DataLocation <- function(x) {
  # load file containing some needed paths
  data(settingsDatabase, envir = environment())
  path <- paste(settings.database$openBIS.data, x$space, x$group, x$experiment, 
                x$plate, sep="/")
  return(path)
}
#' @export
getLocalPath.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}

#' Get the path of a PlateLocation/WellLocation object on openBIS
#'
#' Given a PlateLocation/WellLocation object, return the openBIS location path
#'
#' @param a PlateLocation/WellLocation object
#' 
#' @return The openBIS location path of the PlateLocation/WellLocation object
#'         (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- getOpenBisPath(plate)
#' 
#' @export
getOpenBisPath <- function(x) {
  UseMethod("getOpenBisPath", x)
}
#' @export
getOpenBisPath.DataLocation <- function(x) {
  path <- paste0("/", x$space, "/", x$group, "/", x$experiment, "/", x$plate)
  return(path)
}
#' @export
getOpenBisPath.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}

#' Get path/filename of a PlateLocation/WellLocation object single cell data
#' cache
#'
#' Upon downloading from openBIS and importin to R, all single cell data
#' associated with a PlateLocation/WellLocation is cached in an .rds file. This
#' function gets the filepath of this cache file.
#'
#' @param a PlateLocation/WellLocation object
#' 
#' @return The path/filename of the single cell data cache of
#'         PlateLocation/WellLocation object (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- getScRdsCacheFilename(plate)
#' 
#' @export
getScRdsCacheFilename <- function(x) {
  UseMethod("getScRdsCacheFilename", x)
}
#' @export
getScRdsCacheFilename.WellLocation <- function(x) {
  # load file containing some needed paths
  data(settingsDatabase, envir = environment())
  name <- paste0(x$plate, "_", x$row, x$col, "_sc_all.rds")
  path <- paste(settings.database$openBIS.data, x$space, x$group, x$experiment, 
                x$plate, "WellData", name, sep="/")
  return(path)
}
#' @export
getScRdsCacheFilename.PlateLocation <- function(x) {
  # load file containing some needed paths
  data(settingsDatabase, envir = environment())
  name <- paste0(x$plate, "_sc_all.rds")
  path <- paste(settings.database$openBIS.data, x$space, x$group, x$experiment, 
                x$plate, name, sep="/")
  return(path)
}
#' @export
getScRdsCacheFilename.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}

#' Get the path/filename of the metadata chache for a PlateLocation/WellLocation
#' object
#'
#' Given a PlateLocation/WellLocation object, return the path/filename of the
#' corresponding metadata chache file.
#'
#' @param a PlateLocation/WellLocation object
#' 
#' @return The path/filename of the metadata chache for a 
#'         PlateLocation/WellLocation object (as a string)
#' 
#' @examples
#' plate <- PlateLocation("J101-2C")
#' path <- getMetadataCacheFilename(plate)
#' 
#' @export
getMetadataCacheFilename <- function(x) {
  UseMethod("getMetadataCacheFilename", x)
}
#' @export
getMetadataCacheFilename.DataLocation <- function(x) {
  data(settingsDatabase, envir = environment())
  name <- paste0(x$plate, "_metadata.rds")
  path <- paste(settings.database$openBIS.data, x$space, x$group, x$experiment, 
                x$plate, name, sep="/")
  return(path)
}
#' @export
getMetadataCacheFilename.default <- function(x) {
  stop("can only deal with DataLocation (PlateLocation/WellLocation) objects.")
}

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