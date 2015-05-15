#' Constructor for WellMetadata objects
#' 
#' Given a WellLocation object, the corresponding metadata is compiled
#'
#' @param well A WellLocation object
#'
#' @return A WellMetadata object: a list with slots describing conditions in
#'         a well.
#'
#' @examples
#' meta <- WellMetadata(WellLocation("J107-2C", "H", 6))
#' 
#' @export
WellMetadata <- function(well) {
  # input validation
  if(!any(class(well) == "WellLocation")) {
    stop("can only work with WellLocation objects")
  }
  aggregate <- PlateAggregate(well)
  gen.aggr  <- aggregate$gen[aggregate$gen$WellName == getWellName(well),]
  if(is.null(aggregate$kin)) {
    warning("no kinome aggregate information available: metadata will be ",
            "incomplete for well ", getWellName(well), " on plate ",
            getBarcode(well), ".")
    result <- list(
      barcode             = gen.aggr$Barcode,
      row                 = gen.aggr$WellRow,
      col                 = as.integer(gen.aggr$WellColumn),
      index               = getWellIndex1D(gen.aggr$WellRow,
                                           gen.aggr$WellColumn, NULL),
      well.type           = gen.aggr$WellType,
      control.type        = NULL,
      control.description = NULL,
      gene.id             = ifelse(gen.aggr$ID == "unknown", NA,
                                   as.integer(gen.aggr$ID)),
      gene.name           = gen.aggr$Name,
      sirna.manufacturer  = gen.aggr$LIBRARY,
      sirna.code          = gen.aggr$Catalog_number,
      sirna.sequence      = gen.aggr$Sequence_antisense_5_3
    )
  } else {
    kin.aggr  <- aggregate$kin[aggregate$kin$WellName == getWellName(well),]
    result <- list(
      barcode             = gen.aggr$Barcode,
      row                 = gen.aggr$WellRow,
      col                 = as.integer(gen.aggr$WellColumn),
      index               = getWellIndex1D(gen.aggr$WellRow,
                                           gen.aggr$WellColumn, NULL),
      well.type           = gen.aggr$WellType,
      control.type        = ifelse(kin.aggr$WellType == "siRNA", "none",
                                   kin.aggr$WellType),
      control.description = kin.aggr$ControlType,
      gene.id             = ifelse(gen.aggr$ID == "unknown", NA,
                                   as.integer(gen.aggr$ID)),
      gene.name           = gen.aggr$Name,
      sirna.manufacturer  = gen.aggr$LIBRARY,
      sirna.code          = gen.aggr$Catalog_number,
      sirna.sequence      = gen.aggr$Sequence_antisense_5_3
    )
  }
  return(structure(result, class = c("WellMetadata", "Metadata")))
}

#' Constructor for PlateMetadata objects
#' 
#' Given a PlateLocation object, the corresponding metadata is compiled
#'
#' @param well A PlateLocation object
#'
#' @return A PlateMetadata object: a list with slots describing conditions on
#'         a plate.
#'
#' @examples
#' meta <- PlateMetadata(PlateLocation("J107-2C"))
#' 
#' @export
PlateMetadata <- function(plate) {
  checkIfUnique <- function(x, name) {
    res <- unique(x)
    if(length(res) != 1) warning("non unique within-plate variable: ", name)
    return(res[1])
  }

  # input validation
  if(!any(class(plate) == "PlateLocation")) {
    stop("can only work with PlateLocation objects")
  }
  aggregate <- PlateAggregate(plate)
  gen.aggr  <- aggregate$gen
  result <- list(
    barcode         = checkIfUnique(gen.aggr$Barcode,        "Barcode"),
    space           = checkIfUnique(gen.aggr$Space,          "Space"),
    group           = checkIfUnique(gen.aggr$Group,          "Group"),
    experiment      = checkIfUnique(gen.aggr$Experiment,     "Experiment"),
    experiment.type = checkIfUnique(gen.aggr$ExperimentType, "ExperimentType"),
    pathogen        = checkIfUnique(gen.aggr$PATHOGEN,       "Pathogen"),
    gene.set        = checkIfUnique(gen.aggr$GENESET,        "GeneSet"),
    replicate       = checkIfUnique(gen.aggr$REPLICATE,      "Replicate"),
    manufacturer    = checkIfUnique(gen.aggr$LIBRARY,        "Library"),
    plate.type      = checkIfUnique(gen.aggr$PLATE_TYPE,     "PlateType"),
    batch           = checkIfUnique(gen.aggr$BATCH,          "Batch")
  )
  return(structure(result, class = c("PlateMetadata", "Metadata")))
}

#' Constructor for PlateAggregate objects
#' 
#' Given a PlateLocation/WellLocation object, return the complete and
#' unprocessed genome and kinome plate aggregate information.
#'
#' @param dl A PlateLocation/WellLocation (DataLocation) object
#'
#' @return A list with slots for kinome (kin) and genome (gen) aggregate data.
#'
#' @examples
#' plate <- PlateLocation("J101-2C")
#' aggregate <- PlateAggregate(plate)
#' 
#' @export
PlateAggregate <- function(dl) {
  if(!any(class(dl) == "DataLocation")) {
    stop("can only work woth DataLocation objects")
  }
  if(any(class(dl) == "WellLocation")) {
    dl <- convertToPlateLocation(dl)
  }
  # search for aggregate chache file
  if(file.exists(getCacheFilenameMeta(dl))) {
    cache <- readRDS(getCacheFilenameMeta(dl))
    result <- c(cache, list(met=dl))
    return(structure(result, class = c("PlateAggregate", "Metadata")))
  } else {
    # figure out pathogen name
    patho.name <- getPathogen(dl)
    # get path info
    config <- configGet()
    # get genome and kinome aggreagtes
    gen.file <- list.files(path=config$dataStorage$genome, pattern=patho.name,
                           ignore.case=TRUE, full.names=TRUE)
    kin.file <- list.files(path=config$dataStorage$kinome, pattern="\\.csv$",
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

    # drop all except the dl of interest
    gen.aggr <- gen.aggr[gen.aggr$Barcode == getBarcode(dl),]
    kin.aggr <- kin.aggr[kin.aggr$Barcode == getBarcode(dl),]
    if (nrow(gen.aggr) != 384) {
      warning("found genome aggregate information on ", nrow(gen.aggr),
              " wells instead of 384.")
    }
    if (nrow(kin.aggr) != 384) {
      warning("found kinome aggregate information on ", nrow(kin.aggr),
              " wells instead of 384.")
    }
    if (nrow(kin.aggr) == 0) kin.aggr <- NULL
    result <- structure(list(gen=gen.aggr, kin=kin.aggr, met=dl),
                        class = c("PlateAggregate", "Metadata"))
    saveToCache(result)
    return(result)
  }
}
