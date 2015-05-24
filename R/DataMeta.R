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
  aggregate <- suppressWarnings(PlateAggregate(well))
  if(is.null(aggregate)) {
    empty <- TRUE
  } else {
    if(getBarcode(well) != getBarcode(aggregate$plate)) {
      stop("expecting well and PlateAggregate cache to belong to same plate.")
    }
    aggregate <- aggregate$data[aggregate$data$WellName == getWellName(well),]
    if(nrow(aggregate) == 0) empty <- TRUE
    else empty <- FALSE
  }
  if(empty) {
    warning("no aggregate information available: metadata will be ",
            "incomplete for well ", getWellName(well), " on plate ",
            getBarcode(well), ".")
    return(structure(list(
      plate.barcode       = getBarcode(well),
      plate.quality       = NA,
      experiment.name     = well$experiment,
      experiment.pathogen = getPathogen(well),
      experiment.geneset  = NA,
      experiment.library  = NA,
      well.row            = well$row,
      well.col            = well$column,
      well.index          = well$index,
      well.type           = NA,
      well.quality        = NA,
      gene.name           = NA,
      gene.id             = NA,
      sirna.name          = NA,
      sirna.sequence      = NA,
      sirna.seed          = NA,
      sirna.target        = NA,
      counts.cells        = NA,
      counts.pathogen     = NA,
      counts.infection    = NA
    ), class = c("WellMetadata", "Metadata")))
  } else {
    if("eCount_oPathogen" %in% colnames(aggregate)) {
      n.patho <- aggregate$eCount_oPathogen
    } else if("eCount_oInvasomes" %in% colnames(aggregate)) {
      n.patho <- aggregate$eCount_oInvasomes
    } else {
      n.patho <- NA
    }
    if("dInfectionDT_eCount" %in% colnames(aggregate)) {
      n.inf <- aggregate$dInfectionDT_eCount
    } else {
      n.inf <- NA
    }
    return(structure(list(
      plate.barcode       = aggregate$Barcode,
      plate.quality       = aggregate$PLATE_QUALITY_STATUS,
      experiment.name     = aggregate$Experiment,
      experiment.pathogen = aggregate$PATHOGEN,
      experiment.geneset  = aggregate$GENESET,
      experiment.library  = aggregate$LIBRARY,
      well.row            = aggregate$WellRow,
      well.col            = as.integer(aggregate$WellColumn),
      well.index          = getWellIndex1D(aggregate$WellRow,
                                           aggregate$WellColumn, NULL),
      well.type           = aggregate$WellType,
      well.quality        = aggregate$WELL_QUALITY_STATUS,
      gene.name           = aggregate$Name,
      gene.id             = aggregate$ID_manufacturer,
      sirna.name          = aggregate$ID_openBIS,
      sirna.sequence      = aggregate$Sequence_antisense_5_3,
      sirna.seed          = aggregate$Seed_sequence_antisense_5_3,
      sirna.target        = aggregate$Sequence_target_sense_5_3,
      counts.cells        = aggregate$eCount_oCells,
      counts.pathogen     = n.patho,
      counts.infection    = n.inf
    ), class = c("WellMetadata", "Metadata")))
  }
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
  checkIfUnique <- function(x) {
    res <- unique(x)
    if(length(res) != 1) warning("non unique within-plate variable")
    return(res[1])
  }

  # input validation
  if(!any(class(plate) == "PlateLocation")) {
    stop("can only work with PlateLocation objects")
  }
  aggr <- suppressWarnings(PlateAggregate(plate))
  if(is.null(aggr)) {
    warning("no aggregate information available: metadata will be ",
            "incomplete for plate ", getBarcode(plate), ".")
    return(structure(list(
      plate.barcode        = plate$plate,
      plate.quality        = NA,
      plate.type           = NA,
      experiment.space     = plate$space,
      experiment.group     = plate$group,
      experiment.name      = plate$experiment,
      experiment.pathogen  = getPathogen(plate),
      experiment.geneset   = NA,
      experiment.replicate = NA,
      experiment.library   = NA,
      experiment.batch     = NA
    ), class = c("PlateMetadata", "Metadata")))
  } else {
    if(getBarcode(plate) != getBarcode(aggr$plate)) {
      stop("expecting plate and Plateaggr cache to have same barcode")
    }
    aggr <- aggr$data
    return(structure(list(
      plate.barcode        = checkIfUnique(aggr$Barcode),
      plate.quality        = checkIfUnique(aggr$PLATE_QUALITY_STATUS),
      plate.type           = checkIfUnique(aggr$PLATE_TYPE),
      experiment.space     = checkIfUnique(aggr$Space),
      experiment.group     = checkIfUnique(aggr$Group),
      experiment.name      = checkIfUnique(aggr$Experiment),
      experiment.pathogen  = checkIfUnique(aggr$PATHOGEN),
      experiment.geneset   = checkIfUnique(aggr$GENESET),
      experiment.replicate = checkIfUnique(aggr$REPLICATE),
      experiment.library   = checkIfUnique(aggr$LIBRARY),
      experiment.batch     = checkIfUnique(aggr$BATCH)
    ), class = c("PlateMetadata", "Metadata")))
  }
}

#' Constructor for PlateAggregate objects
#' 
#' Given a PlateLocation/WellLocation object, return the complete and
#' unprocessed plate aggregate information.
#'
#' @param dl A PlateLocation/WellLocation (DataLocation) object
#'
#' @return A list with slots for PlateLocation (plate) and aggregate data
#'        (data).
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
    cache  <- readRDS(getCacheFilenameMeta(dl))
    result <- structure(list(plate=dl, data=cache),
                        class = c("PlateAggregate", "Metadata"))
    return(result)
  }
  # figure out pathogen name
  patho.name <- getPathogen(dl)
  # get path info
  config <- configGet()
  # get genome and kinome aggreagtes
  aggregate.name <- list.files(
    path=paste0(config$dataStorage$metaDir, "/", "Aggregates"),
    pattern=paste0(patho.name, "report_.*\\.csv"),
    ignore.case=TRUE, full.names=TRUE
  )
  if (length(aggregate.name) != 1) {
    stop("found ", length(aggregate.name), " aggregate files instead of 1.")
  }
  # load metadata files
  aggregate.file <- read.delim(aggregate.name, stringsAsFactors=FALSE)

  # drop all except the dl of interest
  aggregate.file <- aggregate.file[aggregate.file$Barcode == getBarcode(dl),]
  if (nrow(aggregate.file) != 384) {
    warning("found genome aggregate information on ", nrow(aggregate.file),
            " wells instead of 384.")
  }
  if(nrow(aggregate.file) > 0) {
    result <- structure(list(plate=dl, data=aggregate.file),
                        class = c("PlateAggregate", "Metadata"))
    saveToCache(result)
    return(result)
  } else {
    warning("could not find any aggregate information for plate ",
            getBarcode(dl))
    return(NULL)
  }
}
