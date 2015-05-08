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
  aggregate <- getCompletePlateAggregate(well)
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
  aggregate <- getCompletePlateAggregate(plate)
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