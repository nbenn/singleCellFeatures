#' Find wells corresponding to a gene
#' 
#' Find all wells on all plates of a given pathogen/experiment that contain
#' the specified gene
#'
#' @param gene       String/integer specifying the gene (case sensitive)
#' @param experiment Name of the experiment to be considered (should start with
#'                   pathogen name rest is optional and used for narrowing down
#'                   the number of results, all converted to uppercase)
#' @param plates     An optional vector of plate names to further restrict the
#'                   results
#' @param verbose    Verbosity argument (FALSE/TRUE)
#'
#' @return A table containing columns for experiment, barcode, well row and well
#'         column
#'
#' @examples
#' set1 <- findPlateWellsFromGene("MTOR", "brucella-du-k")
#' set2 <- findPlateWellsFromGene(2475, "brucella-au-k[1-3]")
#' set3 <- findPlateWellsFromGene("SCRAMBLED", "brucella-du-k")
#'
#' @export

findPlateWellsFromGene <- function(gene, experiment, plates=NULL,
                                   verbose=FALSE) {
  # ensure correct case for pathogen/experiment
  patho.upper <- toupper(unlist(strsplit(experiment, "-"))[1])
  patho.lower <- tolower(patho.upper)
  patho.camel <- paste0(toupper(substring(patho.lower, 1, 1)),
                        substring(patho.lower, 2))
  experiment  <- toupper(experiment)

  dataset.name <- paste0("wellDatabase", patho.camel)
  object.name  <- paste0("well.database.", patho.lower)

  data(plateDatabase, envir=environment())
  data(list=dataset.name, envir=environment())
  well.database <- get(object.name)

  curr.plates <- plate.database[grep(experiment, plate.database$Experiment),]
  if(!is.null(plates)) {
    intersection <- intersect(plates, curr.plates$Barcode)
    if(length(intersection) < 1) {
      stop("no plates found with the given restrictions")
    }
    curr.plates <- curr.plates[which(curr.plates$Barcode %in% intersection),]
  }
  if(is.character(gene)) {
    curr.wells <- well.database[well.database$Name == gene, ]
  } else if(is.integer(gene) | is.numeric(gene)) {
    gene <- as.integer(gene)
    curr.wells <- well.database[well.database$ID == gene, ]
  } else stop("cannot make sense of gene argument (expecting int or char)")

  res.tab <- curr.wells[which(curr.wells$Barcode %in% curr.plates$Barcode),]
  if(nrow(res.tab) == 0) stop("no matching wells found.")
  res.lst <- apply(res.tab, 1, function(row) {
    WellLocation(row[["Barcode"]], row[["WellRow"]], row[["WellColumn"]])
  })
  return(res.lst)
}