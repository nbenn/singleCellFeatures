#' Find wells corresponding to a gene
#' 
#' Find all wells on all plates of a given pathogen/experiment that contain
#' the specified gene
#'
#' @param experiment Name of the experiment to be considered (should start with
#'                   pathogen name, all converted to uppercase)
#'                   to be imported
#' @param gene       String specifying the gene (case sensitive)
#'
#' @return A table containing columns for experiment, barcode, well row and well
#'         column
#'
#' @examples
#' mtor.bruc.du <- findPlateWellsFromGene("brucella-du-k", "MTOR")
#' scra.bruc.du <- findPlateWellsFromGene("brucella-du-k", "SCRAMBLED")
#'
#' @export

findPlateWellsFromGene <- function(experiment, gene) {
  # load file containing paths/passwords
  data(paths, envir = environment())
  # ensure correct case
  patho.upper <- toupper(unlist(strsplit(experiment, "-"))[1])
  experiment  <- toupper(experiment)
  
  # sef filename for the needed genome wide aggregate file (there is one for
  # each pathogen)
  if(patho.upper=="ADENO") {
    gen.name <- paste(filepaths$gen.aggr, "AdenoReport_20141008.csv", sep="/")
  } else if(patho.upper=="BARTONELLA") {
    gen.name <- paste(filepaths$gen.aggr, "BartonellaReport_20141010.csv", 
                      sep="/")
  } else if(patho.upper=="BRUCELLA") {
    gen.name <- paste(filepaths$gen.aggr, "BrucellaReport_20141010.csv",
                      sep="/")
  } else stop("unknown pathogen")

  # set filename for the kinome wide aggregate file (one for all pathogens; not
  # the case right now...)
  kin.name <- paste(filepaths$kin.aggr, 
                    "InfectX\ Kinome\ Data\ -\ Well\ aggregate.csv", sep="/")
  # load the two aggregate files
  gen.data <- read.table(gen.name, header = TRUE, sep = "\t", fill = TRUE, 
                         stringsAsFactors = FALSE, comment.char = "")
  kin.data <- read.table(kin.name, header = TRUE, sep = ";", fill = TRUE, 
                         stringsAsFactors = FALSE, comment.char = "")
  # search for gene and reduce number of columns
  gene.gen <- gen.data[gen.data$Name==gene,
                       c('Experiment', 'Barcode', 'WellRow', 'WellColumn')]
  gene.kin <- kin.data[kin.data$GeneName==gene,
                       c('Experiment', 'Barcode', 'WellRow', 'WellColumn')]
  gene.kin$Experiment <- unlist(lapply(strsplit(gene.kin$Experiment, "/"), 
                                       tail, n=1))
  # search for experiment
  expe.gen <- gene.gen[grep(experiment, gene.gen$Experiment),]
  expe.gen <- expe.gen[order(expe.gen$Barcode),]
  row.names(expe.gen) <- NULL
  expe.kin <- gene.kin[grep(experiment, gene.kin$Experiment),]
  expe.kin <- expe.kin[order(expe.kin$Barcode),]
  row.names(expe.kin) <- NULL
  
  # best case: kinome & genome results are identical as they should be if a
  # kinase is searched for
  if(identical(expe.gen, expe.kin)) return(expe.gen)
  # no results found in genome aggreagte: happens if search is for a control
  # (eg. scrambled), as the type of control is not specified in genome wide
  # aggregates
  else if(nrow(expe.gen)==0) {
    return(expe.kin)
  } else {
    # when the results from a genome wide aggregate differs from the kinome
    # wide, as a workaround the larger of the two is returned
    cat("\nAttention: kinome and genome results differ.\n\nGenome:\n")
    print(expe.gen)
    cat("\nKinome:\n")
    print(expe.kin)
    if(nrow(expe.gen) > nrow(expe.kin)) {
      cat("\nReturning the larger of the two (genome).\n\n")
      return(expe.gen)
    } else {
      cat("\nReturning the larger of the two (kinome).\n\n")
      return(expe.kin)
    }
  }
}