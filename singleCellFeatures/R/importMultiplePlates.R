#' Read single cell feature data for multiple plates
#' 
#' Read all available/the specified set of single cell feature data for all 
#' plates/the specified set below the dir path
#'
#' @param path     Folder path to folder corresponding to an experiment
#' @param plates   Names of plates belonging to the given experiment that are
#'                 to be imported
#' @param features Vector holding strings corresponding to features to be 
#'                 imported
#'
#' @return A nested list with the top level hierarchy corresponding to a plate 
#'         and the level below to a feature
#'
#' @examples
#' path     <- paste("/Users/nbennett/Polybox/MasterThesis/openBISDownload/",
#'                   "INFECTX_PUBLISHED/ADENO_TEAM", sep="")
#' plates   <- c("KB2-03-1I", "KB2-03-1J")
#' features <- c("Cells.AreaShape_Area", "Cells.AreaShape_Eccentricity")
#' dat.exp1 <- importMultiplePlates(path, NULL, features)
#' dat.exp2 <- importMultiplePlates(path, plates, features)
#'
#' @export

importMultiplePlates <- function(path, plates=NULL, features=NULL) {  
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