#' Get all local plates
#'
#' Return plate location objects for each locally present plate.
#' 
#' @return A list of PlateLocation objects.
#' 
#' @examples
#' plates <- getAllLocalPlates()
#' 
#' @export
getAllLocalPlates <- function() {
  config <- configGet()
  all.dirs <- list.dirs(config$dataStorage$dataDir, full.names=TRUE,
                        recursive=TRUE)
  all.dirs <- lapply(all.dirs, function(dir) unlist(strsplit(dir, "/")))
  baselength <- length(all.dirs[[1]])
  barcodes <- lapply(all.dirs, function(dir) dir[baselength + 4])
  barcodes <- unique(unlist(barcodes))
  barcodes <- barcodes[!is.na(barcodes)]
  return(lapply(barcodes, PlateLocation))
}

#' Rebuild all local plates
#'
#' For all local plates, retry downloading missing features and rebuild the
#' corresponding caches.
#' 
#' @return A list booleans indicating whether the plate was successfully
#'         rebuilt. Rebuilding complete plates counts as failure.
#' 
#' @examples
#' res <- rebuildAllMatDataCaches()
#' 
#' @export
rebuildAllMatDataCaches <- function() {
  plates <- getAllLocalPlates()
  lapply(plates, function(plate) {
    tryCatch({
      message("processing plate ", getBarcode(plate))
      temp <- MatData(plate)
      temp <- rebuildCache(temp)
      rm(temp)
      TRUE
    },
    error = function(err) {
      message(err$message)
      return(FALSE)
    })
  })
}