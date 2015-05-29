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
      suppressWarnings(temp <- MatData(plate))
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

#' Remove all local well caches
#'
#' Due to generating well caches being fairly inexpensive and the circumstance
#' that they are downstream to Data object creation (changes in Data object
#' specifications often require all well caches to be rebuilt), this function
#' purges all well caches. All .rds files in well cache directories are deleted
#' and if the direcory is empty afterwards, it is also removed.
#' 
#' @return NULL (invisibly).
#' 
#' @examples
#' deleteAllWellCaches()
#' 
#' @export
deleteAllWellCaches <- function() {
  l_ply(getAllLocalPlates(), function(plate) {
    cacheDir   <- dirname(getCacheFilenameData(WellLocation(plate, "A", 1)))
    cacheFiles <- list.files(path=cacheDir, pattern="\\.rds", full.names=TRUE)
    if(length(cacheFiles) > 0) {
      message("deleting well caches for plate ", getBarcode(plate))
      unlink(cacheFiles)
      if(length(list.files(path=cacheDir, pattern=".*")) == 0) {
        unlink(cacheDir, recursive=TRUE)
      }
    }
  })
  invisible(NULL)
}

#' Remove all metadata caches
#'
#' Due to generating metadata caches being fairly inexpensive and the
#' circumstance that they are downstream to Data object creation (changes in
#' Data/Metadata object specifications often require all metadata caches to be
#' rebuilt), this function purges all metadata caches.
#' 
#' @return NULL (invisibly).
#' 
#' @examples
#' deleteAllWellCaches()
#' 
#' @export
deleteAllMetadataCaches <- function() {
  l_ply(getAllLocalPlates(), function(plate) {
    message("deleting metadata caches for plate ", getBarcode(plate))
    unlink(getCacheFilenameMeta(plate))
  })
  invisible(NULL)
}
