#' A wrapper around saveRDS to enable multicore compression
#' 
#' Piping the datastream to be saved through pigz enables the use of multiple
#' cores for compression.
#'
#' @param object  An object to be saved
#' @param file    The filename of the resulting .rds file
#' @param threads the number of threads to be used for compression (default:
#'                all-1)
#'
#' @return NULL (invisibly). The object is saved to the filesystem.
#'
#' @examples
#' obj <- list(a=1, b=2)
#' saveRDSMC(obj, "~/obj.rds")
#' 
#' @export
saveRDSMC <- function(object, file, threads=getNumCores()) {
  message("using ", threads, " threads for compression.")
  #con <- pipe(paste0("xz -T", threads, " -9 -f > ", file), "wb")
  con <- pipe(paste0("pigz -p ", threads, " -9 -f > ", file), "wb")
  saveRDS(object, file = con)
  on.exit(if(exists("con")) close(con))
}

#' A wrapper around readRDS for faster reading of compressed files
#' 
#' Although decompressing files with pigz is not parallelized, it is faster to
#' run the decompression in a separate process.
#'
#' @param file The file to be read
#'
#' @return The object saved as file
#'
#' @examples
#' data <- readRDSMC("~/obj.rds")

#' @export
readRDSMC <- function(file) {
  #con <- pipe(paste0("xz -d -k -c ", file))
  con <- pipe(paste0("pigz -d -k -c ", file))
  object <- tryCatch({
    readRDS(file = con)
  },
  error = function(err) {
    stop("could not read file\n", file, ":\n", err)
  },
  finally = {
    if(exists("con")) close(con)
  })
  return(object)
}

#' Get the number of available cores
#' 
#' In case, the package is used in an LSF environment, detectCores() will
#' report the number of physical cores insted of the numer allocated by LSF.
#'
#' @return The number of available cores.
#'
#' @examples
#' n.cores <- getNumCores()

#' @export
getNumCores <- function() {
  n.cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC"))
  if(is.na(n.cores)) {
    n.cores <- detectCores()
  }
  if(n.cores > detectCores() | n.cores < 1) {
    stop("illegal number of cores (", n.cores, ") specified.")
  }
  return(n.cores)
}
