.onLoad <- function(libname, pkgname) {
  op <- options()
  op.singleCellFeatures <- list(
    singleCellFeatures.configPath = paste0("~", "/",
                                           ".singleCellFeaturesConfig")
  )
  toset <- !(names(op.singleCellFeatures) %in% names(op))
  if(any(toset)) options(op.singleCellFeatures[toset])
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  wellDatabaseCoverage(verbose=FALSE, is.startup=TRUE)
  config.file <- getOption("singleCellFeatures.configPath")
  if(!file.exists(config.file)) {
    packageStartupMessage("\nthe file ", config.file, "\ndoes not exist.\n",
                          "1. run configPathSet(\"path/to/where/you/want/the/",
                          "configFile.yaml\")\n2. run configInit()\n",
                          "3. edit the file.")
  }
  invisible(NULL)
}

.onUnload <- function(libpath) {
  options(error=NULL)
  invisible(NULL)
}