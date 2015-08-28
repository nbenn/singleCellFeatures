#' Set the config file path
#' 
#' Sets the path to the yaml configuration file for the singleCellFeatures
#' package.
#'
#' @param path The new path of the config file.
#'
#' @return NULL (invisibly)
#'
#' @examples
#' configPathSet("path/to/config.yaml")
#' 
#' @export
configPathSet <- function(path) {
  if(!file.exists(path)) {
    warning("the file ", path, "\ndoes not exist.\nEither retry setting the ",
            "path using configPathSet(\"path/to/config.yaml\"),\nor run ",
            "configInit()  and edit the resulting file.")
  }
  options(singleCellFeatures.configPath = path)
  invisible(NULL)
}

#' Get the config file path
#' 
#' Gets the path to the yaml configuration file for the singleCellFeatures
#' package.
#'
#' @return The config file path
#'
#' @examples
#' path <- configPathGet()
#' 
#' @export
configPathGet <- function() {
  path <- getOption("singleCellFeatures.configPath")
  if(!file.exists(path)) {
    warning("the file ", path, "\ndoes not exist.\nEither retry setting the ",
            "path using configPathSet(\"path/to/config.yaml\"),\nor run ",
            "configInit()  and edit the resulting file.")
  }
  return(path)
}

#' Creates a config file template 
#' 
#' Saves a config file template to the current config file path.
#'
#' @return NULL (invisibly)
#'
#' @examples
#' configInit()
#' 
#' @export
configInit <- function() {
  suppressWarnings(path <- configPathGet())
  if(!dir.exists(dirname(path))) {
    stop("the directory ", dirname(path), " does not exist.\n",
         "please create it first.")
  }
  config <- list(
    dataStorage = list(
      dataDir = "path/to/data/dir",
      metaDir = "path/to/metadata/dir"
    ),
    beeDownloader = list(
      executable  = "path/to/trunk/openBIS/Tools/BeeDataSetDownloader",
      beeSoftsrc = "path/to/trunk"
    ),
    openBIS = list(
      username = "user",
      password = "password"
    ),
    singleCellFeatures = list(
      sourceDir = "path/to/source"
    )
  )
  message("writing config file to ", path)
  write(as.yaml(config), path)
  invisible(NULL)
}

#' Loads the config object
#' 
#' Loads the yaml file saved at the current config file location.
#'
#' @return The config object.
#'
#' @examples
#' config <- configGet()
#' 
#' @export
configGet <- function() {
  suppressWarnings(path <- configPathGet())
  if(!file.exists(path)) {
    stop("the file ", path, "\ndoes not exist.")
  }
  return(yaml.load_file(path))
}

#' Overwrites the config file
#' 
#' Saves the specified config object to the current config file location. Yaml
#' was chosen as config file format due to its human read/writability, so direct
#' editing of the file is possible. This function is available for convenience
#' nevertheless.
#'
#' @param x The new config object.
#'
#' @return NULL (invisibly)
#'
#' @examples
#' conf <- configGet()
#' # edit the list in R
#' configSet(conf)
#' 
#' @export
configSet <- function(config) {
  if(!is.list(config)) stop("expecting a list as config object.")
  suppressWarnings(path <- configPathGet())
  if(!dir.exists(dirname(path))) {
    stop("the directory ", dirname(path), " does not exist.\n",
         "please create it first.")
  }
  message("writing config file to ", path)
  write(as.yaml(config), path)
  invisible(NULL)
}

#' Reload the singleCellFeatures package
#' 
#' Detaches and reloads the singleCellFeatures package.
#'
#' @return NULL (invisibly).
#'
#' @examples
#' reloadSingleCellFeatures()
#' 
#' @export
reloadSingleCellFeatures <- function() {
  detach("package:singleCellFeatures", unload = TRUE, character.only = TRUE)
  library("singleCellFeatures", character.only = TRUE)
  invisible(NULL)
}

#' Get amount of storage used by an object
#' 
#' Prints the amount of storage used by the specified object.
#'
#' @param x The object of interest
#'
#' @return The amount of storage taken up by the object as a string
#'
#' @examples
#' plate <- WellLocation("J101-2C", "B", 15)
#' size(plate)
#' 
#' @export
size <- function(x) {
  format(object.size(x), units="auto")
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
