#' Read all single cell features for wells/plates
#'
#' Read all available single cell feature data for multiple wells and/or
#' multiple plates.
#'
#' @param locations       A list of PlateLocation/WellLocation objects for
#'                        which
#'                        the data is to be returned.
#' @param select.features A vector of strings that are matched with features to
#'                        be kept.
#' @param drop.features   A vector of strings that are matched with features to
#'                        be dropped.
#' @param select.images   A vector of image indices (1:9) indicating which
#'                        images to return for each well.
#'
#' @return If a single PlateLocation/WellLocation is specified as location 
#'         input, the corresponding Data object is returned. For a list of
#'         PlateLocation/WellLocation objects, a list with well/plate data is
#'         returned.
#'
#' @examples
#' plate1 <- PlateLocation("J101-2C")
#' mtor   <- findPlateWellsFromGene("MTOR", "brucella-dp-k[1-3]")
#' plate2 <- PlateLocation("J101-2L")
#' 
#' locations <- c(list(plate1, plate2), mtor)
#' 
#' data <- getSingleCellData(locations)
#' 
#' @export
getSingleCellData <- function(locations, select.features=NULL,
                              drop.features=NULL, select.images=NULL) {
  # check if a single object is requested ot a list of objects
  single <- FALSE
  if(length(locations) == 1) {
    locations <- locations[[1]]
  }
  # validate input
  if(any(class(locations) == "DataLocation")) {
    single <- TRUE
  } else if(!all(sapply(locations,
    function(location) {
      return(any(class(location) == "DataLocation"))
    }
  ))) stop("can only work with a list of plate/well locations")
  if(!(is.null(select.features) | is.vector(select.features,
                                            mode="character"))) {
    stop("exprecting a vector of characters or NULL for select.features")
  }
  if(!(is.null(drop.features) | is.vector(drop.features, mode="character"))) {
    stop("exprecting a vector of characters or NULL for drop.features")
  }
  if(!(is.null(select.images) | is.vector(select.images, mode="numeric"))) {
    stop("exprecting a vector of numeric or NULL for select.images")
  }

  if(single) {
      result <- processSingleCellDataPlate(locations, select.features,
                                           drop.features, select.images)
  } else {
    # get all plates that are of interest
    barcodes <- unique(unlist(lapply(locations, getBarcode)))
    # reorder the list s.t. top level is plates and below that requests per
    # plate
    plates <- lapply(
      barcodes,
      function(barc, locs) {
        res <- lapply(
          locs,
          function(loc, bc) {
            if(getBarcode(loc) == bc) return(loc)
          },
          barc
        )
        return(res[!sapply(res, is.null)])
      },
      locations
    )
    names(plates) <- barcodes
    # process request per plate
    result <- lapply(plates, processSingleCellDataPlate, select.features,
                     drop.features, select.images)
  }
  return(result)
}

#' Read all single cell features for wells/a plate
#'
#' Read all available single cell feature data for multiple wells and/or the
#' whole plate.
#'
#' @param locations A list of PlateLocation/WellLocation objects for which the
#'                  data is to be returned. They all have to lie on the same
#'                  plate.
#' @param select    A vector of strings that are matched with features to be
#'                  kept.
#' @param drop      A vector of strings that are matched with features to be
#'                  dropped.
#' @param images    A vector of image indices (1:9) indicating which images to
#'                  return for each well.
#'              
#' @return A list with well/plate data
#'
#' @examples
#' plate <- PlateLocation("J101-2C")
#' well1 <- WellLocation("J101-2C", "A", 24)
#' well2 <- WellLocation("J101-2C", "B", 12)
#' well3 <- WellLocation("J101-2C", "P", 24)
#' 
#' locations <- list(well1, well2, well3, plate)
#' 
#' data <- processSingleCellDataPlate(locations)
#'
processSingleCellDataPlate <- function(locations, select=NULL, drop=NULL,
                                       images=NULL) {
  # check if a single object is requested ot a list of objects
  single <- FALSE
  if(length(locations) == 1) {
    locations <- locations[[1]]
  }
  # validate input
  if(any(class(locations) == "DataLocation")) {
    single <- TRUE
    barcode <- getBarcode(locations)
  } else if(!all(sapply(locations,
    function(location) {
      return(any(class(location) == "DataLocation"))
    }
  ))) {
    stop("can only work with a list of plate/well locations")
  } else {
    barcodes <- lapply(locations, getBarcode)
    barcode  <- unique(barcodes)
    if(length(barcode) != 1) {
      stop("processSingleCellDataPlate can only handle ",
           "PlateLocation/WellLocation on the same plate.")
    }
  }
  if(!(is.null(select) | is.vector(select, mode="character"))) {
    stop("exprecting a vector of characters or NULL for select")
  }
  if(!(is.null(drop) | is.vector(drop, mode="character"))) {
    stop("exprecting a vector of characters or NULL for drop")
  }
  if(!(is.null(images) | is.vector(images, mode="numeric"))) {
    stop("exprecting a vector of numeric or NULL for images")
  }
  current.plate <- PlateLocation(unlist(barcode))
  # best case: all requested data is cached as well data (cases 1 or 2)
  case <- 0
  if(single) {
    if(file.exists(getScRdsCacheFilename(locations)) &
         any(class(locations) == "WellLocation")) case <- 1
  } else {
    isWell <- sapply(locations, function(location) {
      any(class(location) == "WellLocation")
    })
    if(all(file.exists(sapply(locations, getScRdsCacheFilename))) &
         all(isWell)) case <- 2
  }
  if(case == 1) {
    data <- readRDS(getScRdsCacheFilename(locations))
    message("for plate ", getBarcode(current.plate), " the requested well ",
            "was loaded from cache.")
  } else if(case == 2) {
    data <- lapply(sapply(locations, getScRdsCacheFilename), readRDS)
    message("for plate ", getBarcode(current.plate), " all requested data ",
            "was loaded from cached well files.")
  } else {
    # not all data is in well caches
    if(single) {
      if(any(class(locations) == "WellLocation")) {
        # whole plate has to be imported as well caches have to be written
        data <- PlateData(current.plate)
        case <- 3
        data <- extractWells(data, locations)
        dir  <- dirname(getScRdsCacheFilename(locations))
        if(!dir.exists(dir)) dir.create(dir, recursive=TRUE)
        if(!file.exists(getScRdsCacheFilename(locations))) {
          saveRDS(data, file=getScRdsCacheFilename(locations), compress="xz")
        }
      } else if(any(class(locations) == "PlateLocation")) {
        # plate may only be partially imported (depending on select, drop)
        data <- PlateData(current.plate, select, drop)
      } else stop("can only deal with PlateLocation/WellLocation objects")
    } else {
      if(all(sapply(locations,
        function(location) {
          return(any(class(location) == "PlateLocation"))
        }
      ))) {
        # all list items are plates, therefore no well caches will be written
        data <- PlateData(current.plate, select, drop)
      } else {
        # well caches will be written, import complete plate
        data <- PlateData(current.plate)
        case <- 3
      }
      data <- lapply(
        locations,
        function(location, data) {
          if(any(class(location) == "PlateLocation")) return(data)
          else if (any(class(location) == "WellLocation")) {
            well.data <- extractWells(data, location)
            well.dir  <- dirname(getScRdsCacheFilename(location))
            if(!dir.exists(well.dir)) dir.create(well.dir, recursive=TRUE)
            if(!file.exists(getScRdsCacheFilename(location))) {
              saveRDS(well.data, file=getScRdsCacheFilename(location),
                      compress="xz")
            }
            return(well.data)
          } else stop("can only deal with PlateLocation/WellLocation objects")
        },
        data
      )
    }
  }
  # cases 1 & 2: all data has been read from well caches, features may still
  # have to be dropped
  # case 3: data is still complete, as well caches were written, features may
  # still have to be dropped
  if(case != 0) {
    if(!is.null(select) | !is.null(drop)) {
      if(single) {
        data <- extractFeatures(data, select, drop)
      } else {
        features <- extractFeaturesMatchingHelper(data[1], select, drop)
        data <- lapply(data, extractFeatures, NULL, NULL, features)
      }
    }
  }
  if(!is.null(images)) {
    if(single) data <- extractImages(data, images)
    else data <- lapply(data, extractImages, images)
  }

  return(data)
}