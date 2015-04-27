#' Read all single cell features for wells/plates
#'
#' Read all available single cell feature data for multiple wells and/or
#' multiple plates.
#'
#' @param locations       A list of PlateLocation/WellLocation objects for which
#'                        the data is to be returned.
#' @param select.features A vector of strings that are matched with features to
#'                        be kept.
#' @param drop.features   A vector of strings that are matched with features to
#'                        be dropped.
#' @param select.images   A vector of image indices (1:9) indicating which
#'                        images to return for each well.
#'
#' @return A list with well/plate data
#'
#' @examples
#' plate1 <- PlateLocation("J101-2C")
#' mtor   <- findPlateWellsFromGene("MTOR", "brucella-du-k")
#' plate2 <- PlateLocation("J101-2L")
#' 
#' locations <- c(list(plate1, plate2), mtor)
#' 
#' data <- getSingleCellData(locations)
#' 
#' @export
getSingleCellData <- function(locations, select.features=NULL,
                              drop.features=NULL, select.images=NULL) {
  # vecidate input: chack if it's a list; the type of objects is checked in 
  # getBarcode()
  if(typeof(locations) != "list") {
    stop("expecting a list of plate/well locations")
  }
  if(!(is.vector(features, mode="character") | is.null(features))) {
    stop("exprecting a vector of characters or NULL for features")
  }
  # get all plates that are of interest
  barcodes <- unique(unlist(lapply(locations, getBarcode)))
  # reorder the list s.t. top level is plates and below that requests per plate
  plates <- lapply(barcodes, function(barc, locs) {
    res <- lapply(locs, function(loc, bc) {
      if(getBarcode(loc) == bc) return(loc)
    }, barc)
    return(res[!sapply(res, is.null)])
  }, locations)
  names(plates) <- barcodes
  # process request per plate
  result <- lapply(plates, processSingleCellDataPlate, select.features,
                   drop.features)
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
  # input validation
  barcodes <- lapply(locations, getBarcode)
  barcode  <- unique(barcodes)
  if(length(barcode) != 1) {
    stop("processSingleCellDataPlate can only handle ",
         "PlateLocation/WellLocation on the same plate.") 
  }
  if(!(is.vector(select, mode="character") | is.null(select))) {
    stop("exprecting a vector of characters or NULL for select")
  }
  if(!(is.vector(drop, mode="character") | is.null(drop))) {
    stop("exprecting a vector of characters or NULL for drop")
  }
  current.plate <- PlateLocation(unlist(barcode))
  # best case: all requested data is cached as well data
  if(all(file.exists(sapply(locations, getScRdsCacheFilename))) &
     all(sapply(locations, function(location) {
       any(class(location) == "WellLocation")}))) {
    data <- lapply(sapply(locations, getScRdsCacheFilename), readRDS)
    message("for plate ", getBarcode(current.plate), " all requested data ",
            "was loaded from cached well files.")
  } else {
    # check if plate is cached
    if(!file.exists(getScRdsCacheFilename(current.plate))) {
      # plate has to be downloaded/imported
      path <- getLocalPath(current.plate)
      # search for matlab files
      matlab.files <- list.files(path=path, pattern="\\.mat$", 
                                 full.names=TRUE, recursive=TRUE)
      if (length(matlab.files) == 0) {
        # plate has to be downloaded
        message("downloading plate ", getBarcode(current.plate),
                " from openBIS.")
        data(settingsDatabase, envir = environment())
        # get password from file
        pwd.openBIS <- readRDS(settings.database$openBIS.pwd)
        bee.command <- paste0(
          "cd ", settings.database$bee.dir, "; export BEESOFTSRC='", 
          settings.database$bee.softsrc, "'; ./BeeDataSetDownloader.sh ",
          "--user '", settings.database$openBIS.user, "' --password '",
          pwd.openBIS, "' --outputdir '", settings.database$openBIS.data,
          "' --plateid '^", getOpenBisPath(current.plate), "' --files '.*.mat'",
          " --verbose '10'")
        system(bee.command, ignore.stdout = TRUE)
      }
    }
    data <- importSingleCellDataPlate(current.plate)
    data <- lapply(locations, function(location, data) {
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
    }, data)
  }
  if(!is.null(select) | !is.null(drop)) {
    features <- extractFeaturesMatchingHelper(data[1], select, drop)
    data <- lapply(data, extractFeatures, NULL, NULL, features)
  }
  if(!is.null(images)) {
    data <- lapply(data, extractImages, images)
  }
  return(data)
}

#' Read all single cell features for a single plate
#'
#' Read all available single cell feature data for a plate given by a
#' PlateLocation object
#'
#' @param plate PlateLocation object corresponding to the plate of interest
#'              
#' @return A nested list; for structure see examples section
#'
#' @examples
#' plate <- PlateLocation("J101-2C")
#' data  <- importSingleCellDataPlate(plate)
#'
#' \dontrun{
#'   --> structure of outplut list a PlateData object
#'     |--> metadata for plate
#'     |--> data (list of 384 WellData objects)
#'        |--> A1
#'       ...
#'        |--> P24
#'           |--> metadata for well
#'           |--> data (list of 9 ImageData objects)
#'              |--> img_11
#'             ...
#'              |--> img_33
#' }
#'
importSingleCellDataPlate <- function(plate) {
  
  readFeatureFile <- function(filepath) {
    data <- tryCatch({
      # nesting of list contains no information
      readMat(filepath)[[1]][[1]][[1]][[1]]
    }, error = function(err) {
      # catch errors s.t. the rest of the batch import can continue
      warning("an error reading the file", filepath, "occurred:\n", err)
      return(NULL)
    })
    return(data)
  }
  
  getFeatureName <- function(filepath) {
    # extract the filename from the path
    name <- tail(unlist(strsplit(filepath, "[/]")), n=1)
    # split into substrings separated by "."
    name <- unlist(strsplit(name, "[.]"))
    # remove last substring, hopefully ".mat"
    name <- name[-length(name)]
    name <- paste(name, collapse=".")
    return(name)
  }
  
  # input validation
  if(!any(class(plate) == "PlateLocation")) {
    stop("can only work with PlateLocation objects")
  }
  if(file.exists(getScRdsCacheFilename(plate))) {
    data <- readRDSMC(getScRdsCacheFilename(plate))
  } else {
    # list all files ending in .mat below the input dir
    feature.path <- paste0(getLocalPath(plate),
                           "/HCS_ANALYSIS_CELL_FEATURES_CC_MAT")
    filenames <- list.files(path = feature.path, pattern="\\.mat$",
                            full.names=TRUE, recursive=TRUE)
    if(length(filenames) == 0) stop("no .mat files found.")
    # get the data
    n.cores <- detectCores()-2
    registerDoMC(cores=n.cores)
    message("fetching ", length(filenames), " files; using ", n.cores,
            " cores.")
    data <- foreach(i=1:length(filenames)) %doparMC% {
      readFeatureFile(filenames[i])
    }
    message(paste(data$out, collapse="\n"))
    data <- data$res
    
    import.name <- lapply(filenames, getFeatureName)
    names(data) <- import.name
    
    # flatten imported list hierarchy by one unneeded level
    data <- lapply(data, function(feature) {
      lapply(feature, function(image) return(image[[1]]))
    })
    
    message("successfully read ", length(data), " features.")
    saveRDSMC(data, getScRdsCacheFilename(plate))
    unlink(feature.path, recursive=TRUE)
    if(dir.exists(feature.path)) 
      unlink(feature.path, recursive=TRUE, force=TRUE)
  }
  
  # remove features that are NULL
  ind.rem <- which(sapply(data, function(x) length(x) == 0))
  if(length(ind.rem) > 0) {
    message("removing ", length(ind.rem), " features (length == 0):\n", 
            paste(names(data)[ind.rem], collapse="\n"))
    data <- data[-ind.rem] 
  }
  # remove features contain only 1 entry
  ind.rem <- which(sapply(data, function(x) length(x) == 1))
  if(length(ind.rem) > 0) {
    message("removing ", length(ind.rem), " features (length == 1):\n", 
          paste(names(data)[ind.rem], collapse="\n"))
    data <- data[-ind.rem]
  }
  # add zeros to the end of features which are shorter than 3456
  too.short <- which(sapply(data, function(x) length(x) < 3456))
  if(length(too.short) > 0) {
    message("adding zeros to the end of ", length(too.short), " features:\n", 
          paste(names(data)[too.short], collapse="\n"))
    data <- lapply(data, function(x) {
      if(length(x) == 3456) return(x)
      else return(c(x, lapply(1:(3456-length(x)), function(x) return(0))))
    })
  }
  # check if all features ok
  problems <- which(sapply(data, function(x) length(x) != 3456))
  if(length(problems) > 0) {
    stop(length(problems), " features still are of incorrect length:\n", 
         paste(names(data)[problems], collapse="\n")) 
  }
  
  # for each feature, return the total number of objects or 0 if it's a list of
  # lists
  feature.length <- sapply(data, function(feature) {
    leng <- sum(sapply(feature, length))
    list <- any(sapply(feature, is.list))
    return(leng*(!list))
  })
  # set up groups (defined by equal total length)
  groups <- unique(feature.length)
  # sort all feature.lengths into those groups
  groups <- lapply(groups, function(length, feat.lengths, feat.names) {
    indices <- which(feat.lengths==length)
    length <- length
    return(list(indices=indices, length=length))
  }, feature.length, names(data))
  # get a name for each group by finding the most frequent object
  names(groups) <- lapply(groups, function(group, feat.lengths, feat.names) {
    all.names <- feat.names[feat.lengths==group$length]
    all.names <- sapply(all.names, function(name) {
      res <- unlist(strsplit(name, "[.]"))
      if(length(res) != 2) warning("unexpected feature name format")
      return(res[1])
    })
    counts <- table(all.names)
    return(names(counts)[which.max(counts)])
  }, feature.length, names(data))
  # reorder groups into lists corresponding to how the WellData objects will be
  # constructed
  # the vec group consists of data that has a single value per image, those 
  # single values are merged into a vector per image
  groups.vec <- lapply(groups, function(group) {
    if(group$length == 3456) return(group$indices)
  })
  # the mat group consists of data that is made up of a vector per image, those
  # vectors are sorted into groups of equal length an those groups in turn are
  # saved as matrices
  groups.mat <- lapply(groups, function(group) {
    if(group$length != 0 & group$length != 3456) return(group$indices)
  })
  # the lst group consist of data that is made up of lists (possibly nested) per
  # image
  groups.lst <- lapply(groups, function(group) {
    if(group$length == 0) return(group$indices)
  })
  groups <- list(mat=groups.mat[-(which(sapply(groups.mat,is.null)))],
                 vec=groups.vec[-(which(sapply(groups.vec,is.null)))],
                 lst=groups.lst[-(which(sapply(groups.lst,is.null)))])
  data <- lapply(groups, function(group.type, data) {
    return(lapply(group.type, function(group, data) {
      return(lapply(group, function(member, data) {
        return(data[[member]])        
      }, data))
    }, data))
  }, data)
  # 3456 ImageData objects are built from the complete plate data
  data <- lapply(1:3456, function(ind, data, name) {
    vec <- lapply(data$vec, function(group, ind) {
      return(unlist(lapply(group, function(feature, i) {
        return(feature[i])
      }, ind)))
    }, ind)
    mat <- lapply(data$mat, function(group, ind) {
      return(vapply(group, function(feature, i) {
        return(unlist(feature[i]))
      }, double(length(group[[1]][[ind]])), ind))
    }, ind)
    lst <- lapply(data$lst, function(group, ind) {
      return(lapply(group, function(feature, i) {
        return(feature[i])
      }, ind))
    }, ind)
    return(ImageData(name, ind, vec, mat, lst))
  }, data, getBarcode(plate))

  # construct a matrix with rows corresponding to all wells
  all.wells  <- cbind(rep(getBarcode(plate), 384), 
                      rep(LETTERS[1:16], each=24), rep(1:24, 16))
  img.names  <- c("img_11", "img_12", "img_13", "img_21", "img_22", "img_23",
                  "img_31", "img_32", "img_33")
  # get all well metadata and build 384 WellData objects from the 3456 ImageData
  # objects
  data <- lapply(1:384, function(ind, wells, data) {
    well <- WellLocation(wells[ind,1], wells[ind,2], wells[ind,3])
    imgs <- data[(ind-1)*9+1:9]
    names(imgs) <- img.names
    return(WellData(well, imgs))
  }, all.wells, data)
  names(data) <- paste0(rep(LETTERS[1:16], each=24), rep(1:24, 16))
  # return a PlateData object
  return(PlateData(plate, data))
}