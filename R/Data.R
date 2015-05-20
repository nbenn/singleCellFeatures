#' Constructor for MatData objects
#'
#' Given a PlateLocation object, fetch (from file or from openBIS) the
#' corresponding MatData object.
#'
#' @param plate  PlateLocation object corresponding to the plate of interest
#'              
#' @return A nested list; for structure see examples section
#'
#' @examples
#' plate <- PlateLocation("J101-2C")
#' data  <- MatData(plate)
#'
#' \dontrun{
#'   --> structure of output
#'     |--> meta
#'        |--> metadata for plate
#'     |--> data
#'        |--> Feature_1
#'       ...
#'        |--> Feature_3456
#' }
#' 
#' @export
MatData <- function(plate, force.download=FALSE) {

  # input validation
  if (!any(class(plate) == "PlateLocation")) {
    stop("can only work with a single PlateLocation object")
  }
  # check if plate is cached
  if (file.exists(getCacheFilenameData(plate)) & !force.download) {
    data <- readRDSMC(getCacheFilenameData(plate))
    result <- structure(list(meta = PlateMetadata(plate),
                             data = data), class = c("MatData", "Data"))
  } else {
    data      <- dowloadFeatureHelper(plate, "cc", force.download)
    infection <- dowloadFeatureHelper(plate, "dectree", force.download)

    if(is.null(infection)) {
      infection <- dowloadFeatureHelper(plate, "thresholded", force.download)
    }
    if(is.null(infection)) {
      stop("try something other than 'dectree'/'thresholded'")
    } else {
      names(infection) <- "Cells.Infection_IsInfected"
      data <- c(data, infection)
    }

    result <- structure(list(meta = PlateMetadata(plate),
                             data = data), class = c("MatData", "Data"))
    saveToCache(result, force.write=force.download)
    dirs <- dir(getLocalPath(plate), "^HCS_ANALYSIS_CELL_", full.names=TRUE)
    unlink(dirs, recursive=TRUE)
  }

  # remove features that have zero length
  ind.rem <- which(sapply(result$data, function(x) length(x) == 0))
  if (length(ind.rem) > 0) {
    message("removing ", length(ind.rem), " features (length == 0):\n",
            paste(names(result$data)[ind.rem], collapse="\n"))
    result$data <- result$data[-ind.rem]
  }
  # check if feature set is complete
  checkCompletenessFeature(result)

  return(result)
}

dowloadFeatureHelper <- function(plate, type, force.download=FALSE,
                                 features=".*.mat") {
  # input validation
  if (!any(class(plate) == "PlateLocation")) {
    stop("can only work with a single PlateLocation object")
  }
  if (!type %in% c("cc", "dectree", "thresholded")) {
    stop("unrecognized feature type")
  }
  # plate has to be downloaded/imported
  plate.path <- getLocalPath(plate)
  if(type == "cc") {
    feature.path <- paste0(plate.path, "/", "HCS_ANALYSIS_CELL_FEATURES_CC_MAT")
  } else if(type == "dectree") {
    feature.path <- paste0(plate.path, "/",
                           "HCS_ANALYSIS_CELL_DECTREECLASSIFIER_MAT")
  } else if(type == "thresholded") {
    feature.path <- paste0(plate.path, "/",
                           "HCS_ANALYSIS_CELL_THRESHOLDEDINFECTIONSCORING_MAT")
  } else stop("unrecognized feature type")
  # search for matlab files
  filenames <- list.files(path = feature.path, pattern="\\.mat$",
                          full.names=TRUE, recursive=TRUE)
  if (length(filenames) == 0 | force.download) {
    # plate has to be downloaded
    message("downloading '", type, "' for plate ", getBarcode(plate),
            " from openBIS.")
    config <- configGet()
    bee.command <- paste0(
      "cd ", config$beeDownloader$executable, "; export BEESOFTSRC='",
      config$beeDownloader$beeSoftsrc, "'; ./BeeDataSetDownloader.sh ",
      "--user '", config$openBIS$username, "' --password '",
      config$openBIS$password, "' --outputdir '", config$dataStorage$dataDir,
      "' --plateid '^", getOpenBisPath(plate), "' --files '", features, "'",
      " --verbose '10'")
    if(type == "dectree") {
      bee.command <- paste0(bee.command,
                            " --type 'HCS_ANALYSIS_CELL_DECTREECLASSIFIER_MAT'")
    } else if(type == "thresholded") {
      bee.command <- paste0(bee.command, " --type 'HCS_ANALYSIS_CELL_",
                            "THRESHOLDEDINFECTIONSCORING_MAT'")
    }
    system(bee.command, ignore.stdout = TRUE)
  }
  # list all files ending in .mat below the input dir
  filenames <- list.files(path=feature.path, pattern="\\.mat$",
                          full.names=TRUE, recursive=TRUE)

  if(type != "cc" & length(filenames) != 1) {
    message("a total of ", length(filenames), " .mat files found. expecting 1.")
    return(NULL)
  } else if (length(filenames) == 0) stop("no .mat files found.")

  res <- readMatFeatureHelper(feature.path)
  return(res)
}

readMatFeatureHelper <- function(path) {

  readFeatureFile <- function(filepath) {
    data <- tryCatch({
      # nesting of list contains no information
      readMat(filepath)[[1]][[1]][[1]][[1]]
    },
    error = function(err) {
      # catch errors s.t. the rest of the batch import can continue
      warning("an error reading the file", filepath, "occurred:\n", err)
      errfile <- unlist(strsplit(filepath, "/"))
      errfile[length(errfile)-2] <- "Errors"
      errfile <- paste0(errfile[-(length(errfile)-1)],
                       collapse="/")
      errdir  <- dirname(errfile)
      if(!dir.exists(errdir)) dir.create(errdir, recursive=TRUE)
      file.rename(filepath, errfile)
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

  # list all files ending in .mat below the input dir
  filenames <- list.files(path=path, pattern="\\.mat$",
                          full.names=TRUE, recursive=TRUE)
  if (length(filenames) == 0) stop("no .mat files found.")
  # get the data
  message("fetching ", length(filenames), " files")
  data <- lapply(filenames, readFeatureFile)
  import.name <- lapply(filenames, getFeatureName)
  names(data) <- import.name

  # flatten imported list hierarchy by one unneeded level
  data <- lapply(data, function(feature) {
    lapply(feature, function(image) return(image[[1]]))
  })

  message("successfully read ", length(data), " features.")
  return(data)
}

#' Constructor for PlateData objects
#'
#' Given a PlateLocation object, fetch and/or restructure the corresponding
#' MatData object into a PlateLocation object.
#'
#' @param plate  PlateLocation object corresponding to the plate of interest
#' @param data   MatData object corresponding to the plate of interest
#' @param select A vector of strings that are matched with features to be kept.
#' @param drop   A vector of strings that are matched with features to be
#'               dropped.
#'              
#' @return A nested list; for structure see examples section
#'
#' @examples
#' plate <- PlateLocation("J101-2C")
#' data  <- PlateData(plate)
#'
#' \dontrun{
#'   --> structure of output
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
#' @export
PlateData <- function(plate, select=NULL, drop=NULL, data=NULL) {
  # input validation
  if (!any(class(plate) == "PlateLocation")) {
    stop("can only work with a single PlateLocation object")
  }
  # get data
  if (is.null(data)) {
    data <- MatData(plate)
  } else {
    if (!any(class(data) == "MatData")) {
      stop("can only work with a single MatData object")
    }
  }
  # determine if wells contain 9 or 6 images
  n.imgs <- getNoImgPerWell(data)
  tot.nimgs <- n.imgs * 384
  # remove specified features
  if (!is.null(select) | !is.null(drop)) {
    data <- extractFeatures(data, select, drop)
  }
  # remove features contain only 1 entry
  ind.rem <- which(sapply(data$data, function(x) length(x) == 1))
  if (length(ind.rem) > 0) {
    message("removing ", length(ind.rem), " features (length == 1):\n",
          paste(names(data$data)[ind.rem], collapse="\n"))
    data$data <- data$data[-ind.rem]
  }
  # add zeros to the end of features which are shorter than 2304/3456
  too.short <- which(sapply(data$data, function(x) length(x) < tot.nimgs))
  if (length(too.short) > 0) {
    message("adding zeros to the end of ", length(too.short), " features:")
    data$data[too.short] <- lapply(too.short, function(ind, dat) {
      n.zero <- tot.nimgs - length(dat[[ind]])
      message("  adding ", n.zero, " zeros to ", names(dat[ind]))
      return(c(dat[[ind]], lapply(1:n.zero, function(x) return(0))))
    },data$data)
  }
  # check if all features ok
  problems <- which(sapply(data$data, function(x) length(x) != tot.nimgs))
  if (length(problems) > 0) {
    stop(length(problems), " features still are of incorrect length:\n",
         paste(names(data$data)[problems], collapse="\n"))
  }

  # for each feature, return the total number of objects or 0 if it's a list of
  # lists
  feature.length <- sapply(data$data, function(feature) {
    leng <- sum(sapply(feature, length))
    list <- any(sapply(feature, is.list))
    return(leng * (!list))
  })
  # set up groups (defined by equal total length)
  groups <- unique(feature.length)
  # sort all feature.lengths into those groups
  groups <- lapply(groups, function(length, feat.lengths, feat.names) {
    indices <- which(feat.lengths == length)
    length <- length
    return(list(indices=indices, length=length))
  }, feature.length, names(data$data))
  # get a name for each group by finding the most frequent object
  names(groups) <- lapply(groups, function(group, feat.lengths, feat.names) {
    all.names <- feat.names[feat.lengths == group$length]
    all.names <- sapply(all.names, function(name) {
      res <- unlist(strsplit(name, "[.]"))
      if (length(res) != 2) warning("unexpected feature name format")
      return(res[1])
    })
    counts <- table(all.names)
    return(names(counts)[which.max(counts)])
  }, feature.length, names(data$data))
  # reorder groups into lists corresponding to how the WellData objects will be
  # constructed
  # the vec group consists of data that has a single value per image, those
  # single values are merged into a vector per image
  groups.vec <- lapply(groups, function(group) {
    if (group$length == tot.nimgs) return(group$indices)
  })
  # the mat group consists of data that is made up of a vector per image, those
  # vectors are sorted into groups of equal length an those groups in turn are
  # saved as matrices
  groups.mat <- lapply(groups, function(group) {
    if (group$length != 0 & group$length != tot.nimgs) return(group$indices)
  })
  # the lst group consist of data that is made up of lists (possibly nested)
  # per image
  groups.lst <- lapply(groups, function(group) {
    if (group$length == 0) return(group$indices)
  })
  # remove empty lists
  rem.mat <- which(sapply(groups.mat,is.null))
  rem.vec <- which(sapply(groups.vec,is.null))
  rem.lst <- which(sapply(groups.lst,is.null))
  if (length(rem.mat) > 0) groups.mat <- groups.mat[-rem.mat]
  if (length(rem.vec) > 0) groups.vec <- groups.vec[-rem.vec]
  if (length(rem.lst) > 0) groups.lst <- groups.lst[-rem.lst]
  # PercentTouchingNeighbors requires IdentityOfNeighbors to be of use
  neigh.ident <- grep("^Neighbors.IdentityOfNeighbors_", names(groups.lst[[1]]))
  neigh.touch <- grep("^Neighbors.PercentTouchingNeighbors_",
                      names(groups.lst[[1]]))
  if(length(neigh.touch) > 0) {
    if(length(neigh.ident) != length(neigh.touch)) {
      warning("having \"Neighbors.PercentTouchingNeighbors\" features",
              " without corresponding\n  \"Neighbors.IdentityOfNeighbors\"",
              " features makes little sense.")
    } else {
      l_ply(names(groups.lst[[1]][neigh.touch]), function(touch, ident) {
        object <- unlist(strsplit(touch, 
                                  "Neighbors.PercentTouchingNeighbors_"))[2]
        if(!any(grepl(paste0(object, "$"), ident))) {
          warning("having \"Neighbors.PercentTouchingNeighbors\" features",
                  " without corresponding\n  \"Neighbors.IdentityOfNeighbors\"",
                  " features makes little sense: could not\n  find ", touch)
        }
      }, names(groups.lst[[1]][neigh.ident]))      
    }
  }
  # manually modify the lst group
  if(length(neigh.ident) > 0 & length(neigh.touch) > 0) {
    rest <- groups.lst[[1]][-c(neigh.ident, neigh.touch)]
    iden <- groups.lst[[1]][neigh.ident]
    touc <- groups.lst[[1]][neigh.touch]
    if(length(rest) > 0) {
      groups.lst <- list("IdentityOfNeighbors" = iden,
                         "PercentTouchingNeighbors" = touc,
                         "OtherFeatures" = rest)
    } else {
      groups.lst <- list("IdentityOfNeighbors" = iden,
                         "PercentTouchingNeighbors" = touc)
    }
  } else if(length(neigh.ident) > 0) {
    rest <- groups.lst[[1]][-neigh.ident]
    iden <- groups.lst[[1]][neigh.ident]
    if(length(rest) > 0) {
      groups.lst <- list("IdentityOfNeighbors" = iden,
                         "OtherFeatures" = rest)
    } else {
      groups.lst <- list("IdentityOfNeighbors" = iden)
    }
  } else if(length(neigh.touch) > 0) {
    rest <- groups.lst[[1]][-neigh.touch]
    touc <- groups.lst[[1]][neigh.touch]
    if(length(rest) > 0) {
      groups.lst <- list("PercentTouchingNeighbors" = touc,
                         "OtherFeatures" = rest)
    } else {
      groups.lst <- list("PercentTouchingNeighbors" = touc)
    }
  }

  groups <- list(vec=groups.vec,
                 mat=groups.mat,
                 lst=groups.lst)

  data$data <- lapply(groups, function(group.type, data) {
    return(lapply(group.type, function(group, data) {
      return(lapply(group, function(member, data) {
        return(data[[member]])
      }, data))
    }, data))
  }, data$data)

  # the following step will need a lot of memory
  gc()
  # tot.nimgs ImageData objects are built from the complete plate data
  data$data <- llply(1:tot.nimgs, function(ind, data, name, n.imgs) {
    vec <- lapply(data$vec, function(group, ind) {
      return(unlist(lapply(group, function(feature, i) {
        return(feature[i])
      }, ind)))
    }, ind)
    mat <- lapply(data$mat, function(group, ind) {
      n.rows <- length(group[[1]][[ind]])
      grp <- vapply(group, function(feature, i) {
        return(unlist(feature[i]))
      }, double(n.rows), ind)
      names <- colnames(grp)
      if(is.null(names)) names <- names(grp)
      n.cols <- length(names)
      dim(grp) <- c(n.rows, n.cols)
      colnames(grp) <- names
      rownames(grp) <- NULL
      return(grp)
    }, ind)
    lst <- mapply(function(group, gname, ind) {
      if(gname == "IdentityOfNeighbors") {
        return(lapply(group, function(feature, i) {
          # build sparse adjacency matrices
          l <- length(feature[[i]])
          if(l > 1) {
            p <- c(0, cumsum(sapply(feature[[i]], function(x) length(x[[1]]))))
            j <- unlist(feature[[i]])
            return(sparseMatrix(j=j, p=p, dims=c(l, l)))
          } else return(NULL)
        }, ind))
      } else if(gname == "PercentTouchingNeighbors") {
        return(lapply(group, function(feature, i) return(unlist(feature[[i]])),
                      ind))
      } else {
        return(lapply(group, function(feature, i) return(feature[i]), ind))
      }
    }, data$lst, names(data$lst), list(ind=ind), SIMPLIFY = FALSE)
    if(!is.null(lst$PercentTouchingNeighbors) & 
       !is.null(lst$IdentityOfNeighbors)) {
      lst$PercentTouchingNeighbors <- mapply(function(feat, fname, mats) {
        object <- unlist(strsplit(fname, 
                                  "Neighbors.PercentTouchingNeighbors_"))[2]
        mat <- mats[[grep(paste0("Neighbors.IdentityOfNeighbors_", object),
                          names(mats))]]
        if(is.null(mat)) {
          return(feat)
        } else {
          return(sparseMatrix(j=mat@i, p=mat@p, x=feat, dims=dim(mat),
                 index1=FALSE))
        }
      }, lst$PercentTouchingNeighbors, names(lst$PercentTouchingNeighbors),
         list(mats=lst$IdentityOfNeighbors), SIMPLIFY = FALSE)
    }
    return(ImageData(name, ind, n.imgs, vec, mat, lst))
  }, data$data, getBarcode(plate), n.imgs, .progress = "text")
  # free no longer needed memory
  gc()

  # construct a matrix with rows corresponding to all wells
  all.wells  <- cbind(rep(getBarcode(plate), 384),
                      rep(LETTERS[1:16], each=24), rep(1:24, 16))
  if(n.imgs == 9) {
    img.names  <- c("img_11", "img_12", "img_13", "img_21", "img_22", "img_23",
                    "img_31", "img_32", "img_33")
  } else {
    img.names  <- c("img_11", "img_12", "img_21", "img_22", "img_31", "img_32")
  }
  # get all well metadata and build 384 WellData objects from the 3456
  # ImageData objects
  data$data <- lapply(1:384,
    function(ind, wells, data) {
      well <- WellLocation(wells[ind,1], wells[ind,2], wells[ind,3])
      imgs <- data[(ind - 1) * n.imgs + 1:n.imgs]
      names(imgs) <- img.names
      return(WellData(well, imgs))
    },
  all.wells, data$data)
  names(data$data) <- paste0(rep(LETTERS[1:16], each=24), rep(1:24, 16))
  # return a PlateData object
  return(structure(data, class = c("PlateData", "Data")))
}

#' Constructor for WellData objects
#' 
#' Given a WellLocation object, and a list of ImageData objects, a WellData
#' object is built up.
#'
#' @param well A WellLocation object for the target well
#' @param data An optional list of ImageData objects corresponding to the
#'             target well
#'
#' @return A WellData object: a list with slots for metadata and data, with the
#'         data slot holding upt to 9 ImageData objects
#'
#' @examples
#' data <- WellData(WellLocation("J107-2C", "D", 15))
#' 
#' \dontrun{
#'   --> structure of output
#'     |--> metadata for well
#'     |--> data (list of 6/9 ImageData objects)
#'       |--> img_11
#'      ...
#'       |--> img_33/32
#' }
#'
#' @export
WellData <- function(well, data=NULL) {
  if(is.null(data)) {
    if(file.exists(getCacheFilenameData(well))) {
      result <- readRDS(getCacheFilenameData(well))
    } else {
      data   <- PlateData(convertToPlateLocation(well))
      result <- extractWells(data, well, FALSE)
      saveToCache(result)
    }
  } else {
    if(!is.list(data)) stop("expecting a LIST of ImageData objects")
    if(length(data) < 1 | length(data) > 9) {
      stop("expecting a list of ImageData objects with length in 1:9")
    }
    if(!all(sapply(data, function(x) any(class(x) == "ImageData")))) {
      stop("expecting a list of IMAGEDATA objects")
    }
    result <- structure(list(meta = WellMetadata(well),
                             data = data), class = c("WellData", "Data"))
  }

  return(result)
}

#' Constructor for ImageData objects
#' 
#' Due to performance requirements (when a whole plate is imported, 3456
#' ImageData objects are created), there are two ways this function is executed:
#' Either the data is supplied (at leas one of vec, mat, lst is not NULL), in
#' which case barcode is expected to be a string, index the linear plate index
#' and n.img in c(6, 9), the number of images per well. Nothing is checked to
#' save time.
#' If no data is supplied, it is fetched. In this case, barcode is expected to
#' be a WellLocation object, index the image index inside the well and n.img to
#' be NULL, as this has already been determined.
#'
#' @param barcode Either the barcode of the plate or a WellLocation object
#'                corresponding to the image's well
#' @param index   Either the linearized index of the image within the whole
#'                plate or the image index within the well
#' @param n.img   Number of images per well (either 6 or 9)
#' @param vec     All data than contains only a single value per feature
#' @param mat     All data which consists of a vector of values per feature
#' @param lst     All data than can only be expressed as a nested list per
#'                feature
#'
#' @return An ImageData object: a list with slots for plate barcode, well
#'         location within the plate (index and row/col), image location within
#'         the well (index, row/col), number of images per well and image data
#'
#' @examples
#' data <- ImageData(WellLocation("J107-2C", "D", 15), 5)
#' 
#' @export
ImageData <- function(barcode, index, n.img=NULL, vec=NULL, mat=NULL,
                      lst=NULL) {
  if(is.null(vec) & is.null(mat) & is.null(lst)) {
    if(!is.null(n.img)) warning("ignoring parameter n.img")
    if(!any(class(barcode) == "WellLocation")) {
      stop("in case no data is supplied, the first argument is expected to ",
           "be of type WellLocation")
    }
    data.wel <- WellData(barcode)
    data.img <- extractImages(data.wel, index, keep.well=FALSE)
    return(data.img)
  } else {
    ind <- getWellIndex2D(index, n.img)
    result <- list(plate       = barcode,
                   well.index  = index,
                   well.row    = ind$wel.row,
                   well.col    = ind$wel.col,
                   image.index = ind$img.ind,
                   image.row   = ind$img.row,
                   image.col   = ind$img.col,
                   image.total = n.img,
                   data.vec    = vec,
                   data.mat    = mat,
                   data.lst    = lst)
    return(structure(result, class = c("ImageData", "Data")))
  }
}