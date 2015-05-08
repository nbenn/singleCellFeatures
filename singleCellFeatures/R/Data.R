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
  readFeatureFile <- function(filepath) {
    data <- tryCatch({
      # nesting of list contains no information
      readMat(filepath)[[1]][[1]][[1]][[1]]
    },
    error = function(err) {
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
  if (!any(class(plate) == "PlateLocation")) {
    stop("can only work with a single PlateLocation object")
  }
  # check if plate is cached
  if (file.exists(getScRdsCacheFilename(plate))) {
    data <- readRDSMC(getScRdsCacheFilename(plate))
  } else {
    # plate has to be downloaded/imported
    plate.path <- getLocalPath(plate)
    feature.path <- paste0(plate.path, "/", "HCS_ANALYSIS_CELL_FEATURES_CC_MAT")
    # search for matlab files
    filenames <- list.files(path = feature.path, pattern="\\.mat$",
                            full.names=TRUE, recursive=TRUE)
    if (length(filenames) == 0) {
      # plate has to be downloaded
      message("downloading plate ", getBarcode(plate),
              " from openBIS.")
      data(settingsDatabase, envir = environment())
      # get password from file
      pwd.openBIS <- readRDS(settings.database$openBIS.pwd)
      bee.command <- paste0(
        "cd ", settings.database$bee.dir, "; export BEESOFTSRC='",
        settings.database$bee.softsrc, "'; ./BeeDataSetDownloader.sh ",
        "--user '", settings.database$openBIS.user, "' --password '",
        pwd.openBIS, "' --outputdir '", settings.database$openBIS.data,
        "' --plateid '^", getOpenBisPath(plate), "' --files '.*.mat'",
        " --verbose '10'")
      system(bee.command, ignore.stdout = TRUE)
    }

    # list all files ending in .mat below the input dir
    filenames <- list.files(path = feature.path, pattern="\\.mat$",
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
    saveRDSMC(data, getScRdsCacheFilename(plate))
    unlink(feature.path, recursive=TRUE)
    if (dir.exists(feature.path))
      unlink(feature.path, recursive=TRUE, force=TRUE)
  }

  result <- list(meta = PlateMetadata(plate),
                 data = data)
  return(structure(result, class = c("MatData", "Data")))
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
  lengths <- sapply(data$data, length)
  max.length <- max(lengths)
  ful.length <- lengths[lengths == max.length]
  full.fract <- length(ful.length) / length(lengths)
  if (max.length > 2304) n.imgs <- 9
  else if (full.fract > 0.5) n.imgs <- 6
  else {
    n.imgs <- 9
    warning("unsure about the number of images per well.")
  }
  message("assuming ", n.imgs, " images per well:\nmax legnth: ", max.length,
          ", fraction of full length features: ", full.fract)
  tot.nimgs <- n.imgs * 384
  # check for completenes of single cell data
  pathogen <- getPathogen(plate)
  # load feature database
  data(featureDatabase, envir=environment())
  feat.exp <- feature.database[[pathogen]]
  feat.dat <- getFeatureNames(data)
  if(is.null(feat.exp)) {
    warning("for ", pathogen, ", currently no feature list is available. ",
            "Please add one using updateDatabaseFeatures.")
  }
  superfl <- setdiff(feat.exp, feat.dat)
  missing <- setdiff(feat.dat, feat.exp)
  if(length(missing) > 0) {
    warning("found ", length(missing), " missing features (", getBarcode(plate),
            "):\n  ", paste(missing, collapse="\n  "))
  }
  if(length(superfl) > 0) {
    warning("found ", length(superfl), " superfluous features (",
            getBarcode(plate), "):\n  ", paste(superfl, collapse="\n  "))
  }
  # remove specified features
  if (!is.null(select) | !is.null(drop)) {
    data <- extractFeatures(data, select, drop)
  }
  # remove features that are NULL
  ind.rem <- which(sapply(data$data, function(x) length(x) == 0))
  if (length(ind.rem) > 0) {
    message("removing ", length(ind.rem), " features (length == 0):\n",
            paste(names(data$data)[ind.rem], collapse="\n"))
    data$data <- data$data[-ind.rem]
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
    data$data[too.short] <- lapply(
      too.short,
      function(ind, dat) {
        n.zero <- tot.nimgs - length(dat[[ind]])
        message("  adding ", n.zero, " zeros to ", names(dat[ind]))
        return(c(dat[[ind]], lapply(1:n.zero, function(x) return(0))))
      },
      data$data
    )
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
  groups <- lapply(
    groups,
    function(length, feat.lengths, feat.names) {
      indices <- which(feat.lengths == length)
      length <- length
      return(list(indices=indices, length=length))
    },
    feature.length, names(data$data)
  )
  # get a name for each group by finding the most frequent object
  names(groups) <- lapply(
    groups,
    function(group, feat.lengths, feat.names) {
      all.names <- feat.names[feat.lengths == group$length]
      all.names <- sapply(all.names, function(name) {
        res <- unlist(strsplit(name, "[.]"))
        if (length(res) != 2) warning("unexpected feature name format")
        return(res[1])
      })
      counts <- table(all.names)
      return(names(counts)[which.max(counts)])
    },
    feature.length, names(data$data)
  )
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
  groups <- list(mat=groups.mat,
                 vec=groups.vec,
                 lst=groups.lst)
  data$data <- lapply(
    groups,
    function(group.type, data) {
      return(lapply(group.type,
        function(group, data) {
          return(lapply(group,
            function(member, data) {
              return(data[[member]])
            },
          data))
        },
      data))
    },
    data$data
  )

  #the following step will need a lot of memory
  gc()
  # tot.nimgs ImageData objects are built from the complete plate data
  data$data <- lapply(
    1:tot.nimgs,
    function(ind, data, name, n.imgs) {
      vec <- lapply(data$vec,
        function(group, ind) {
          return(unlist(lapply(group,
            function(feature, i) {
              return(feature[i])
            },
          ind)))
        },
      ind)
      mat <- lapply(data$mat,
        function(group, ind) {
          return(vapply(group,
            function(feature, i) {
              return(unlist(feature[i]))
            },
          double(length(group[[1]][[ind]])), ind))
        },
      ind)
      lst <- lapply(data$lst,
        function(group, ind) {
          return(lapply(group,
            function(feature, i) {
              return(feature[i])
            },
          ind))
        },
      ind)
      return(ImageData(name, ind, n.imgs, vec, mat, lst))
    },
    data$data, getBarcode(plate), n.imgs
  )

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
#' @export
WellData <- function(well, data) {
  result <- list(meta = WellMetadata(well),
                 data = data)
  return(structure(result, class = c("WellData", "Data")))
}

#' Constructor for ImageData objects
#' 
#' @export
ImageData <- function(barcode, index, n.img, vec, mat, lst) {
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