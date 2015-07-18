#' Move feature(s) to different group
#'
#' A specified feature (or a vector of features) is moved to the desired group.
#' If the target group has more rows than the source group, the features are
#' replicated according to the respective pattern and if the target group has
#' fewer rows, the feature is summarized with the supplied function.
#'
#' @param data      The data (the result of a call to meltData) to be acted
#'                  upon.
#' @param from      The feature name (or vector of feature names) to be added. 
#' @param to        The group to add them to.
#' @param aggregate If the target group is smaller, this function is used to
#'                  aggregate the data.
#' 
#' @return The data.frame specified with the parameter to with the chosen
#'         features added.
#' 
#' @examples
#' data   <- PlateData(PlateLocation("J107-2C"))
#' molten <- meltData(data)
#' test   <- moveFeatures(molten, from="Cells.AreaShape_Area", to="Image")
#' 
#' @export
moveFeatures <- function(data, from, to, aggregate="mean") {
  df.to <- lapply(data, function(types, to) {
    if(to %in% names(types)) {
      return(types)
    } else {
      return(NULL)
    }
  }, to)
  not.null <- which(!sapply(df.to, is.null))[1]
  df.to  <- df.to[[not.null]]
  df.to  <- df.to[[match(to, names(df.to))[1]]]
  if("Image.Index_Plate" %in% names(df.to)) {
    image <- tabulate(df.to$Image.Index_Plate, nbins=3456)
  } else image <- NULL
  if("Well.Index" %in% names(df.to)) {
    well <- tabulate(df.to$Well.Index, nbins=384)
  } else well <- NULL
  if("Plate.Barcode" %in% names(df.to)) {
    plate <- table(df.to$Plate.Barcode)
  } else plate <- NULL
  if(is.null(image) & is.null(well) & is.null(plate)) {
    stop("need one of \"Image.Index_Plate\", \"Well.Index\", or ",
         "\"Plate.Barcode\" in \"from\" data.frame.")
  }
  add <- lapply(from, function(feat, dat, dest, pattern) {
    source <- lapply(dat, function(types, feat) {
      res <- lapply(types, function(group, feat) {
        if(feat %in% names(group)) {
          return(group)
        } else {
          return(NULL)
        }
      }, feat)
      not.null <- which(!sapply(res, is.null))[1]
      return(res[[not.null]])
    }, feat)
    not.null <- which(!sapply(source, is.null))[1]
    source   <- source[[not.null]]
    if("Image.Index_Plate" %in% names(source) & !is.null(pattern[[1]])) {
      test    <- "Image.Index_Plate"
      pattern <- pattern[[1]]
    } else if("Well.Index" %in% names(source) & !is.null(pattern[[2]])) {
      test    <- "Well.Index"
      pattern <- pattern[[2]]
    } else if("Plate.Barcode" %in% names(source) & !is.null(pattern[[3]])) {
      test    <- "Plate.Barcode"
      pattern <- pattern[[3]]
    } else {
      stop("need one of \"Image.Index_Plate\", \"Well.Index\", or ",
           "\"Plate.Barcode\" in \"to\" data.frame.")
    }
    if(nrow(dest) > nrow(source)) {
      if(!all(rep(source[[test]], pattern) == dest[[test]])) {
        stop("source and dest ordering are not identical and enforcing this",
             "\n  is currently not implemented (>).")
      }
      res      <- rep(source[[feat]], pattern)
      dim(res) <- c(length(res), 1)
      return(res)
    } else if(nrow(dest) == nrow(source)) {
      if(!all(source[[test]] == dest[[test]])) {
        stop("source and dest ordering are not identical and enforcing this",
             "\n  is currently not implemented (==).")
      }
      res      <- source[[feat]]
      dim(res) <- c(length(res), 1)
      return(res)
    } else if(nrow(dest) < nrow(source)) {
      if(test == "Image.Index_Plate") {
        pattern <- tabulate(source$Image.Index_Plate, nbins=3456)
      } else if(test == "Well.Index") {
        pattern <- tabulate(source$Well.Index, nbins=384)
      } else if(test == "Plate.Barcode") {
        pattern <- table(source$Plate.Barcode)
      } else stop("?!?")
      pattern <- cbind(c(1, cumsum(pattern[-length(pattern)]) + 1),
                       cumsum(pattern))
      test.dat <- apply(pattern, 1, function(bounds, data) {
        if(bounds[1] >= bounds[2]) return(NA)
        else return(unique(data[bounds[1]:bounds[2]]))
      }, source[[test]])
      test.dat[which(is.na(test.dat))] <- which(is.na(test.dat))
      if(!all(dest[[test]] == test.dat)) {
        stop("source and dest ordering are not identical and enforcing this",
             "\n  is currently not implemented (<).")
      }
      fun <- get(aggregate, mode="function")
      res <- apply(pattern, 1, function(bounds, data, func) {
        if(bounds[1] >= bounds[2]) return(NA)
        else return(func(data[bounds[1]:bounds[2]]))
      }, source[[feat]], fun)
      dim(res) <- c(length(res), 1)
      return(res)
    } else stop("?!?")
  }, data, df.to, list(image, well, plate))
  add <- do.call(cbind, add)
  res <- cbind(df.to, add, stringsAsFactors=FALSE)
  names(res) <- c(names(df.to), from)
  return(res)
}

#' Combine all data in a given Data object
#'
#' Given a PlateData/WellData/ImageData object, rbind all data contained in the
#' hierarchy below.
#'
#' @param x The MatData/PlateData/WellData/ImageData object of interest.
#' 
#' @return A list of data frames for the combined objects.
#' 
#' @examples
#' plate    <- PlateData(PlateLocation("J101-2C"))
#' plate.df <- meltData(plate)
#' 
#' @export
meltData <- function(x) {
  UseMethod("meltData", x)
}

#' @export
meltData.PlateData <- function(x) {
  grps.vec <- names(x$data[[1]]$data[[1]]$data.vec)
  grps.vec <- c(grps.vec, "Well")
  grps.mat <- names(x$data[[1]]$data[[1]]$data.mat)
  grps.lst <- names(x$data[[1]]$data[[1]]$data.lst)
  progress.bar <- getOption("singleCellFeatures.progressBars")
  if(progress.bar != "none") {
    message("extracting features from wells/images:")
  }
  all      <- unlist(llply(x$data, meltData, .progress=progress.bar),
                     recursive=FALSE)
  if(progress.bar != "none") {
    message("building vec feature data frame.")
  }
  plate    <- data.frame("Plate.Barcode"=x$meta$plate.barcode,
                         "Plate.Quality"=x$meta$plate.quality,
                         "Plate.Type"=x$meta$plate.type,
                         "Experiment.Name"=x$meta$experiment.name,
                         "Experiment.Replicate"=x$meta$experiment.replicate,
                         "Experiment.Batch"=x$meta$experiment.batch,
                         stringsAsFactors=FALSE)
  res.vec  <- llply(grps.vec, function(group, data) {
    regexp <- paste0("^[A-P]([1-9]|1[0-9]|2[0-4])\\.vec\\.", group, "$")
    res <- data[grep(regexp, names(data))]
    res <- rbind.fill(res)
    return(res)
  }, unlist(all[grep("^[A-P]([1-9]|1[0-9]|2[0-4])\\.vec$", names(all))],
            recursive=FALSE))
  res.vec  <- c(res.vec, list(plate))
  if(progress.bar != "none") {
    message("building mat feature data frame.")
  }
  res.mat  <- llply(grps.mat, function(group, data) {
    regexp <- paste0("^[A-P]([1-9]|1[0-9]|2[0-4])\\.mat\\.", group, "$")
    res <- data[grep(regexp, names(data))]
    res <- rbind.fill(res)
    return(res)
  }, unlist(all[grep("^[A-P]([1-9]|1[0-9]|2[0-4])\\.mat$", names(all))],
            recursive=FALSE))
  if(progress.bar != "none") {
    message("building lst feature matrices.")
  }
  res.lst  <- llply(grps.lst, function(group, data) {
    regexp <- paste0("^[A-P]([1-9]|1[0-9]|2[0-4])\\.lst\\.", group, "$")
    res <- data[grep(regexp, names(data))]
    if(group == "OtherFeatures") {
      names(res) <- sapply(names(res), function(x) {
        return(unlist(strsplit(x, "[.]"))[1])
      })
    } else {
      fnames <- names(res[[1]])
      res <- lapply(fnames, function(name, data) {
        regexp <- paste0("^[A-P]([1-9]|1[0-9]|2[0-4])\\.lst\\..*", name, "$")
        dat <- lapply(data[grep(regexp, names(data))], function(x) {
          if(is.null(x)) return(matrix(nrow=0, ncol=0))
          else return(x)
        })
        result <- bdiag(dat)
        dims <- sapply(dat, function(x) return(dimnames(x)[1]))
        dimname <- unlist(lapply(1:length(dims), function(i, d) {
          if(length(d[[i]]) > 0) {
            wellname <- getWellIndex2D(i)
            return(paste0(wellname$wel.row, wellname$wel.col, "_", d[[i]]))
          } else return(NULL)
        }, dims))
        dimnames(result) <- list(dimname, NULL)
        return(result)
      }, unlist(res, recursive=FALSE))
      names(res) <- fnames
    }
    return(res)
  }, unlist(all[grep("^[A-P]([1-9]|1[0-9]|2[0-4])\\.lst$", names(all))],
            recursive=FALSE))
  names(res.vec) <- c(grps.vec, "Plate")
  names(res.mat) <- grps.mat
  names(res.lst) <- grps.lst
  return(list(vec=res.vec, mat=res.mat, lst=res.lst))
}

#' @export
meltData.WellData <- function(x) {
  grps.vec <- names(x$data[[1]]$data.vec)
  grps.mat <- names(x$data[[1]]$data.mat)
  grps.lst <- names(x$data[[1]]$data.lst)
  all      <- unlist(lapply(x$data, meltData), recursive=FALSE)
  well     <- data.frame("Plate.Barcode"=x$meta$plate.barcode,
                         "Well.Index"=x$meta$well.index,
                         "Well.Row"=x$meta$well.row,
                         "Well.Col"=x$meta$well.col,
                         "Well.Count_Cells"=x$meta$counts.cells,
                         "Well.Count_Pathogens"=x$meta$counts.pathogen,
                         "Well.Count_Infected"=x$meta$counts.infection,
                         "Well.Type"=x$meta$well.type,
                         "Well.Quality"=x$meta$well.quality,
                         "Well.Gene_Name"=x$meta$gene.name,
                         "Well.Gene_ID"=x$meta$gene.id,
                         "Well.siRNA_Name"=x$meta$sirna.name,
                         "Well.siRNA_Sequence"=x$meta$sirna.sequence,
                         "Well.siRNA_Seed"=x$meta$sirna.seed,
                         "Well.siRNA_Target"=x$meta$sirna.target,
                         stringsAsFactors=FALSE)
  res.vec  <- lapply(grps.vec, function(group, data) {
    regexp <- paste0("^img_[1-3]{2}\\.vec\\.", group, "$")
    res <- data[grep(regexp, names(data))]
    res <- do.call(rbind, res)
    return(res)
  }, unlist(all[grep("^img_[1-3]{2}\\.vec$", names(all))], recursive=FALSE))
  res.vec  <- c(res.vec, list(well))
  res.mat  <- lapply(grps.mat, function(group, data) {
    regexp <- paste0("^img_[1-3]{2}\\.mat\\.", group, "$")
    res <- data[grep(regexp, names(data))]
    res <- do.call(rbind, res)
    return(res)
  }, unlist(all[grep("^img_[1-3]{2}\\.mat$", names(all))], recursive=FALSE))
  res.lst  <- lapply(grps.lst, function(group, data) {
    regexp <- paste0("^img_[1-3]{2}\\.lst\\.", group, "$")
    res <- data[grep(regexp, names(data))]
    if(group == "OtherFeatures") {
      names(res) <- sapply(names(res), function(x) {
        return(unlist(strsplit(x, "[.]"))[1])
      })
    } else {
      fnames <- names(res[[1]])
      res <- lapply(fnames, function(name, data) {
        regexp <- paste0("^img_[1-3]{2}\\.lst\\..*", name, "$")
        dat <- lapply(data[grep(regexp, names(data))], function(x) {
          if(is.null(x)) return(matrix(nrow=0, ncol=0))
          else return(x)
        })
        result <- bdiag(dat)
        dims <- sapply(dat, function(x) return(dim(x)[1]))
        dimname <- unlist(lapply(1:length(dims), function(i, d) {
          return(rep(i, d[i]))
        }, dims))
        dimnames(result) <- list(dimname, NULL)
        return(result)
      }, unlist(res, recursive=FALSE))
      names(res) <- fnames
    }
    return(res)
  }, unlist(all[grep("^img_[1-3]{2}\\.lst$", names(all))], recursive=FALSE))
  names(res.vec) <- c(grps.vec, "Well")
  names(res.mat) <- grps.mat
  names(res.lst) <- grps.lst
  return(list(vec=res.vec, mat=res.mat, lst=res.lst))
}

#' @export
meltData.ImageData <- function(x) {
  well.index  <- getWellIndex1D(x$well.row, x$well.col, NULL, x$image.total)
  image.index <- getWellIndex1D(x$well.row, x$well.col, x$image.index,
                                x$image.total)
  well.name   <- paste0(x$well.row, x$well.col)
  
  res.vec <- lapply(x$data.vec, function(group) {
    if(is.null(group)) {
      return(NULL)
    } else {
      if (nrow(group) > 0) {
        newcols <- data.frame(x$image.index, image.index, well.index,
                              well.name, x$plate, stringsAsFactors=FALSE)
        colnames(newcols) <- c("Image.Index", "Image.Index_Plate",
                               "Well.Index", "Well.Name", "Plate.Barcode")
        return(cbind(group, newcols))
      } else {
        return(NULL)
      }
    }
  })

  res.mat <- lapply(x$data.mat, function(group) {
    if(is.null(group)) {
      return(NULL)
    } else {
      n.rows <- nrow(group)
      if(n.rows > 0) {
        image.ind        <- rep(x$image.index, n.rows)
        dim(image.ind)   <- c(n.rows, 1)
        imgwell.ind      <- rep(image.index, n.rows)
        dim(imgwell.ind) <- c(n.rows, 1)
        well.ind         <- rep(well.index, n.rows)
        dim(well.ind)    <- c(n.rows, 1)
        well.nme         <- rep(well.name, n.rows)
        dim(well.nme)    <- c(n.rows, 1)
        plate            <- rep(x$plate, n.rows)
        dim(plate)       <- c(n.rows, 1)

        info             <- data.frame(image.ind, imgwell.ind, well.ind,
                                       well.nme, plate, stringsAsFactors=FALSE)
        colnames(info)   <- c("Image.Index", "Image.Index_Plate", "Well.Index",
                              "Well.Name", "Plate.Barcode")

        return(cbind(group, info))
      } else {
        return(NULL)
      }
    }
  })
  return(list(vec=res.vec, mat=res.mat, lst=x$data.lst))
}

#' @export
meltData.default <- function(x) {
  stop("can only deal with Data (ImageData/WellData/PlateData) objects.")
}