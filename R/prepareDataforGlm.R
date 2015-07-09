#' Prepare a dataset for GLM analysis
#' 
#' Given two design matrices (one corresponding to a knockdown and one to a
#' control well), create a concatenated design matrix with a column "Response"
#' encoding for original membership. This unified data matrix is then divided
#' into 90% training and 10% testing. To avoid missingness problems, all
#' variables containing NA/NaN are dropped.
#'
#' @param active    Data coming from a well where a gene knockdown occurred.
#' @param control   Data belonging to a control well.
#' @param drop.feat A vector of strings of column names that will be dropped.
#' @param drop.sep  Drop variables that separate data. 
#' @param scale     Scale all matched features to the unit interval.
#' @param test      The fraction of rows to be used for testing is 1/test. If
#'                  NULL is supplied, all data is used for training.
#' @param verbose   Whether to list all modified/dropped features
#'
#' @return A list with entries "test" and "train" each holding a data frame
#'         containing a design matrix and a response vector.
#'
#' @examples
#' # get gene locations
#' mtor.loc <- findWells(experiments="brucella-du-k1", contents="MTOR")
#' scra.loc <- findWells(plates=sapply(mtor.loc, getBarcode),
#'                       contents="SCRAMBLED", well.names="G23")
#' # combine for faster fetching
#' data     <- getSingleCellData(list(mtor.loc[[1]], scra.loc[[1]]))
#' mtor.dat <- meltData(cleanData(data[[1]]$H6))
#' scra.dat <- meltData(cleanData(data[[1]]$G23))
#' # prepare data for glm
#' data <- prepareDataforGlm(mtor.dat$mat$Cells, scra.dat$mat$Cells)
#' data <- makeRankFull(data)
#' # run glm
#' model <- glm("Response ~ .", binomial, data$train)
#' 
#' @export
prepareDataforGlm <- function(active, control, drop.feat=NULL,
                              drop.sep=FALSE, scale="none", test=10,
                              verbose=FALSE) {

  if(!(is.data.frame(active) & is.data.frame(control))) {
    stop("expecting data frames as aguments active and control.")
  }

  active$Response <- factor(rep("active", nrow(active)),
                                levels=c("active", "control"))
  control$Response <- factor(rep("control", nrow(control)),
                                levels=c("active", "control"))

  data.all <- rbind(active, control)

  drop.feat <- c(drop.feat, "Image.Index", "Well.Index", "Well.Name", "Plate.Barcode")
  data.all <- data.all[, !names(data.all) %in% drop.feat]

  complete <- complete.cases(t(data.all))
  if(sum(!complete) > 0) {
    message("removing ", sum(!complete), " incomplete features")
    if(verbose) {
      l_ply(names(data.all[!complete]), function(x) message("  ", x))
    }
    data.all <- data.all[, complete]
    warning("removed ", sum(!complete), " variables containing Na/NaN.")
  }

  if(drop.sep) {
    drop       <- suppressMessages(analyzeSeparation(data.all, max.dim=1,
                                                     all.dim=FALSE))
    drop.names <- rownames(drop$single)[drop$single[,1] <= 0]
    if(length(drop.names) > 0) {
      message("dropping ", length(drop.names), " separating features")
      if(verbose) {
        l_ply(drop.names, function(x) message("  ", x))
      }
      data.all   <- data.all[, !names(data.all) %in% drop.names]
    }
  }
  if(scale != "none") {
    scale.ind <- grepl(scale, names(data.all), ignore.case=TRUE)
    if(length(scale.ind) > 0) {
      message("scaling ", length(scale.ind), " features to unit interval")
      if(verbose) {
        l_ply(names(data.all[scale.ind]), function(x) message("  ", x))
      }
      feat.scal <- apply(data.all[,scale.ind], 2, function(x) {
        return((x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)))
      })
      data.all <- cbind(data.all[,!scale.ind], data.frame(feat.scal))
      data.all <- data.all[,order(names(data.all))]
    }
  }

  if(!is.null(test)) {
    if(!is.numeric(test)) stop("expecting a numeric agrument for test.")
    set.seed(7)
    test.ind <- sort(sample.int(nrow(data.all), nrow(data.all) / test))
    testing <- data.all[test.ind,]
    trainin <- data.all[-test.ind,]
    res <- list(test=testing, train=trainin)
  } else {
    res <- list(test=NULL, train=data.all)
  }
}

#' Enforce full rank design matrix
#' 
#' Given a design matrix (if the response is also included it has to be named 
#' "Response"), rank deficiency is corrected by dropping zero variance
#' variables and one of each pair of highly correlated variables to make the
#' design matrix invertible.
#'
#' @param data Matrix/data frame holding the design matrix and optionally a
#'             response vector named "Response"
#'
#' @return A subset of the original data structure (some of the variables are
#'         removed)
#'
#' @examples
#' # get gene locations
#' mtor.loc <- findWells(experiments="brucella-du-k1", contents="MTOR")
#' scra.loc <- findWells(plates=sapply(mtor.loc, getBarcode),
#'                       contents="SCRAMBLED", well.names="G23")
#' # combine for faster fetching
#' data     <- getSingleCellData(list(mtor.loc[[1]], scra.loc[[1]]))
#' mtor.dat <- meltData(cleanData(data[[1]]$H6))
#' scra.dat <- meltData(cleanData(data[[1]]$G23))
#' # prepare data for glm
#' data <- prepareDataforGlm(mtor.dat$mat$Cells, scra.dat$mat$Cells)
#' data <- makeRankFull(data)
#' # run glm (enforcing full rank might be advisbale)
#' model <- glm("Response ~ .", binomial, data$train)
#' 
#' @export
makeRankFull <- function(data) {
  # input validation
  if(is.list(data)) {
    if(!is.data.frame(data)) {
      if(!is.null(data$train)) {
        temp <- data$train
      } else {
        stop("if data is a list, it is expected to contrain a slot \"train\".")
      }
    } else temp <- data
  } else stop("expecting a list or a data.frame for data.")
  if(!"Response" %in% colnames(temp)) {
    stop("expecting a column called \"Response\" in data")
  }

  temp <- temp[, !names(temp) %in% "Response"]
  remove <- NULL
  orig.rank <- qr(temp)$rank
  if(orig.rank < ncol(temp)) {
    # rank deficiency
    message("rank: ", orig.rank, " full: ", ncol(temp))
    col.var <- apply(temp, 2, var)
    remove <- names(which(col.var == 0))
    if(length(remove) > 0) {
      message("removing zero variance variables:\n  ",
              paste(remove, collapse="\n  "))
      temp <- temp[, !names(temp) %in% remove]
      warning("removed ", length(remove), " zero variance variables.")
    }
    new.rank <- qr(temp)$rank
    if(new.rank < ncol(temp)) {
      # still rank deficient
      message("rank: ", new.rank, " full: ", ncol(temp))
      one.cor <- which(cor(temp) > 0.9999, arr.ind = TRUE)
      off.diag <- one.cor[one.cor[,1] != one.cor[,2],]
      remo.ind <- unique(pmax(off.diag[,1], off.diag[,2]))
      if(length(remo.ind) > 0) {
        remo.nam <- names(temp)[remo.ind]
        message("removing one of highly correlated (>0.9999) variable pairs:",
                "\n  ", paste(remo.nam, collapse="\n  "))
        temp <- temp[, !names(temp) %in% remo.nam]
        warning("removed ", length(remo.ind), " variables due to highly ",
                "correlation (>0.9999) ")
        remove <- c(remove, remo.nam)
      }
    }
    new.rank <- qr(temp)$rank
    if(new.rank < ncol(temp)) {
      # still rank deficient
      stop("rank: ", new.rank, " full: ", ncol(temp))
    }
  }
  if(length(remove) > 0 & is.data.frame(data)) {
    data <- data[, !names(data) %in% remove]
  } else if (length(remove) > 0 & is.list(data)) {
    data$train <- data$train[, !names(data$train) %in% remove]
    data$test  <- data$test[, !names(data$test) %in% remove]
  }
  return(data)
}

#' Compare result of binary classification to truth
#'
#' Given two vectors of factors (two identical levels) of equal length,
#' calculate confusion matrix marginals
#'
#' @param estim Vector holding a two level factor with the classification
#'              result.
#' @param truth Vector holding a two level factor with the truth.
#'
#' @return A list with various key charactereistics of the resulting confusion
#'         matrix
#'
#' @examples
#' # get gene locations
#' mtor.loc <- findWells(experiments="brucella-du-k1", contents="MTOR")
#' scra.loc <- findWells(plates=sapply(mtor.loc, getBarcode),
#'                       contents="SCRAMBLED", well.names="G23")
#' # combine for faster fetching
#' data     <- getSingleCellData(list(mtor.loc[[1]], scra.loc[[1]]))
#' mtor.dat <- meltData(cleanData(data[[1]]$H6))
#' scra.dat <- meltData(cleanData(data[[1]]$G23))
#' # prepare data for glm
#' data <- prepareDataforGlm(mtor.dat$mat$Cells, scra.dat$mat$Cells)
#' data <- makeRankFull(data)
#' # run glm
#' model <- glm("Response ~ .", binomial, data$train)
#' # compare to testing data
#' predi <- as.factor(round(predict(model, newdata=data$test,
#'                                  type="response")))
#' levels(predi) <- c("active", "control")
#' comparison <- compareModeltoTruth(predi, data$test$Response)
#'
#' @export
compareModeltoTruth <- function(estim, truth) {
  if(length(estim) != length(truth)) {
    stop("comparing vectors of unequal length")
  }
  if(!identical(levels(truth), levels(estim))) {
    stop("levels don't match")
  }
  po <- levels(estim)[1]
  ne <- levels(estim)[2]
  tp <- sum(estim == po & truth == po)
  tn <- sum(estim == ne & truth == ne)
  fp <- sum(estim == po & truth == ne)
  fn <- sum(estim == ne & truth == po)
  tpr <- tp / (tp + fn)
  tnr <- tn / (fp + tn)
  fpr <- fp / (fp + tn)
  fnr <- fn / (fn + tp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  fdr <- fp / (fp + tp)
  fur <- fn / (fn + tn)
  f1s <- 2 * tp / (2 * tp + fp + fn)
  acc <- (tp + tn) / (tp + tn + fp + fn)
  n <- tn + tp + fn + fp
  s <- (tp + fn) / n
  p <- (tp + fp) / n
  mcc <- (tp / n - s * p) / sqrt(p * s * (1 - s) * (1 - p))
  return(list(POS=po, NEG=ne, TP=tp, TN=tn, FP=fp, FN=fn, TPR=tpr, TNR=tnr,
              FPR=fpr, FNR=fnr, PPV=ppv, NPV=npv, FDR=fdr, FOR=fur, F1S=f1s,
              ACC=acc, MCC=mcc)
  )
}

#' Find a hyperplane separating the sample points 
#'
#' Given a data frame holding feature data and a binary response column, find a
#' hyperplane separating the two groups of data points. The core code solving
#' the linear program is taken from the package safeBinaryRegression.
#'
#' @param data    A data.frame containing a column called Response, holding a
#'                binary vector.
#' @param max.dim The maximum number of features to simultaneously be
#'                considered for allowing data to be separated. Hence, if
#'                max.dim == 2, all pairs of features are screened. 
#' @param all.dim Analyze all features simulatneously.
#' @param tol1    The tolerance for considering the points separated: the sum
#'                of all betas has to be larger than tol1.
#' @param tol2    The tolerance for considering an individual feature to be
#'                involved.
#'
#' @return A list with two slots: a logical value saying if the data points are
#'         separated and a vector of distances to the separating hyperplane.
#'
#' @examples
#' # get gene locations
#' wells <- findWells(plates=c("J101-2C", "J107-2D", "J110-2L"),
#'                    well.names=c("H2", "H6"))
#' data  <- unlist(getSingleCellData(wells), recursive=FALSE)
#' h21 <- list(meta=data[["J101-2C.H2"]]$meta,
#'             data=meltData(cleanData(data[["J101-2C.H2"]], "lower")))
#' h61 <- list(meta=data[["J101-2C.H6"]]$meta,
#'             data=meltData(cleanData(data[["J101-2C.H6"]], "lower")))
#' dat1 <- prepareDataforGlm(h61$data$mat$Cells, h21$data$mat$Cells, test=NULL)
#' sep1 <- analyzeSeparation(dat1)
#'
#' @export
analyzeSeparation <- function(data, max.dim=2, all.dim=TRUE,
                              tol1=1e-3, tol2=1e-09) {
  separation <- function(x, y, tol1, purpose="find") {
    n <- dim(x)[1]
    p <- dim(x)[2]

    columnnames <- colnames(x)
    dimnames(x) <- NULL

    y.bar <- -sign(y - 0.5)
    x.bar <- y.bar * x

    lp <- make.lp(n, p)
    for(j in 1:p) status <- set.column(lp, j, x.bar[, j])
    status <- set.rhs(lp, rep(0.0, n))
    status <- set.constr.type(lp, rep(1, n))
    status <- set.objfn(lp, -colSums(x.bar))
    if(purpose == "test") {
      status <- set.bounds(lp, lower = rep(-Inf, p), upper = rep(Inf, p))
    } else if(purpose == "find") {
      status <- set.bounds(lp, lower=rep(-1, p), upper=rep(1, p))
    }
    control <- lp.control(lp, pivoting="firstindex", sense="max",
                          simplextype=c("primal", "primal"))
    status <- solve(lp)

    if(purpose == "test") {
      if(status == 0) return(0)
      else if(status == 3) return(1)
      else return(-1)
    } else if(purpose == "find") {
      if(status != 0) stop("unexpected result from lpSolveAPI for primal test")

      beta <- get.variables(lp)
      names(beta) <- columnnames

      if(sum(abs(beta)) > tol1) separation <- TRUE
      else separation <- FALSE
      return(list(separation=separation, beta=beta))
    }
  }

  if(is.list(data)) {
    if(!is.data.frame(data)) {
      if(!is.null(data$train)) {
        data <- data$train
      } else {
        stop("if data is a list, it is expected to contrain a slot \"train\".")
      }
    }
  } else stop("expecting a list or a data.frame for data.")
  if(!"Response" %in% colnames(data)) {
    stop("expecting a column called \"Response\" in data")
  }
  if(max.dim < 1) stop("max.dim expected to be >= 1")

  # first check each coordinate individually
  one.dim <- apply(data[, !names(data) %in% "Response"], 2, function(col, y) {
    acti <- col[y == "active"]
    ctrl <- col[y == "control"]
    min.a <- min(acti)
    max.a <- max(acti)
    min.c <- min(ctrl)
    max.c <- max(ctrl)
    sorted  <- sort(c(min.a, max.a, min.c, max.c))
    if((max.a < min.c & min.a < min.c & max.a > max.c) |
       (max.c < min.a & min.c < min.a & max.c > max.a)) {
      # negative overlap as completely separated intervals
      overlap <- sorted[2] - sorted[3]
      shared.act <- 0
      shared.ctr <- 0
    } else if((max.a == min.c & min.a <= min.c & max.a >= max.c) |
              (max.c == min.a & min.c <= min.a & max.c >= max.a)) {
      # zero overlap as quasi-separated
      overlap <- 0
      shared.act <- sum(acti == sorted[2]) / length(acti)
      shared.ctr <- sum(ctrl == sorted[2]) / length(ctrl)
    } else {
      # positive overlap as acutally overlapping
      overlap <- sorted[3] - sorted[2]
      shared.act <- sum(acti >= sorted[2] & acti <= sorted[3]) / length(acti)
      shared.ctr <- sum(ctrl >= sorted[2] & ctrl <= sorted[3]) / length(ctrl)
    }
    if((max(col) - min(col)) == 0) coverag <- 0
    else coverag <- overlap / (max(col) - min(col))
    return(c(coverag, shared.act, shared.ctr, min.a, max.a, min.c, max.c))
  }, data$Response)
  one.dim <- as.data.frame(t(one.dim))
  colnames(one.dim) <- c("coverage", "shared.act", "shared.ctr", "min.act",
                         "max.act", "min.ctrl", "max.ctrl")
  if(sum(one.dim$coverage < 0.3) > 0) {
    completely  <- one.dim[one.dim$coverage < 0,]
    quasi       <- one.dim[one.dim$coverage == 0,]
    interesting <- one.dim[as.logical((one.dim$coverage > 0) * 
                                      (one.dim$coverage < 0.3)),]
    col1  <- max(nchar(c(rownames(completely), rownames(quasi),
                         rownames(interesting))))
    if(nrow(completely) > 0) {
      message("for ", nrow(completely), " features, complete separation ",
              "occurs:\n  ", stri_pad_right("feature", col1), "  covarge ",
              "shrd.ac shrd.ct\n  ",
              paste(stri_pad_right(rownames(completely), col1),
                    formatC(completely$coverage, digits=4, width=6,
                            format="f"),
                    formatC(completely$shared.act, digits=4, width=6,
                            format="f"),
                    formatC(completely$shared.ctr, digits=4, width=6,
                            format="f"),
                    sep="  ", collapse="\n  "))
    }
    if(nrow(quasi) > 0) {
      message("for ", nrow(quasi), " features, quasi-separation ",
              "occurs:\n  ", stri_pad_right("feature", col1), "  covarge ",
              "shrd.ac shrd.ct\n  ",
              paste(stri_pad_right(rownames(quasi), col1),
                    formatC(quasi$coverage, digits=4, width=6,
                            format="f"),
                    formatC(quasi$shared.act, digits=4, width=6,
                            format="f"),
                    formatC(quasi$shared.ctr, digits=4, width=6,
                            format="f"),
                    sep="  ", collapse="\n  "))
    }
    if(nrow(interesting) > 0) {
      message("for ", nrow(interesting), " features, low coverage ",
              "occurs:\n  ", stri_pad_right("feature", col1), "  covarge ",
              "shrd.ac shrd.ct\n  ",
              paste(stri_pad_right(rownames(interesting), col1),
                    formatC(interesting$coverage, digits=4, width=6,
                            format="f"),
                    formatC(interesting$shared.act, digits=4, width=6,
                            format="f"),
                    formatC(interesting$shared.ctr, digits=4, width=6,
                            format="f"),
                    sep="  ", collapse="\n  "))
    }
  }

  if(max.dim >= 2 | all.dim) {
    # onl needed for further analysis
    if(sum(one.dim$coverage <= 0) > 0) {
      remove <- colnames(data)[one.dim$coverage <= 0]
      message("removing features from analysis:\n  ",
              paste(remove, collapse="\n  "))
      data <- data[, !colnames(data) %in% remove]
    }

    data <- makeRankFull(data)
    x <- as.matrix(data[, !colnames(data) %in% "Response"])
    y <- as.integer(data$Response) - 1
    if(length(unique(y)) != 2) stop("expecting binary response.")
  }

  if(max.dim >= 2) {
    dims <- 2:min(ncol(x), max.dim)
    n.dim <- lapply(dims, function(n.dim) {
      # every unordered n-tuple of indices is looked at
      if(choose(ncol(x), n.dim) > 2.5e5) {
        message("skipping ", n.dim, " dim step du to large number of ",
                "combinations (", choose(ncol(x), n.dim), ")")
        return(NULL)
      } else {
        indices <- t(combn(1:ncol(x), n.dim))
        res <- unlist(alply(indices, 1, function(ind) {
          return(separation(x[,ind], y, tol1, "test"))
        }, .progress=getOption("singleCellFeatures.progressBars")))
        crashed <- res < 0
        if(sum(crashed) > 0) {
          message("a total of ", sum(crashed), " errors were encountered")
          if(sum(crashed) < 10000) {
            a_ply(indices[crashed,], 1, function(row, names) {
              message("  ", names[row[1]], "\n    & ", names[row[2]])
            }, colnames(x))
          } else message("  not printing that much...")
        }
        separated <- res == 1
        if(sum(separated) > 0) {
          message("a total of ", sum(separated), " pairs with linear 2D ",
                  "(quasi-)separation were found")
          if(sum(separated) < 10000) {
            if(sum(separated) == 1) {
              row <- indices[separated,]
              message("  ", paste0(colnames(x)[row], collapse="\n    & "))
            } else {
              a_ply(indices[separated,], 1, function(row, names) {
                message("  ", paste0(colnames(x)[row], collapse="\n    & "))
              }, colnames(x))
              
              counts <- sort(table(indices[separated]), decreasing=TRUE)
              names  <- colnames(x)[as.integer(names(counts))]
              width  <- max(nchar(names))
              output <- paste(stri_pad_right(names, width),
                              format(counts), sep="  ")
              message("in summary, the feature x occurs n times:")
              l_ply(output, function(str) message("  ", str))
            }
          } else message("  not printing that much...")
        }
        res <- cbind(res, indices)
        colnames(res) <- c("separated", paste0("feat.", 1:n.dim))
        return(res)
      }
    })
    names(n.dim) <- paste0(dims, ".dim")
  } else {
    n.dim <- NULL
  }

  if(all.dim) {
    all.res <- separation(x, y, tol1)

    if(all.res$separation) {
      message("a separating hyperplane has been detected.")
      if(sum(abs(all.res$beta) > tol2) > 0) {
        beta.keep <- all.res$beta[abs(all.res$beta) < 1 - tol2]
        beta.sort <- beta.keep[order(abs(beta.keep))]
        sep.terms <- names(beta.sort)[abs(beta.sort) > tol2]
        sep.betas <- beta.sort[abs(beta.sort) > tol2]
        colwidth  <- max(nchar(sep.terms))
        output <- paste(stri_pad_right(sep.terms, colwidth),
                        format(sep.betas, digits=5), sep="  ")
        message("All ", length(all.res$beta), " features are involved (of ",
                ncol(data) - 1, ") and the following ", length(sep.terms),
                "\n  terms have abs(.) < ", 1 - tol2, ":")
        l_ply(output, function(string) message("  ", string))
      }
    }
  } else {
    all.res <- NULL
  }

  return(c(list(single=one.dim), n.dim, list(all=all.res)))
}
