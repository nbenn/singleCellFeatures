#' Analyze stability of top n coefficients
#' 
#' Resamples the data n.rep times, fits a glm modle using the function glm.fun
#' and records the top n.hit coefficient names. Returns a table with
#' frequencies of features that were among the top n.hit coefficients.
#'
#' @param well.a      A WellData object for glm analysis.
#' @param well.b      A second WellData object for glm analysis.
#' @param n.rep       The number of bootstrap repeats.
#' @param n.hit       Record this many top coefficients.
#' @param frac.sample Size of the subsample as fraction of all available data.
#' @param seed        Seed for reproducible sampling.
#' @param glm.fun     The glm function to use.
#' @param norm.feat   Which features to normalize. Either a regular
#'                    expression or "all".
#' @param norm.method How to normalize (see function normalizeData for
#'                    details). 
#' @param norm.sep    Separate data into wells form normalization? 
#' @param select.inf  Run stability analysis only on "infected", "uninfected"
#'                    or on "all"?
#' @param n.cores     The number of cores to use for bootstrapping.
#' @param ...         Further arguments to be passed to the glm function. For
#'                    example using glmnet, the alpha value.
#'
#' @return A vector of frequencies for features that were among the top n.hit
#'         coefficients.
#'
#' @examples
#' wells <- findWells(plates=c("J110-2C", "J104-2D", "J107-2L"),
#'                    well.names=c("H2", "H6"))
#' data  <- unlist(getSingleCellData(wells), recursive=FALSE)
#' a <- cleanData(data[["J110-2C.H2"]], "lower")
#' b <- cleanData(data[["J110-2C.H6"]], "lower")
#' test <- glmBootstrapStability(a, b, n.rep=20)
#' 
#' @export
glmBootstrapStability <- function(well.a, well.b, n.rep=100, n.hit=20,
                                  frac.sample=0.7, seed=7, glm.fun="glmnet",
                                  norm.feat="all", norm.method="zScore",
                                  norm.sep=FALSE, select.inf="all",
                                  n.cores=getNumCores(), ...) {

  if(any(class(well.a) == "WellData") & any(class(well.b) == "WellData")) {
    data.a <- meltData(well.a)
    data.b <- meltData(well.b)
    data <- prepareDataforGlm(data.a$mat$Cells, data.b$mat$Cells, test=NULL)
    if(norm.method != "none") {
      data$train <- normalizeData(data$train, norm.feat, norm.method, norm.sep)
    }
  } else {
    if(class(well.a) == "list" & class(well.b) == "list") {
      class.a <- sapply(well.a, function(x) any(class(x) == "WellData"))
      class.b <- sapply(well.b, function(x) any(class(x) == "WellData"))
      if(!all(class.a) & !all(class.b)) {
        stop("expecting WellData or lists of WellData objects for well.a/",
             "well.b")
      }
      data.a <- do.call(rbind, lapply(well.a, function(x) {
        meltData(x)$mat$Cells
      }))
      data.b <- do.call(rbind, lapply(well.b, function(x) {
        meltData(x)$mat$Cells
      }))
      data <- prepareDataforGlm(data.a, data.b, test=NULL)
      if(norm.method != "none") {
        data$train <- normalizeData(data$train, norm.feat, norm.method,
                                    norm.sep)
      }
    } else {
      stop("expexing lists or WellData objects for well.a/well.b")
    }
  }

  if(select.inf != "all") {
    infected <- as.logical(data$train$Cells.Infection_IsInfected)
    feature  <- !(names(data$train) %in% "Cells.Infection_IsInfected")
    if(select.inf == "infected") {
      message("dropping ", length(infected) - sum(infected), " non-infected ",
              "cells.")
      data$train <- data$train[infected,feature]
    } else if (select.inf == "uninfected") {
      message("dropping ", length(infected) - sum(!infected), " infected ",
              "cells.")
      data$train <- data$train[!infected,feature]
    } else stop("unrecognized string for param select.inf")
  }

  if(glm.fun == "glmnet") {
    dat.x <- as.matrix(data$train[,!(names(data$train) %in% "Response")])
    dat.y <- data$train$Response
  }

  registerDoParallel(cores=n.cores)
  if(glm.fun == "glmnet") {
    res <- foreach(i=1:n.rep, .combine=cbind) %dopar% {
      set.seed(i + seed)
      sample <- sample.int(nrow(dat.x), nrow(dat.x) * frac.sample)
      model  <- glmnet(dat.x[sample,], dat.y[sample], family="binomial", ...)
      last   <- model$beta[,ncol(model$beta)]
      last   <- last[order(abs(last), decreasing=TRUE)]
      names(last[1:n.hit])
    }
  } else if(glm.fun == "step") {
    res <- foreach(i=1:n.rep, .combine=cbind) %dopar% {
      set.seed(i + seed)
      sample <- sample.int(nrow(dat), nrow(dat) * frac.sample)
      subset <- dat[sample,]
      full  <- do.call("glm", list("Response ~ .", family=binomial,
                                   data=subset))
      empty <- do.call("glm", list("Response ~ 1", family=binomial,
                                   data=subset))
      model <- step(empty, scope=list(lower=formula(empty),
                    upper=formula(full)), direction="forward", trace=0,
                    steps=25)
      coeff <- model$coefficients
      coeff <- coeff[order(abs(coeff), decreasing=TRUE)]
      names(coeff[1:n.hit])
    }
  } else {
    res <- foreach(i=1:n.rep, .combine=cbind) %dopar% {
      set.seed(i + seed)
      sample <- sample.int(nrow(dat), nrow(dat) * frac.sample)
      subset <- dat[sample,]
      fun <- match.fun(glm.fun)
      model <- fun("Response ~ .", binomial, subset, ...)
      coeff <- model$coefficients
      coeff <- coeff[order(abs(coeff), decreasing=TRUE)]
      names(coeff[1:n.hit])
    }
  }

  counts <- as.data.frame(table(res))
  counts <- counts[order(counts$Freq, decreasing=TRUE),]
  rownames(counts) <- counts$res
  counts$res <- NULL

  output <- counts[1:min(20, nrow(counts)), , drop=FALSE]
  width  <- max(nchar(rownames(output)))
  output <- paste(stri_pad_right(rownames(output), width), output$Freq,
                  sep="  ")
  message("the feature x occurs n times:")
  l_ply(output, function(str) message("  ", str))

  return(counts)
}