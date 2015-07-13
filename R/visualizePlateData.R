#' Plot a plate heatmap
#' 
#' Given plate data, plot a heatmap of the specified feature using a function
#' to aggregate data within wells and an optional function to transform the
#' resulting vector of length 384.
#'
#' @param plate.dat     Data used to generate heatmap. Expecting a PlateData
#'                      object.
#' @param feature       Name of the feature to extract and aggreagte within
#'                      wells for heatmap plot. Expecting a string.
#' @param fun.aggregate The function used to aggregate data within wells.
#'                      Expecting either a function object or a string.
#' @param fun.transform An optional function used to transform the data before
#'                      plotting. Expecting either a function object or a
#'                      string.
#'
#' @return A ggplot2 plot object ready for printing.
#'
#' @examples
#' plate <- PlateLocation("J101-2C")
#' data  <- PlateData(plate)
#' plateHeatmap(data, "VoronoiCells.Location_Center_Y", median, sqrt)
#'
#' @export

plateHeatmap <- function(plate.dat, feature, fun.aggregate="mean",
                         fun.transform="identity") {
  if(!any(class(plate.dat) == "PlateData")) {
    stop("expecting PlateData for parameter \"plate.dat\".")
  }
  plate.vals <- sapply(plate.dat$data, function(well, feat, fun) {
    well.vals <- sapply(well$data, function(img, feat) {
      res.vec <- lapply(img$data.vec, function(grp, feat) {
        if(feat %in% colnames(grp)) {
          return(grp[[feat]])
        }
      }, feat)
      res.mat <- lapply(img$data.mat, function(grp, feat) {
        if(feat %in% colnames(grp)) {
          return(grp[,feat])
        }
      }, feat)
      res <- unlist(list(res.vec, res.mat), recursive=FALSE)
      res <- res[!sapply(res, is.null)]
      if(length(res) != 1) {
        stop("expecting to find exactly 1 feature, instead found ",
             length(res), ".")
      }
      return(unlist(res))
    }, feat)
    res <- unlist(well.vals)
    func <- match.fun(fun)
    return(func(res))
  }, feature, fun.aggregate)

  if(length(plate.vals) != 384) {
    stop("expecting 384 wells, instead got ", length(plate.vals), ".")
  }

  well.type <- sapply(plate.dat$data, function(well) {
    return(well$meta$well.type)
  })

  res <- cbind(plate.vals, rep(LETTERS[1:16], each=24), rep(1:24, 16))
  colnames(res) <- c("value", "rows", "cols")
  res <- data.frame(res)
  transf <- match.fun(fun.transform)
  res$value <- as.numeric(levels(res$value))[res$value]
  res$value <- transf(res$value)
  res$rows <- factor(res$rows, levels=rev(LETTERS[1:16]), ordered=TRUE)
  res$cols <- factor(res$cols, levels=1:24, ordered=TRUE)

  myPalette <- colorRampPalette(brewer.pal(9, "Reds")[3:9], space="Lab")
  if(as.character(substitute(fun.transform)) == "identity") {
    legend.title <- as.character(substitute(fun.aggregate))
  } else {
    legend.title <- paste(as.character(substitute(fun.aggregate)),
                          as.character(substitute(fun.transform)), sep="\n")
  }

  frames <- data.frame(cbind(well.type, rep(1:24, 16),
                             rev(rep(1:16, each=24))), stringsAsFactors=FALSE)
  names(frames) <- c("well.type", "cols", "rows")
  frames$cols <- as.integer(frames$cols)
  frames$rows <- as.integer(frames$rows)
  frames$color <- ifelse(frames$well.type == "SIRNA", "white", "black")

  heat <- ggplot(data=res) +
    geom_raster(aes(x=cols, y=rows, fill=value)) +
    geom_rect(data=frames, size=0.5, fill=NA, colour=frames$color,
              aes(xmin=cols - 0.475, xmax=cols + 0.475,
                  ymin=rows - 0.475, ymax=rows + 0.475)) +
    scale_x_discrete(name="") +
    scale_y_discrete(name="") +
    theme_bw() +
    theme(axis.ticks=element_blank()) +
    ggtitle(feature) +
    scale_fill_gradientn(colours = myPalette(100),
                         name=legend.title) +
    coord_fixed()
  return(heat)
}