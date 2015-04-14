## Download data
# cd /Users/nbennett/ETH/MasterThesis/software/tginfectx/trunk/openBIS/Tools/
#   BeeDataSetDownloader
# export BEESOFTSRC="/Users/nbennett/ETH/MasterThesis/software/tginfectx/trunk"
# ./BeeDataSetDownloader.sh 
#   --user "nbennett@student.ethz.ch"
#   --password "****"
#   --outputdir "/Users/nbennett/Polybox/MasterThesis/openBISDownload"
#   --plateid "^/INFECTX_PUBLISHED/SALMONELLA_TEAM/SALMONELLA-AU-K1/KB2-04-1A"
#   --files ".*.mat"
#   --verbose "10"

## Load functions
library(singleCellFeatures)

## Define functions
makeDataFrame <- function(dataList) {
  if(all(sapply(dataList, function(x, first) {
    ifelse(length(x) == length(first), TRUE, FALSE)
  }, dataList[[1]]))) {
    df <- do.call(cbind, lapply(dataList, data.frame, stringsAsFactors=FALSE))
    names(df) <- sapply(names(dataList), function(x) {
      unlist(strsplit(x, "[.]"))[2]
    })
    return(df)
  } else {
    print("not all same length")
    return(NULL)
  }
}

## Setup data location
team       <- "SALMONELLA_TEAM"
experiment <- "SALMONELLA-AU-K1"
plate      <- "KB2-04-1A"
path       <- paste("/Users/nbennett/Polybox/MasterThesis/openBISDownload",
              "INFECTX_PUBLISHED", team, experiment, plate, sep="/")
setwd(path)

## Loading data
# import the plate with all available features
bacteria <- importCompletePlate(path, "^Bacteria\\.")
cells    <- importCompletePlate(path, "^Cells\\.")
nuclei   <- importCompletePlate(path, "^Nuclei\\.")
perinuc  <- importCompletePlate(path, "^PeriNuclei\\.")
voronoi  <- importCompletePlate(path, "^VoronoiCells\\.")
  
bacteria.i.8 <- extractPartialPlate(bacteria, cols=8, rows="I", 
                                    imgs=5)$data$I8$data$img_22
cells.i.8    <- extractPartialPlate(cells, cols=8, rows="I", 
                                    imgs=5)$data$I8$data$img_22
nuclei.i.8   <- extractPartialPlate(nuclei, cols=8, rows="I", 
                                    imgs=5)$data$I8$data$img_22
perinuc.i.8  <- extractPartialPlate(perinuc, cols=8, rows="I", 
                                    imgs=5)$data$I8$data$img_22
voronoi.i.8  <- extractPartialPlate(voronoi, cols=8, rows="I", 
                                    imgs=5)$data$I8$data$img_22

## Location plots
# using ggplot
library(ggplot2)
library(grid)
library(png)

cells_all.png  <- readPNG("I8/2-2_cells.png")
cells_all.grob <- rasterGrob(img, interpolate=TRUE)

cells.i.8.df <- makeDataFrame(cells.i.8)

ggplot(cells.i.8.df, aes(x=Location_Center_X, y=-Location_Center_Y, 
                         size=AreaShape_Area, label=rownames(cells.i.8.df)), 
       guide=FALSE) +
  annotation_custom(cells_all.grob, xmin=0, xmax=1410, ymin=-1040, ymax=0) +
  geom_point(colour="white", fill="red", shape=21) +
  scale_size_area(max_size = 15) +
  geom_text(size=4) +
  theme_bw()

## Scatterplot matrix
library(car)
bacteria.i.8.df <- makeDataFrame(bacteria.i.8[1:19])
scatterplotMatrix(bacteria.i.8.df, diagonal="histogram", smooth=FALSE)


ggplot(cells.i.8.df) +
  geom_star(aes(x=Location_Center_X, y=-Location_Center_Y, r=AreaShape_Area, angle = date,
                fill = mean(temperature)), r.zero = FALSE)

hist(i.8$data$I8$data$img_22$Cells.Intensity_StdIntensity_CorrActin)

beanplot(i.8$data$I8$data$img_22$Cells.AreaShape_Area,
         i.8$data$I8$data$img_22$Cells.AreaShape_Perimeter,
         i.8$data$I8$data$img_22$Cells.Intensity_IntegratedIntensity_CorrActin,
         i.8$data$I8$data$img_22$Cells.Intensity_IntegratedIntensity_CorrDNA,
         i.8$data$I8$data$img_22$Cells.Intensity_IntegratedIntensity_CorrPathogen)

plot(i.8$data$I8$data$img_22$Cells.Location_Center_X,
     i.8$data$I8$data$img_22$Cells.Location_Center_Y)

plot(i.8$data$I8$data$img_22$Bacteria.Location_Center_X,
     i.8$data$I8$data$img_22$Bacteria.Location_Center_Y)


sum(sapply(dat.plat$data$B10$data, function(x) length(x[[1]])))

subset.img <- extractPartialPlate(dat.plat)
subset.wel <- extractPartialPlate(dat.plat, cols=c(3,4), rows=c("A","B"), 
                                  wels=c(333, 21), imgs=NULL)

# import whole experiment
path     <- paste("/Users/nbennett/Polybox/MasterThesis/openBISDownload/",
                  "INFECTX_PUBLISHED/ADENO_TEAM", sep="")
plates   <- c("KB2-03-1I", "KB2-03-1J")
features <- c("Cells.AreaShape_Area", "Cells.AreaShape_Eccentricity")
dat.expe <- importMultiplePlates(path, NULL, features)

