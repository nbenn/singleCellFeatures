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

## Loading data
# import the plate with all available features
path  <- paste("/Users/nbennett/Polybox/MasterThesis/openBISDownload/",
               "INFECTX_PUBLISHED/SALMONELLA_TEAM/SALMONELLA-AU-K1/KB2-04-1A", 
               sep="")
plate <- importCompletePlate(path, NULL)

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
