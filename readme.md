# singleCellFeatures
This R package serves to acquire and analyse single cell feature data as generated and compiled by TargetInfectX. The data is produced by [iBRAIN2](http://www.infectx.org/about_us/iBRAIN2/) based image analysis workflow and is hosted in an [openBIS](http://www.targetinfectx.ch/about_us/openBIS/) instance.

## Installation
The easiest way is with help of the devtools package:

```R
install.packages('devtools')
library(devtools)

install_github('nbenn/singleCellFeatures')
```

Alternatively it can be [downloaded](https://github.com/nbenn/singleCellFeatures/archive/master.zip) and installed manually:

```bash
unzip ~/Downloads/singleCellFeatures-master.zip
R CMD INSTALL --no-multiarch --with-keep.source ~/Downloads/singleCellFeatures-master
```

## Configuration
