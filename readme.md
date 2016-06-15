# singleCellFeatures
This R package serves to acquire and analyse single cell feature data as generated and compiled by TargetInfectX. The data is produced by [iBRAIN2](http://www.infectx.org/about_us/iBRAIN2/) based image analysis workflow and is hosted in an [openBIS](http://www.targetinfectx.ch/about_us/openBIS/) instance.

## Installation
The easiest way is throught the `install_github` functionality of the devtools package:

```R
install.packages("devtools")
library(devtools)

install_github("nbenn/singleCellFeatures")
```

Alternatively the package can be [downloaded](https://github.com/nbenn/singleCellFeatures/archive/master.zip) and installed manually:

```bash
unzip ~/Downloads/singleCellFeatures-master.zip
R CMD INSTALL --no-multiarch --with-keep.source ~/Downloads/singleCellFeatures-master
```

## Configuration
Some setup dependent information has to be provided, all of which is stored in a yaml file. The default location of this config file is `~/.singleCellFeaturesConfig`. This can be changed on a per-session basis using the function `configPathSet()` of more permanently, using an `.Rprofile` file.

```R
## if no config file is present, set one up
# set the config file location
configPathSet("path/to/where/you/want/your/config.yaml")
# create a template file
configInit()
# using a text editor, modify this file for your system

## for inter-session persistence, add the following to your .Rprofile
options(singleCellFeatures.configPath = "path/to/your/config.yaml")
```

The config file has the following structure:

```yaml
dataStorage:
  dataDir: "/path/to/data/dir"
  metaDir: "/path/to/metadata/dir"
beeDownloader:
  executable: "/path/to/repo/openBIS/Tools/BeeDataSetDownloader"
openBIS:
  username: "user"
  password: "password"
singleCellFeatures:
  sourceDir: "/path/to/source"
```

The two entries under `dataStorage` should be located on a volume with a couple of GB of free storage, to be able to hold a couple of plates. A complete plate uses 1-2 GB of storage, so having upwards of 50 GB available is recommended. The two entries under `beeDownloader` are concerned with `BeeDataSetDownloader`. More information available [here](https://wiki.systemsx.ch/pages/viewpage.action?title=InfectX+Single+Cell+Data+Access&spaceKey=InfectXRTD). The final section holds the path to the local source of this package. It is only used to update the databases in `/data` (see the `utilDatabase.R` file).
