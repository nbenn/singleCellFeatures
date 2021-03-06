#' Find PlateLocations
#' 
#' Find all plates using the following set of parameters:
#'
#' @param pathogens       Narrow down the results by selecting pathogens
#'                        (exact, matches, ignoring case).
#' @param experiments     Name of the experiment to be considered (regular
#'                        expression, case insensitive).
#' @param plates          A vector of regular expressions for plates.
#' @param well.types      A vector of regular expressions that is matched with
#'                        the well.types column in well databases.
#' @param contents        This vector of strings/integers is matched against
#'                        the three well database columns id.openBIS,
#'                        id.infx and name via case insensitive regular
#'                        expressions. In case it is specified, the next three
#'                        parameters are ignored.
#' @param id.openBIS      A vector of regular expressions that is matched
#'                        against the id.openBIS column.
#' @param id.infx         A vector of regular expressions that is matched
#'                        against the id.infx column.
#' @param name            A vector of regular expressions that is matched
#'                        against the name column.
#'
#' @return A list of PlateLocation objects.
#'
#' @examples
#' set1 <- findPlates(contents="MTOR", experiments="brucella-du-k")
#' set2 <- findPlates(contents=2475, experiments="brucella-au-k[1-3]")
#' set3 <- findPlates(contents="SCRAMBLED", experiments="brucella-du-k")
#'
#' @export

findPlates <- function(pathogens=NULL, experiments=NULL, plates=NULL,
                       well.types=NULL, contents=NULL, id.openBIS=NULL,
                       id.infx=NULL, name=NULL, verbose=FALSE) {

  data(plateDatabase, envir=environment())
  curr.plates <- plate.database
  if(verbose) message("starting with ", nrow(curr.plates), " plates.")
  # if pathogens specified, exclude others
  curr.plates <- processPathogens(pathogens, curr.plates, verbose)
  # if experiments specified, exclude others
  curr.plates <- processExperiments(experiments, curr.plates, verbose)
  # if plates specified, exclude others
  curr.plates <- processPlates(plates, curr.plates, verbose)
  # print current plate set
  printPlateSet(verbose, curr.plates)
  if(any(c(!is.null(well.types), !is.null(contents), !is.null(id.openBIS),
           !is.null(id.infx), !is.null(name)))) {
    # find all well databases to be loaded
    well.db <- getWells(curr.plates, verbose)
    # if well types specified, exclude others
    well.db <- processWellTypes(well.types, well.db, verbose)
    # if contents specified, exclude others
    well.db <- processContents(contents, id.openBIS, id.infx, name, well.db,
                               verbose)
    # remove wells on plates that have no single cell features available
    keep <- well.db$barcode %in% plate.database$Barcode
    if(verbose & sum(!keep) > 0) {
      message("removing ", sum(!keep), " wells because they have no single",
              " cell features available.")
    }
    well.db <- well.db[keep,]
    # turn wells back into plates
    result <- unique(well.db$barcode)
  } else {
    result <- unique(curr.plates$Barcode)
  }

  if(length(result) == 0) {
    stop("no matching plates found.")
  } else {
    message("there are ", length(result), " plates remaining:\n  ",
            paste(result, collapse=", "))
    res.lst <- lapply(result, function(bc) return(PlateLocation(bc)))
    return(res.lst)
  }
}

#' Find WellLocations
#' 
#' Find all wells on all plates using the following set of parameters:
#'
#' @param pathogens       Narrow down the results by selecting pathogens
#'                        (exact, matches, ignoring case).
#' @param experiments     Name of the experiment to be considered (regular
#'                        expression, case insensitive).
#' @param plates          A vector of regular expressions for plates.
#' @param well.rows       All rows to be included (a vector of characters,
#'                        exact match, case insensitive)
#' @param well.cols       The columns to be included (a vector of integers)
#' @param well.names      A vector of well names (e.g. B14), expecting strings,
#'                        case insensitive matching.
#' @param well.types      A vector of regular expressions that is matched with
#'                        the well.types column in well databases.
#' @param contents        This vector of strings/integers is matched against
#'                        the three well database columns id.openBIS,
#'                        id.infx and name via case insensitive regular
#'                        expressions. In case it is specified, the next three
#'                        parameters are ignored.
#' @param id.openBIS      A vector of regular expressions that is matched
#'                        against the id.openBIS column.
#' @param id.infx         A vector of regular expressions that is matched
#'                        against the id.infx column.
#' @param name            A vector of regular expressions that is matched
#'                        against the name column.
#'
#' @return A list of WellLocation objects.
#'
#' @examples
#' set1 <- findWells(contents="MTOR", experiments="brucella-du-k")
#' set2 <- findWells(contents=2475, experiments="brucella-au-k[1-3]")
#' set3 <- findWells(contents="SCRAMBLED", experiments="brucella-du-k")
#'
#' @export

findWells <- function(pathogens=NULL, experiments=NULL, plates=NULL,
                      well.rows=NULL, well.cols=NULL, well.names=NULL,
                      well.types=NULL, contents=NULL, id.openBIS=NULL,
                      id.infx=NULL, name=NULL, verbose=FALSE) {

  data(plateDatabase, envir=environment())
  curr.plates <- plate.database
  if(verbose) message("starting with ", nrow(curr.plates), " plates.")
  # if pathogens specified, exclude others
  curr.plates <- processPathogens(pathogens, curr.plates, verbose)
  # if experiments specified, exclude others
  curr.plates <- processExperiments(experiments, curr.plates, verbose)
  # if plates specified, exclude others
  curr.plates <- processPlates(plates, curr.plates, verbose)
  # print current plate set
  printPlateSet(verbose, curr.plates)
  # find all well databases to be loaded
  well.db <- getWells(curr.plates, verbose)
  # if well rows specified, exclude others
  well.db <- processWellRows(well.rows, well.db, verbose)
  # if well columns specified, exclude others
  well.db <- processWellCols(well.cols, well.db, verbose)
  # if well names specified, exclude others
  well.db <- processWellNames(well.names, well.db, verbose)
  # if well types specified, exclude others
  well.db <- processWellTypes(well.types, well.db, verbose)
  # if contents specified, exclude others
  well.db <- processContents(contents, id.openBIS, id.infx, name, well.db,
                             verbose)
  # remove wells on plates that have no single cell features available
  keep <- well.db$barcode %in% plate.database$Barcode
  if(verbose & sum(!keep) > 0) {
    message("removing ", sum(!keep), " wells because they have no single cell",
            " features available.")
  }
  well.db <- well.db[keep,]

  message("there are ", nrow(well.db), " wells remaining:")
  out <- cbind(
    stri_pad_right(well.db$barcode, max(nchar(well.db$barcode))),
    stri_pad_right(paste0(well.db$well.row, well.db$well.col), 3),
    stri_pad_right(well.db$well.type, max(nchar(well.db$well.type))),
    stri_pad_right(well.db$id.openBIS, max(nchar(well.db$id.openBIS))),
    stri_pad_right(well.db$id.infx,
                   max(nchar(well.db$id.infx))),
    stri_pad_right(well.db$name, max(nchar(well.db$name)))
  )
  apply(out, 1, function(row) {
    message("  ", paste(row, collapse="  "))
  })

  if(nrow(well.db) == 0) stop("no matching wells found.")
  res.lst <- apply(well.db, 1, function(row) {
    WellLocation(row[["barcode"]], row[["well.row"]], row[["well.col"]])
  })
  return(res.lst)
}

processPathogens <- function(pathogens, curr.plates, verbose) {
  if(!is.null(pathogens)) {
    if(!is.character(pathogens)) {
      stop("expecting a vector of characters for pathogens")
    }
    group.names <- sapply(pathogens, function(pathogen) {
      if(tolower(pathogen) == "mock") return(toupper(pathogen))
      else return(paste0(toupper(pathogen), "_TEAM"))
    })
    curr.plates <- curr.plates[curr.plates$Group %in% group.names,]
    if(verbose) {
      message("after applying pathogens, ", nrow(curr.plates), " plates ",
              "remaining.")
    }
  }
  return(curr.plates)
}

processExperiments <- function(experiments, curr.plates, verbose) {
  if(!is.null(experiments)) {
    if(!is.character(experiments)) {
      stop("expecting a vector of characters for experiments")
    }
    found <- unlist(sapply(experiments, function(exp) {
          return(grep(exp, curr.plates$Experiment, ignore.case=TRUE))
    }))
    curr.plates <- curr.plates[unique(found),]
    if(verbose) {
      message("after applying experiments, ", nrow(curr.plates), " plates ",
              "remaining.")
    }
  }
  return(curr.plates)
}

processPlates <- function(plates, curr.plates, verbose) {
  if(!is.null(plates)) {
    if(!is.character(plates)) {
      stop("expecting a vector of characters for plates")
    }
    found <- unlist(sapply(plates, function(plate) {
          return(grep(plate, curr.plates$Barcode, ignore.case=TRUE))
    }))
    curr.plates <- curr.plates[unique(found),]
    if(verbose) {
      message("after applying plates, ", nrow(curr.plates), " plates ",
              "remaining.")
    }
  }
  return(curr.plates)
}

printPlateSet <- function(verbose, curr.plates) {
  if(verbose) {
    message("there are ", nrow(curr.plates), " plates remaining:")
    out <- cbind(
      stri_pad_right(curr.plates$Barcode, max(nchar(curr.plates$Barcode))),
      stri_pad_right(curr.plates$Space, max(nchar(curr.plates$Space))),
      stri_pad_right(curr.plates$Group, max(nchar(curr.plates$Group))),
      stri_pad_right(curr.plates$Experiment,
                     max(nchar(curr.plates$Experiment))),
      stri_pad_right(curr.plates$DataID, max(nchar(curr.plates$DataID)))
    )
    apply(out, 1, function(row) {
      message("  ", paste(row, collapse="  "))
    })
  }
}

getWells <- function(curr.plates, verbose) {
  well.db <- lapply(unique(curr.plates$Group), function(pathogen) {
    if(pathogen == "MOCK") {
      db.name  <- "wellDatabaseMock"
      obj.name <- "well.database.mock"
    }
    else {
      db.name <- unlist(strsplit(pathogen, "_"))
      if(length(db.name) != 2 | db.name[2] != "TEAM") {
        stop("unexpected group name")
      }
      db.lower <- tolower(db.name[1])
      db.name  <- paste0(toupper(substring(db.lower, 1, 1)),
                        substring(db.lower, 2))
      db.name  <- paste0("wellDatabase", db.name)
      obj.name <- paste0("well.database.", db.lower)
    }
    data(list=db.name, envir=environment())
    return(get(obj.name))
  })
  # combine into large well db
  well.db <- do.call(rbind, well.db)
  # select only the plates found above
  well.db <- well.db[well.db$barcode %in% curr.plates$Barcode,]
  if(verbose) {
    message("using the current set of plates, ", nrow(well.db),
            " wells remain.")
  }
  return(well.db)
}

processWellRows <- function(well.rows, well.db, verbose) {
  if(!is.null(well.rows)) {
    if(!is.character(well.rows)) {
      stop("expecting a vector of characters for well.rows")
    }
    well.db <- well.db[well.db$well.row %in% well.rows,]
    if(verbose) {
      message("after applying well.rows, ", nrow(well.db), " wells ",
              "remaining.")
    }
  }
  return(well.db)
}

processWellCols <- function(well.cols, well.db, verbose) {
  if(!is.null(well.cols)) {
    if(!is.integer(well.cols)) {
      stop("expecting a vector of integers for well.cols")
    }
    well.db <- well.db[well.db$well.col %in% well.cols,]
    if(verbose) {
      message("after applying well.cols, ", nrow(well.db), " wells ",
              "remaining.")
    }
  }
  return(well.db)
}

processWellNames <- function(well.names, well.db, verbose) {
  if(!is.null(well.names)) {
    if(!is.character(well.names)) {
      stop("expecting a vector of characters for well.names")
    }
    all.names <- paste0(well.db$well.row, well.db$well.col)
    well.db <- well.db[all.names %in% well.names,]
    if(verbose) {
      message("after applying well.names, ", nrow(well.db), " wells ",
              "remaining.")
    }
  }
  return(well.db)
}

processWellTypes <- function(well.types, well.db, verbose) {
  if(!is.null(well.types)) {
    if(!is.character(well.types)) {
      stop("expecting a vector of characters for well.types")
    }
    found <- unlist(sapply(well.types, function(type) {
      return(grep(type, well.db$well.type, ignore.case=TRUE))
    }))
    well.db <- well.db[unique(found),]
    if(verbose) {
      message("after applying well.types, ", nrow(well.db), " wells ",
              "remaining.")
    }
  }
  return(well.db)
}

processContents <- function(contents, id.openBIS, id.infx, name, well.db,
                            verbose) {
  if(!is.null(contents)) {
    if(!(is.character(contents) | is.integer(contents))) {
      stop("expecting a vector of characters or integers for contents")
    }
    if(!(is.null(id.openBIS) | is.null(id.infx) | is.null(name))) {
      warning("when using contents, id.openBIS, id.infx and name are ",
              "ignored.")
    }
    found <- unlist(sapply(contents, function(content) {
      found1 <- grep(content, well.db$id.openBIS, ignore.case=TRUE)
      found2 <- grep(content, well.db$id.infx, ignore.case=TRUE)
      found3 <- grep(content, well.db$name, ignore.case=TRUE)
      return(unique(c(found1, found2, found3)))
    }))
    well.db <- well.db[unique(found),]
    if(verbose) {
      message("after applying contents, ", nrow(well.db), " wells ",
              "remaining.")
    }
  } else {
    # if id.openBIS specified, exclude others
    if(!is.null(id.openBIS)) {
      if(!(is.character(id.openBIS) | is.integer(id.openBIS))) {
        stop("expecting a vector of characters or integers for id.openBIS")
      }
      found <- unlist(sapply(id.openBIS, function(x) {
        return(grep(x, well.db$id.openBIS, ignore.case=TRUE))
      }))
      well.db <- well.db[unique(found),]
      if(verbose) {
        message("after applying id.openBIS, ", nrow(well.db), " wells ",
                "remaining.")
      }
    }
    # if id.infx specified, exclude others
    if(!is.null(id.infx)) {
      if(!(is.character(id.infx) | is.integer(id.infx))) {
        stop("expecting a vector of characters or integers for id.infx")
      }
      found <- unlist(sapply(id.infx, function(x) {
        return(grep(x, well.db$id.infx, ignore.case=TRUE))
      }))
      well.db <- well.db[unique(found),]
      if(verbose) {
        message("after applying id.infx, ", nrow(well.db), " wells ",
                "remaining.")
      }
    }
    # if name specified, exclude others
    if(!is.null(name)) {
      if(!(is.character(name) | is.integer(name))) {
        stop("expecting a vector of characters or integers for name")
      }
      found <- unlist(sapply(name, function(x) {
        return(grep(x, well.db$name, ignore.case=TRUE))
      }))
      well.db <- well.db[unique(found),]
      if(verbose) {
        message("after applying name, ", nrow(well.db), " wells ",
                "remaining.")
      }
    }
  }
  return(well.db)
}
