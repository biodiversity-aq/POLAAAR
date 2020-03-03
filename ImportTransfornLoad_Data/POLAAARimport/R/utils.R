#==============================================================
# POLAAARimport package
#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2020-03-03)
# file encdong UTF-8
#
#==============================================================
# This script contains the different functions that support the working of metadata.to.polaaar(), which is the main import function to write data to POLAAAR.
#
#==============================================================

#' extract information from an onine EML document
#'
#' @author Maxime Sweetlove ccBY 4.0 2019
#'
#' @description get.polaaar.EML.data extracts the information needed for the POLAAAR database from a given onine EML document
#' @param EML.url the URL to the EML document
#' @return a list with the data needed for the POLAAAR database. This includes: abstract, start_date, end_date, bounding_box, is_public and associated_references
#' @export
#'
get.polaaar.EML.data <- function(EML.url){
  # requires EML, emld

  eml_file <- emld::as_emld(EML.url, from = "xml")

  dataset_name <- unname(unlist(eml_get(eml_file,"title", from="xml"))["title"])

  # abstract
  if(any(grepl("abstract", as.character(eml_file)))){
    abstract <- eml_get(eml_file,"abstract", from="xml")$para
    if(is.null(abstract)){
      abstract <- ""
    }
  }else{
    abstract <- ""
  }

  # start_date | end_date
  if(any(grepl("formationPeriod", as.character(eml_file)))){
    dates <- eml_get(eml_file,"formationPeriod", from="xml")
    dates <- unlist(dates)[1]
    start_date  <- strsplit(dates, " ")[[1]][1]
    end_date  <- strsplit(dates, " ")[[1]][2]
  } else if(any(grepl("temporalCoverage", as.character(eml_file)))){
    dates <- eml_get(eml_file,"temporalCoverage", from="xml")
    if("rangeOfDates" %in% names(dates)){
      start_date  <- gsub('\n| |\t', '', dates$rangeOfDates$beginDate)
      end_date  <- gsub('\n| |\t', '', dates$rangeOfDates$endDate)
    } else if("singleDateTime" %in% names(dates)){
      start_date  <- gsub('\n| |\t', '', dates$singleDateTime)
      end_date <- start_date
    }else{
      start_date  <- ""
      end_date <- ""
    }
  }

  # bounding_box
  if(any(grepl("boundingCoordinates", as.character(eml_file)))){
    bounding_box <- unlist(EML::eml_get(eml_file, "boundingCoordinates", from = "xml"))
    bounding_box <- paste("SRID=4326;POLYGON ((",
                          bounding_box["northBoundingCoordinate"], ", ",
                          bounding_box["eastBoundingCoordinate"], ", ",
                          bounding_box["southBoundingCoordinate"], ", ",
                          bounding_box["westBoundingCoordinate"], "))", sep="")
  }else{
    bounding_box<-""
  }

  # is_public
  if(any(grepl("pubDate", as.character(eml_file)))){
    is_public <- TRUE
  }else{
    is_public <- FALSE
  }

  # associated_references
  if(any(grepl("bibliography", as.character(eml_file)))){
    associated_references <- unlist(EML::eml_get(eml_file, "bibliography", from = "xml"))["citation"]
    associated_references <- unname(associated_references)
  }else{
    associated_references <- ""
  }

  output<-list(name = dataset_name,
               abstract = abstract,
               start_date = start_date,
               end_date = end_date,
               bounding_box = bounding_box,
               is_public = is_public,
               associated_references = associated_references
  )

  return(output)
}

#' preprocess metadata.MIxS for POLAAAR
#'
#' @author Maxime Sweetlove ccBY 4.0 2019
#'
#' @description preprocess metadata.MIxS for POLAAAR
#' @param metadata.object a metadata.MIxS object
#' @return a metadata.MIxS object that has been prepared for metadata.to.polaaar()
#'
preprocess.polaaar.MIxS <- function(metadata.object){
  ## note: this is a preprocessing function, no QC is executed
  #' @param metadata.object a metadata.MIxS object
  #' @return a metadata.MIxS object that has been prepared for metadata.to.polaaar()

  metaunits <- metadata.object@units
  metasection <- metadata.object@section
  dataset <- metadata.object@data

  # a. samplingProtocol
  if("samp_collect_device" %in% colnames(dataset)){
    dataset$samplingProtocol <- dataset$samp_collect_device
  }

  newTerms <- setdiff(colnames(dataset), names(metaunits))
  for(trm in newTerms){
    metaunits[trm] <- as.character(TermsLib[TermsLib$name=="samplingProtocol",]$expected_unit)
    metasection[trm] <- as.character(TermsLib[TermsLib$name=="samplingProtocol",]$MIxS_section)
  }

  New_metadata <- new("metadata.MIxS",
                      data   = dataset,
                      section = metasection,
                      units      = metaunits,
                      env_package = metadata.object@env_package,
                      type = metadata.object@type,
                      QC = metadata.object@QC
  )

  return(New_metadata)
}

#' preprocess DarwinCore data for POLAAAR
#'
#' @author Maxime Sweetlove ccBY 4.0 2019
#'
#' @description preprocess DarwinCore data for POLAAAR
#' @param dataset a data.frame formated with DarwinCore terms
#' @return a data.frame that has been prepared for metadata.to.polaaar()
#'
preprocess.polaaar.DwC.core <- function(dataset){
  ## note: this is a preprocessing function, no QC is executed

  # a. create a collection_date and time field
  if("eventDate" %in% colnames(dataset)){
    dataset$collection_date <- dataset$eventDate
  }
  if("eventTime" %in% colnames(dataset)){
    dataset$collection_time <- dataset$eventTime
  }

  # b. geo_loc_name
  geoTerms <- intersect(colnames(dataset), c("continent", "country", "state_province",
                                             "waterBody", "islandGroup", "island", "locality"))
  if(length(geoTerms)!=0){
    for(gt in geoTerms){
      dataset[,gt]<-paste(gt, "=", dataset[,gt], sep="")
      dataset[,gt]<-gsub(".+=$", "",dataset[,gt], fixed=FALSE)
    }
    dataset$geo_loc_name <- apply(dataset[,geoTerms], 1, function(x) paste(x[!is.na(x) & x != ""], collapse = ";", sep=""));
    dataset$geo_loc_name<-gsub(";;", ";",dataset$geo_loc_name, fixed=TRUE)
  }

  # c. dynamicProperties to individualcolumns
  # doing this cell by cell, because the content can be anything...
  # expected format: {"colNameWithUnits"=value}, {...}
  if("dynamicProperties" %in% colnames(dataset)){
    for(i in 1:nrow(dataset)){
      dc<-dataset[i,]$dynamicProperties
      if(grepl("\\{", dc)){
        # format brackets to one type
        dc<-gsub("\\[", "\\{", dc)
        dc<-gsub("\\]", "\\}", dc)
        dc<-gsub("\\(", "\\{", dc)
        dc<-gsub("\\)", "\\}", dc)

        #remove quotes
        dc<-gsub("\"", "", dc)
        dc<-gsub("\'", "", dc)

        #try to split fieds
        dc<-gsub("\\} \\{", "---", dc)
        dc<-gsub("\\}, \\{", "---", dc)
        dc<-gsub("\\},\\{", "---", dc)
        dc<-gsub("\\}\\{", "---", dc)
        dc<-gsub("\\{", "", dc)
        dc<-gsub("\\}", "", dc)
        dc<-strsplit(dc, "---")[[1]]

        for(dx in dc){
          name<-strsplit(dx, "=")[[1]][1]
          name<-gsub(" ", "", name)
          value<-strsplit(dx, "=")[[1]][2]
          value<-gsub(" ", "", value)
          if(!name %in% colnames(dataset)){
            dataset[,name]<-""
            dataset[i,name]<-value
          }else{
            dataset[i,name]<-value
          }
        }

      }

    }

  }

  # get associatedSequences

  ## associatedSequences to emof??

  return(dataset)
}
