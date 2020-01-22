###library() the polartools package for the QC functions
MarsLib<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/MarsLibrary.csv", 
                  header=TRUE)

polaaarLib<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/PolaaarLibrary.csv", 
                                 header=TRUE)



write.metadata.MIxS.as.mars <- function(metadata.object, name.prefix=NULL, dest.dir=getwd(),
                                        add.missing.data.columns=TRUE, ask.input=TRUE){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description write.metadata.MIxS.as.mars takes a metadata.MIxS object and outputs it as a MiMARKS and SequenceSet table to upload to mARS
  #' @param metadata.object a metadata.MIxS class object. The object to be submitted to mARS
  #' @param name a character string naming a file. Either a full path or a file name to write to the working directory. Invalid input will cause output to be printed to the console.
  #' @param checklist_accession character. The name of a MIxS environmetal package or it's ENA checklist accession number.
  #' @details ENA uses its own variant of MIxS, and has
  #' write.metadata.MIxS.as.ENA assumes the dataset has already been subjected to process the process.metadata function to standardize in input. If this is not the case, "garbage in, garbage out" is applicable. 
  #' @return tab separated .txt file
  #' @example 
  #' 
  
  warningmessages<-c()
  
  # 1. check input
  if(check.valid.metadata.MIxS(metadata.object)){
    metaunits <- metadata.object@units
    metasection <- metadata.object@section
    metapackage <- metadata.object@env_package
    metadata <- metadata.object@data
  }else{
    stop("metadata.object argument must be a metadata.MIxS class object.
         Common dataframes should be converted to metadata.MIxS with 
         the data.frame.to.metadata.MIxS function.")
  } 
  # 2. check output destination
  if(!is.character(name.prefix) | c(NULL,NA) %in% name.prefix | length(name.prefix)>1){
    stop("Need input for the name.prefix argument.")
  } 
  if(!dir.exists(file.path(dest.dir))){
    stop("Could not find the directory provided in the dest.dir argument.")
  }
  mimarks.name <- paste(dest.dir, "/MiMARKS_", name.prefix, ".csv", sep="")
  seqset.name <- paste(dest.dir, "/SeqSet_", name.prefix, ".csv", sep="")
  
  # 3. create MiMARKS tabble
  # 3.1 first row
  Row01 <- paste("section", "Structured Comment Name", "units", 
                 paste(rownames(metadata), collapse=","), sep=",")
  # 3.2 general terms
  mimarksTerms <- as.character(MarsLib[MarsLib$mars_mimarks==TRUE,]$name) #the basic terms that need to be in every dataset
  #associated units
  mimarksUnits <- as.character(TermsLib[match(setdiff(mimarksTerms, names(metaunits)), TermsLib$name),]$expected_unit)
  if(length(mimarksUnits)!=0){
    names(mimarksUnits) <- setdiff(mimarksTerms, names(metaunits))
  }
  mimarksUnits <- c(mimarksUnits, metaunits)
  #assocated MIxS sections
  mimarksSection <- as.character(TermsLib[match(setdiff(mimarksTerms, names(metasection)), TermsLib$name),]$section)
  if(length(mimarksSection)!=0){
    names(mimarksSection) <- setdiff(mimarksTerms, names(metasection))
  }
  mimarksSection <- c(mimarksSection, metasection)
  Row02<-character()
  #some additional steps for the primer field (to split forward and reverse primer in the seqset file)
  investigation_type<-character() 
  fw_primerName <- character() 
  fw_primerSeq <- character() 
  rv_primerName <- character()
  rv_primerSeq <- character()
  mimarks_data <- metadata # save data for the seqset round
  for(tx in 1:length(mimarksTerms)){
    Row0x_section <- as.character(mimarksSection[names(mimarksSection)==mimarksTerms[tx]])
    Row0x_section <- paste("MiMARKS_", Row0x_section,sep="")
    Row0x_name <- mimarksTerms[tx]
    Row0x_units <- as.character(mimarksUnits[names(mimarksUnits)==mimarksTerms[tx]])
    if(mimarksTerms[tx] %in% colnames(mimarks_data)){
      Row0x_data <- mimarks_data[,colnames(mimarks_data)==mimarksTerms[tx]]
      #remove the used data to see what columns will remain at the end
      mimarksUnits <- mimarksUnits[!names(mimarksUnits)==mimarksTerms[tx]] 
      mimarksSection <- mimarksSection[!names(mimarksSection)==mimarksTerms[tx]] 
      mimarks_data <- mimarks_data[,!colnames(mimarks_data)==mimarksTerms[tx],drop=FALSE] 
      if(mimarksTerms[tx]=="primer"){
        Row0x_data
        fw_primerSeq <- character() 
        rv_primerSeq <- character()
      }
    } else if(ask.input){
      if(mimarksTerms[tx] == "primer" && "mimarks-survey" %in% investigation_type){ ###special case
        cat(paste("no primer information found:\n\tPlease provide the forward primer name, or hit enter to leave blank.\n", sep=""))
        fw_primerName <- readline() 
        cat(paste("no primer information found:\n\tPlease provide the forward primer sequence, or hit enter to leave blank.\n", sep=""))
        fw_primerSeq <- readline() 
        cat(paste("no primer information found:\n\tPlease provide the reverse primer name, or hit enter to leave blank.\n", sep=""))
        rv_primerName <- readline() 
        cat(paste("no primer information found:\n\tPlease provide the reverse primer sequence, or hit enter to leave blank.\n", sep=""))
        rv_primerSeq <- readline() 
        Row0x_data <- rep(paste("(", fw_primerName, "):",fw_primerSeq, " (", rv_primerName, "):",rv_primerSeq, sep=""), nrow(mimarks_data))
      } else{
        cat(paste("no data found for ",mimarksTerms[tx],":\n\tPlease type the info to fill the cells, or hit enter to leave blank.\n", sep=""))
        user_input <- readline() 
        Row0x_data <- rep(user_input, nrow(mimarks_data))
      }
    }else{
      Row0x_data <- rep("", nrow(mimarks_data))
    }
    Row0x_data[is.na(Row0x_data)]<-""
    Row0x_data <- gsub("^NA$", "", Row0x_data)
    if(Row0x_name =="investigation_type"){
      investigation_type<-unique(Row0x_data)
    }
    Row0x_data <- paste(Row0x_data, collapse=",")
    Row0x <- paste(Row0x_section, Row0x_name, Row0x_units, Row0x_data, sep=",")
    Row02 <-  paste(Row02, Row0x, sep="\n")
    
  }
  
  # 3.3 remaining dataset-specific terms
  doubleTerms <- c("decimalLatitude", "decimalLongitude")
  mimarksUnits <- mimarksUnits[!names(mimarksUnits) %in% doubleTerms] 
  mimarksSection <- mimarksSection[!names(mimarksSection) %in% doubleTerms] 
  mimarks_data <- mimarks_data[,!colnames(mimarks_data) %in% doubleTerms, drop=FALSE]
  if(ncol(mimarks_data)>0){
    Row03 <- character()
    for(cx in 1:ncol(mimarks_data)){
      #assuming the data has gone through the QC of process.metadata()
      Row0x_section <- mimarksSection[colnames(mimarks_data[cx])]
      if(!Row0x_section %in% c("package", "miscellaneous")){
        Row0x_section <- paste("MiMARKS_", Row0x_section,sep="")
      }
      Row0x_name <- colnames(mimarks_data[cx])
      Row0x_units <- mimarksUnits[colnames(mimarks_data[cx])]
      Row0x_data <- mimarks_data[,cx]
      Row0x_data[is.na(Row0x_data)]<-""
      Row0x_data <- gsub("^NA$", "", Row0x_data)
      Row0x_data <- paste(Row0x_data, collapse=",")
      Row0x <- paste(Row0x_section, Row0x_name, Row0x_units, Row0x_data, sep=",")
      Row03 <- paste(Row03, Row0x, sep="\n")
    }
  }
  mimarks.data <- paste(Row01, Row02, Row03, sep="")
  
  # 4. create SeqSet tabble
  # 4.1 first row
  Row01 <- paste("unique_sequence_set_id", 
                 paste(rownames(metadata), collapse=","), sep=",")
  # 4.2 other rows
  seqsetTerms <- as.character(MarsLib[MarsLib$mars_seqset==TRUE,]$name)
  Row02<-character()
  for(tx in 1:length(seqsetTerms)){
    Row0x_name <- seqsetTerms[tx]
    if(seqsetTerms[tx] %in% colnames(metadata)){
      Row0x_data <- metadata[,colnames(metadata)==seqsetTerms[tx]]
    } else{
      mimarksEquiv <- as.character(MarsLib[MarsLib$name==seqsetTerms[tx],]$mimarks_equivalent)[1]
      if(mimarksEquiv == "investigation_type"){ ###special case
        Row0x_data <- gsub("mimarks-survey", "marker gene", Row0x_data)
      }
      if(nchar(mimarksEquiv)==0){ #case it had no mimarks equivalent or wasn't in the colnames of metadata
        if(ask.input){
          cat(paste("no data found for ",seqsetTerms[tx],":\n\tPlease type the info to fill the cells, or hit enter to leave blank.\n", sep=""))
          user_input <- readline() 
          if(user_input == ""){
            Row0x_data <- rep("", nrow(metadata))
          } else{
            Row0x_data <- rep(user_input, nrow(metadata))
          }
        }else{
          Row0x_data <- rep("", nrow(metadata))
        }
      }else{
        Row0x_data <- metadata[,colnames(metadata)==mimarksEquiv]
      }
    }
    Row0x_data[is.na(Row0x_data)]<-""
    Row0x_data <- gsub("^NA$", "", Row0x_data)
    Row0x_data <- paste(Row0x_data, collapse=",")
    Row0x <- paste(Row0x_name, Row0x_data, sep=",")
    Row02 <-  paste(Row02, Row0x, sep="\n")
  }
  seqset.data <- paste(Row01, Row02, sep="")
  
  # 5. finalize
  cat(paste("The data had been written to ", dest.dir, "\n",sep=""))
  
  write.table(mimarks.data, file=mimarks.name, 
              col.names = FALSE, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  write.table(seqset.data, file=seqset.name, 
              col.names = FALSE, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  }



###database acces info

#DATABASES = {
#  'default': {
#    'ENGINE': 'django.contrib.gis.db.backends.postgis',
#    'NAME': 'd86hngmf7e303f',
#    'USER': 'hwyoarfegnhrwn',
#    'PASSWORD': '286e142ad8f2d8aa4d5ad529c43bbe847f4037ddc4f065f51bf18366cfac4cd5',
#    'HOST': 'ec2-54-217-234-157.eu-west-1.compute.amazonaws.com',
#    'PORT': '5432'
#  }
#}

library(RPostgreSQL)
con <- dbConnect(PostgreSQL(), host="ec2-54-217-234-157.eu-west-1.compute.amazonaws.com", 
                 user= "hwyoarfegnhrwn", 
                 password="286e142ad8f2d8aa4d5ad529c43bbe847f4037ddc4f065f51bf18366cfac4cd5", 
                 dbname="d86hngmf7e303f", port=5432)

dbListTables(conn)

table_name<- dbReadTable(conn, "tablen_name")


dbGetQuery(con,
           "SELECT table_name FROM information_schema.tables
           WHERE table_schema='sch2014'")


dbDisconnect(con)



#EML.url <- "https://api.gbif.org/v1/dataset/ebe73e19-11eb-48d2-9263-fb0caa3b7b5a/document"
#EML.url <- "https://ipt.biodiversity.aq/eml.do?r=airborne_bacteria_miers_antarctica&v=1.2"

get.polaaar.EML.data <- function(EML.url){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description get.polaaar.EML.data extracts the information needed for the POLAAAR database from a given onine EML document
  #' @param EML.url the URL to the EML document
  #' @return a list of the data needed for the POLAAAR database, including: abstract, start_date, end_date, bounding_box, is_public and associated_references

  require(EML)
  require(emld)

  eml_file <- emld::as_emld(EML.url, from = "xml")
  
  dataset_name <- unname(unlist(eml_get(eml_file,"title", from="xml"))["title"])
  
  # abstract
  if(any(grepl("abstract", as.character(eml_file)))){
    abstract <- eml_get(eml_file,"abstract", from="xml")$para
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




metadata.to.polaaar <- function(metadata.object, EML.url=NA, user_ID="", 
                                dest.dir=getwd(), create.dir=FALSE){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description metadata.to.polaaar writes a metadata.MIxS or metadata.DwC object to the polaaar database
  #' @param metadata.object a metadata.MIxS class object or metadata.DwC. The object to be formatted to fit the polaaar database schema
  #' @param EML.url a named vector. URL to the EML (IPT) associated with the dataset, written as a named vector, e.g. c(project1="http://url1", project2="http://url1"). The names in the vector must correspond to the different projects listed under "project_name" in the metadata object. If only one link is provided and is unnamed, than this will automatically be the parent project to all the samples in the metadata object.
  #' @param user_ID ID of the user that send the data
  #' @param dest.dir a character vector. The complete path to the directory where the CSV files should be written to. Default is the working diectory. If create.dir is TRUE, the given directory will be created if it does not exist.
  #' @param create.dir boolean. If TRUE, the directory path in the variable dest.dir will be created if it does not already exist. Default FALSE
  #' @details 
  #' @return data is written to polaaar

  # important! This function assumes data has been QC'd by process.metadata()
  
  require(RCurl)

  # 1. check input
  #****************************************************************************
  if(check.valid.metadata.MIxS(metadata.object)){
    metaformat <- "MIxS"
    metaunits <- metadata.object@units
    metasection <- metadata.object@section
    metapackage <- metadata.object@env_package
    metadata <- metadata.object@data
  }else if(check.valid.metadata.DwC(metadata.object)){
    ### experimantal part, still to be develloped
    metaformat <- "DwC"
    metapackage <- "not_specified"
    coredata <- metadata.object@core
    metadata <- metadata.object@emof
  }else{
    stop("metadata.object argument must be a metadata.MIxS or metadata.DwC class object.\n\tCommon dataframes should be converted to metadata.MIxS or metadata.DwC with \n\tthe data.frame.to.metadata.MIxS function.")
  } 
  
  # the data must be associated with an EML reccord, otherwise it should not be added to the database
  if(is.na(EML.url)){
    stop("No EML.url provided.\n\tAny dataset must be associated with a project, of which the project metadata has been\n\tsubmitted to the IPT.")
  }
  # multiple EML URLs must be named
  if(length(EML.url)>1 & !is.null(names(EML.url))){
    stop("If multiple EML URLs provided, they must be named with names corresponding to\n\tproject names in the \"project_name\" column of the metadata object")
  }

  # check if all EML.urls lead to an existing web page
  for(emlurl in EML.url){
    if(!RCurl::url.exists(emlurl)){
      StopMessage <- paste("The following EML page does not seem to exist:\n\t", emlurl, sep="")
      stop(StopMessage)
    }
  }
  
  # to keep track of columns in the input data that have been sent to the polaaarDB
  ColsDone <- c()
  num_recs <- nrow(metadata)
 #IDtable <- do we need a table to keep track of IDs?

  # check output destination
  if(!dir.exists(file.path(dest.dir))){
    if(create.dir){
      dir.create(dest.dir)
    }else{
      stop("the output directory does not exist.\n\tCheck for typo's.")
    }
  }
  
  
  # 2. filling the Polaaar database (PDB) tables
  #****************************************************************************
  #------------------------------------------
  #### create Event table ####
  #------------------------------------------
  # Event $ parent_event [ref: ParentEvent $ id]
  # Event $ metadata [ref: Metadata $ id]
  # Event $ occurence [ref: Occurence $ id]
  # Event $ environment [ref: Environment $ id]
  
  EventTerms <- c("id", "footprintWKT", "eventRemarks", "sample_name", "collection_date", 
                  "collection_time", "parent_event", "parent_sample", 
                  "samplingProtocol", "occurrence", "metadata", "environment", "metadata_exists",
                  "occurrence_exists", "environment_exists")
  Event <- data.frame(matrix(nrow=num_recs, ncol=length(EventTerms), data=""))
  colnames(Event) <- EventTerms
  rownames(Event) <- rownames(metadata)
  
  # Event $ id (should be generated by the PDB)
  Event$id <- 1:num_recs
  
  # Event $ footprintWKT = coordinates in well-known text format
  if("footprintWKT" %in% colnames(metadata)){
    Event$footprintWKT<-metadata$footprintWKT
    ColsDone <- c(ColsDone, "footprintWKT")
  }else if("decimalLatitude" %in% colnames(metadata) &
           "decimalLongitude" %in% colnames(metadata)){
    Event$footprintWKT <- paste("SRID=4326; POINT (", 
                                     as.character(metadata$decimalLatitude),
                                     ", ", as.character(metadata$decimalLongitude), ")", 
                                     sep="")
  }else if("lat_lon" %in% colnames(metadata)){
    Event$footprintWKT <- paste("SRID=4326; POINT (", 
                                     gsub(" ", ", ", as.character(metadata$lat_lon)), ")", 
                                     sep="")
  }
  
  Event$footprintWKT <- gsub("POINT (, )", NA, Event$footprintWKT)
  Event$footprintWKT <- gsub("POINT (NA, NA)", NA, Event$footprintWKT)
  
  # Event $ eventRemarks / collection_date / collection_time / samplingProtocol
  for(tm in c("eventRemarks", "collection_date", "collection_time", "samplingProtocol")){
    if("tm" %in% colnames(metadata)){
      Event$eventRemarks<-metadata$eventRemarks
      ColsDone <- c(ColsDone, "eventRemarks")
    }
  }
  
  # Event $ sample_name
  if("original_name" %in% colnames(metadata)){
    Event$sample_name<-metadata$original_name
    ColsDone <- c(ColsDone, "original_name")
  } else if("INSDC_SampleID" %in% colnames(metadata)){
    Event$sample_name<-metadata$INSDC_SampleID
    ColsDone <- c(ColsDone, "INSDC_SampleID")
  } else if("occurenceID" %in% colnames(metadata)){
    Event$sample_name<-metadata$occurenceID
    ColsDone <- c(ColsDone, "occurenceID")
  } else{
    Event$sample_name <- colnames(metadata) 
  }
  
  #------------------------------------------
  #### create ParentEvent table ####
  #------------------------------------------
  # starts with same number of rows as dataset, but will later be reduced to unique eventIDs
  # if there are no eventIDs, the project(s) will assume the role of parent event to all the samples
  # ParentEvent $ event_type [ref: EventType $ id]
  # ParentEvent $ project [ref: ProjectMetadata $ id]
  # ParentEvent $ id [ref: Event $ parent_event]
  # hierarchy for event_type: sample [=row in metadata] > event [groups samples] > parentEvent [groups events] > project [groups samples, events or parentEvents]
  ParentEventTerms <- c("id", "parent_event_name", "event_type", "description", "parent_event",
                        "event_creator", "created_on", "updated_on", "project")
  ParentEvent <- data.frame(matrix(nrow=num_recs, ncol=length(ParentEventTerms), data=""))
  colnames(ParentEvent) <- ParentEventTerms
  
  #add eventID
  if("eventID" %in% colnames(metadata)){
    if(length(unique(metadata$eventID))!=num_recs){
      ParentEvent$parent_event_name <- metadata$eventID
      ParentEvent$event_type <- rep("event", num_recs)
      # link to Event table (first name, change to id later)
      Event$parent_event <- ParentEvent$parent_event_name
      # link to the project (first name, change to id later)
      PRJ <- "case1"
      if("project_name" %in% colnames(metadata)){
        ParentEvent$project <- metadata$project_name
      }else if(!is.null(names(EML.url)) & length(EML.url)==1){
        ParentEvent$project <- rep(names(EML.url), nrow(ParentEvent))
      } else if("bioproject" %in% colnames(metadata)){
        ParentEvent$project <- metadata$bioproject
      }else{
        ParentEvent$project <- rep(1, num_recs)
      }
      
      # look for higher level events (i.e. parent events)
      if("parentEventID" %in% colnames(metadata)){
        # ParentEvent $ parent_event
        ParentEvent$parent_event <- metadata$parentEventID
        # add new rows for the parent events
        PE_dict <- unique(metadata[,c("eventID", "parentEventID")])
        num_PE <- length(unique(PE_dict$parentEventID))
        PE.temp <- data.frame(matrix(nrow=num_PE, ncol=length(ParentEventTerms), data=""))
        colnames(PE.temp) <- ParentEventTerms
        PE.temp$parent_event_name <- unique(PE_dict$parentEventID)
        PE.temp$event_type <- rep("parentEvent", num_PE)
        PE.temp <- PE.temp[!PE.temp$parent_event_name %in% c("", "NA", NA),] #remove NA and ""
      
        #project details of the parentEvents
        PRJ_dict <- unique(ParentEvent[,c("parent_event", "project")])
        PE.temp$project <- unlist(sapply(PE.temp$parent_event_name, 
                                         FUN = function(x){
                                           gsub(x,PRJ_dict[PRJ_dict$parent_event==x,]$project,x)
                                           }))
        PE.temp$parent_event <- PE.temp$project

        ParentEvent <- rbind(ParentEvent, PE.temp)
        
      }else{
        #project details of the parentEvents
        PRJ_dict <- unique(ParentEvent[,c("parent_event", "project")])
        ParentEvent$parent_event <- ParentEvent$project
      }
      
      #add project row
      num_PRJ <- length(unique(PRJ_dict$project))
      PRJ.temp <- data.frame(matrix(nrow=num_PRJ, ncol=length(ParentEventTerms), data=""))
      colnames(PRJ.temp) <- ParentEventTerms
      PRJ.temp$parent_event_name <- unique(PRJ_dict$project)
      PRJ.temp$project <- PRJ.temp$parent_event_name
      PRJ.temp$event_type <- rep("project", num_PRJ)
      PRJ.temp <- PRJ.temp[!PRJ.temp$parent_event_name %in% c("", "NA", NA),] #remove NA and ""
      
      ParentEvent <- rbind(ParentEvent, PRJ.temp)
      
      #remove duplicated events
      ParentEvent <- unique(ParentEvent)
      # ParentEvent $ id
      ParentEvent$id <- 1:nrow(ParentEvent)
      # Event $ parent_event [ref: ParentEvent $ id]
      # replace name with ID in Event table
      Event$parent_event <- unlist(sapply(Event$parent_event, FUN = function(x){gsub(x,ParentEvent[ParentEvent$parent_event_name==x,]$id,x)}))
      
      # ParentEvent $ parent_event [ref: ParentEvent $ id]
      ParentEvent$parent_event <- unlist(sapply(ParentEvent$parent_event, 
                                                FUN = function(x){
                                                  if(x!=""){
                                                    gsub(x,ParentEvent[ParentEvent$parent_event_name==x,]$id,x)
                                                    }else{x<-""}
                                                  }))
      
    } else if(length(unique(metadata$eventID))==num_recs &&
              "parentEventID" %in% colnames(metadata)){
      #each sample is a unique event, so the events are meaningless, but there are parentEvents
      #the parentEvent will take the place of the event
      # ParentEvent $ parent_event
      ParentEvent$parent_event_name <- metadata$parentEventID
      ParentEvent$event_type <- rep("event", num_recs)
      # link to Event table (first name, change to id later)
      Event$parent_event <- ParentEvent$parent_event_name
      # link to the project (first name, change to id later)
      PRJ <- "case1"
      if("project_name" %in% colnames(metadata)){
        ParentEvent$project <- metadata$project_name
      }else if(!is.null(names(EML.url)) & length(EML.url)==1){
        ParentEvent$project <- rep(names(EML.url), nrow(ParentEvent))
      }else if("bioproject" %in% colnames(metadata)){
        ParentEvent$project <- metadata$bioproject
      }else{
        ParentEvent$project <- rep(1, num_recs)
      }
      ParentEvent$parent_event <- ParentEvent$project
      
      PRJ_dict <- unique(ParentEvent[,c("parent_event", "project")])
      #add project row
      num_PRJ <- length(unique(PRJ_dict$project))
      PRJ.temp <- data.frame(matrix(nrow=num_PRJ, ncol=length(ParentEventTerms), data=""))
      colnames(PRJ.temp) <- ParentEventTerms
      PRJ.temp$parent_event_name <- unique(PRJ_dict$project)
      PRJ.temp$project <- PRJ.temp$parent_event_name
      PRJ.temp$event_type <- rep("project", num_PRJ)
      PRJ.temp <- PRJ.temp[!PRJ.temp$parent_event_name %in% c("", "NA", NA),] #remove NA and ""
      
      ParentEvent <- rbind(ParentEvent, PRJ.temp)
      
      #remove duplicated events
      ParentEvent <- unique(ParentEvent)
      # ParentEvent $ id
      ParentEvent$id <- 1:nrow(ParentEvent)
      # Event $ parent_event [ref: ParentEvent $ id]
      # replace name with ID in Event table
      Event$parent_event <- unlist(sapply(Event$parent_event, FUN = function(x){gsub(x,ParentEvent[ParentEvent$parent_event_name==x,]$id,x)}))

      ParentEvent$parent_event <- unlist(sapply(ParentEvent$parent_event, 
                                                FUN = function(x){
                                                  if(x!=""){
                                                    gsub(x,ParentEvent[ParentEvent$parent_event_name==x,]$id,x)
                                                  }else{x<-""}
                                                }))
      
    } else if(length(unique(metadata$eventID))==num_recs &&
             !"parentEventID" %in% colnames(metadata)){
      #each sample is a unique event, so the events are meaningless
      PRJ <- "case2"
          }
  } else if(!"eventID" %in% colnames(metadata)){
    PRJ <- "case2"
  }
  
  ### add rows for the projects
  if(PRJ=="case1"){ #->project is the highest level event
    
  }else{ #PRJ=="case2" ->project is the only higher level event
    # link to the project (first name, change to id later)
    if("project_name" %in% colnames(metadata)){
      ParentEvent$parent_event_name <- metadata$project_name
    }else if(!is.null(names(EML.url)) & length(EML.url)==1){
      ParentEvent$project <- rep(names(EML.url), nrow(ParentEvent))
    }else if("bioproject" %in% colnames(metadata)){
      ParentEvent$parent_event_name <- metadata$bioproject
    }else{
      ParentEvent$parent_event_name <- rep("unnamed_project", num_recs)
    }
    Event$parent_event <- ParentEvent$parent_event_name
    
    ParentEvent$project <- ParentEvent$parent_event_name
    ParentEvent$event_type <- rep("project", num_recs)
    #remove duplicated events
    ParentEvent <- unique(ParentEvent)
    # ParentEvent $ id
    ParentEvent$id <- 1:nrow(ParentEvent)
    # Event $ parent_event [ref: ParentEvent $ id]
    # replace name with ID in Event table
    Event$parent_event <- unlist(sapply(Event$parent_event, FUN = function(x){gsub(x,ParentEvent[ParentEvent$parent_event_name==x,]$id,x)}))
  }
  
  # ParentEvent $ description
  # no standard term for event description yet...
  
  # ParentEvent $ created_on / updated_on
  ParentEvent$created_on <- rep(Sys.Date(), nrow(ParentEvent))
  ParentEvent$updated_on <- ParentEvent$created_on
  
  # ParentEvent $ event_creator
  # take user ID of the user that send the data, or main author of the paper if the data was harvested that way
  ParentEvent$event_creator <- rep(user_ID, nrow(ParentEvent))

  #------------------------------------------
  #### create ProjectMetadata table ####
  #------------------------------------------
  # preferably 1 dataset == 1 project, if so it can be perfectly linked to the EML
  # project_qaqc is a boolean to indicate if a dataset has been QCd.
  # ProjectMetadata $ id [ref: ParentEvent $ project]
  # ProjectMetadata $ associatedReferences [ref: Reference $ id]
  # project_creator = user id of person that emailed the spreadsheet
  # boundingbox notation:
  # SRID=4326;POLYGON ((N, E, S, W))
  ProjectMetadataTerms <- c("id", "project_name", "start_date", "end_date", "EML_URL", "abstract", 
                            "bounding_box", "is_public", "associated_references", 
                            "associated_media", "created_on", "updated_on", 
                            "project_creator", "project_qaqc")
  ProjectMetadata <- data.frame(matrix(nrow=length(EML.url), ncol=length(ProjectMetadataTerms), data=""), stringsAsFactors = FALSE)
  colnames(ProjectMetadata) <- ProjectMetadataTerms
  
  # some checks first
  if(length(unique(ParentEvent$project))>0 &
     !is.null(names(EML.url)) &&
     !all(names(EML.url) %in% unique(ParentEvent$project))){  
    #EML URL names must correspond to project names (but project names can have no corresponding EML name)
    #otherwise there will be records created in the ProjectMetadata table that correspond to nothing
    stop("There are EML.url names that do not correspond to project names in the metadata object.")
  }
  
  # ProjectMetadata $ id
  ProjectMetadata$id <- 1:length(EML.url) 
  
  # ProjectMetadata $ EML_URL
  ProjectMetadata$EML_URL <- EML.url
  
  # ProjectMetadata $ project_creator
  ProjectMetadata$project_creator <- rep(user_ID, nrow(ProjectMetadata))
  
  # ProjectMetadata $ project_qaqc
  # assume this to be true if it is being added to the database
  ProjectMetadata$project_qaqc <- rep(TRUE, nrow(ProjectMetadata))
  
  # ProjectMetadata $ project_name
  if(!is.null(names(EML.url))){
    # the EML.url has names
    # all EML names correspond to project names (has been checked)
    # EML names get priority, project names that do not correspond to an EML name will be droped (only in the ParentEvent$project column)
    ProjectMetadata$project_name <- names(EML.url) 
    # revisit ParentEvent $ project and remove all project names with no link to an EML
    ParentEvent$project <- unlist(sapply(ParentEvent$project, 
                                         FUN = function(x){
                                           if(x %in% ProjectMetadata$project_name){
                                             gsub(x,ProjectMetadata[ProjectMetadata$project_name==x,]$id,x)
                                           }else{
                                             x<-""
                                           }
                                           }))
  }else{
    # no EML names, can only happen if there was just one URL provided (was checked in the beginning)
    # in that case all samples will be linked to this EML, regardless what project name has been provided
    ProjectMetadata$project_name <- "unnamed_project"
    # revisit ParentEvent $ project and remove all project names with no link to an EML
    ParentEvent$project <- rep(1, nrow(ParentEvent))
  }
  
  
  #### getting data from EML webpage
   for(emlurl in EML.url){
     eml_row <- as.numeric(row.names(ProjectMetadata[match(emlurl,ProjectMetadata$EML_URL),]))
     eml_data <- get.polaaar.EML.data(emlurl)
     #ProjectMetadata $ abstract
     ProjectMetadata[eml_row,"abstract"] <- eml_data$abstract
     #ProjectMetadata $ start_date
     ProjectMetadata[eml_row,"start_date"] <- eml_data$start_date
     #ProjectMetadata $ end_date
     ProjectMetadata[eml_row,"end_date"] <- eml_data$end_date
     #ProjectMetadata $ bounding_box
     ProjectMetadata[eml_row,"bounding_box"] <- eml_data$bounding_box
     #ProjectMetadata $ is_public
     ProjectMetadata[eml_row,"is_public"] <- eml_data$is_public
     #ProjectMetadata $ associated_references
     ProjectMetadata[eml_row,"associated_references"] <- eml_data$associated_references
     
     #if no name for the project was provided, change it to the dataset name found in the EML
     if(ProjectMetadata[eml_row,"project_name"] == "unnamed_project"){
       ProjectMetadata[eml_row,"project_name"] <- eml_data$name
     }
   }
  
  if("associatedReferences" %in% colnames(metadata)){
    warning("There are additional associated references in the data table that need to be checked by hand")
  }
  
  # ProjectMetadata $ created_on / updated_on
  ProjectMetadata$created_on <- rep(Sys.Date(), nrow(ProjectMetadata))
  ProjectMetadata$updated_on <- ProjectMetadata$created_on 
  
  # ProjectMetadata $ associated_media
  ### not in EML, most likely in a DwC metadata object, still need to figure out how to put this in...
  
  ### change references (if present) to IDs of Reference table
  if(any(ProjectMetadata$associated_references != "")){
    #------------------------------------------
    #### create Reference table ####
    #------------------------------------------
    # Reference $ id [ref: ProjectMetadata $ associated_references]
    # short author list = "et al." -version
    ReferenceTerms <- c("id", "authors_list", "doi", "short_authors", "title", "journal", 
                        "year", "occurences")
    Reference <- data.frame(matrix(nrow=0, ncol=length(ReferenceTerms), data=""), stringsAsFactors = FALSE)
    colnames(Reference) <- ReferenceTerms

    for(ref in ProjectMetadata$associated_references){
      if(ref != ""){
        # split different references, assume multiple references are separated by "|"
        all_refs <- c(unlist(strsplit(ref, "|", fixed=TRUE)))
        for(ref in all_refs){
          
          ######################
          ##
          ##
          ##   PART UNRESOLVED
          ##
          ##
          ######################
          print(ref)
          #### parsing text to citation is a major issue
          #see utils=> parse.citation(ref)
          ###

        }
      }
    }
  }



  #------------------------------------------
  #### create event_type table ####
  #------------------------------------------
  # event_type $ id [ref: ParentEvent $ event_type]
  event_type <- data.frame(matrix(nrow=length(unique(ParentEvent$event_type)), ncol=2, data=""), stringsAsFactors = FALSE)
  colnames(event_type) <- c("id", "name")

  event_type$name <- unique(ParentEvent$event_type)
  event_type$id <- 1:nrow(event_type)

  # change ParentEvent $ event_type to ID
  ParentEvent$event_type <- unlist(sapply(ParentEvent$event_type, 
                                       FUN = function(x){
                                         gsub(x,event_type[event_type$name==x,]$id,x)
                                       }))


  #------------------------------------------
  #### create Metadata table ####
  #------------------------------------------
  # Metadata $ id [ref: Event $ metadata]
  # Metadata $ sequence [ref: Sequences $ id]
  # Metadata $ env_biome [ref: Biome $ id]

  MetadataTerms <- c("id", "metadata_tag", "md_created_on", "metadata_creator", "license", 
                     "continent", "country", "state_province", "waterBody", "islandGroup", 
                     "island", "location", "geo_loc_name", "additional_info", "env_biome", 
                     "env_package", "env_feature", "env_material", "institutionID", 
                     "nucl_acid_amp", "nucl_acid_ext", "ref_biomaterial", "rel_to_oxygen", 
                     "rightsHolder", "samp_collect_device", "samp_store_dur", "samp_store_loc",
                     "samp_store_temp", "samp_vol_we_dna_ext", "samplingProtocol", 
                     "source_mat_id", "submitted_to_insdc", "investigation_type",
                     "isol_growth_condt", "lib_size", "sequence")
  Metadata <- data.frame(matrix(nrow=num_recs, ncol=length(MetadataTerms), data=""), stringsAsFactors = FALSE)
  colnames(Metadata) <- MetadataTerms
  
  # Metadata $ md_created_on
  Metadata$md_created_on <- rep(Sys.Date(), nrow(Metadata))
  
  # Metadata $ metadata_creator
  Metadata$metadata_creator <- rep(user_ID, nrow(Metadata))
  
  # Metadata $ license
  if("license" %in% colnames(metadata)){
    Metadata$license <- metadata$license
  }else{ # use CC BY 4.0 as default (Creative Commons BY attibution, version 4.0)
    Metadata$license <- rep("CC BY 4.0", num_recs)
  }
  
  ### geographical info
  ## remark: make this a hierarchical table, like Biome?
  # Metadata $ continent, country, state_province, waterBody, islandGroup, island, location
  geoTerms <- intersect(colnames(metadata), c("continent", "country", "state_province", 
                                              "waterBody", "islandGroup", "island", "locality"))
  if(length(geoTerms)!=0){
    for(gt in geoTerms){
      if(gt != "locality"){
        Metadata[,gt] <- metadata[,gt]
      }else{
        Metadata$location <- metadata$locality
      }
    }
  }
  # Metadata $ geo_loc_name
  if("geo_loc_name" %in% colnames(metadata)){
    Metadata$geo_loc_name <- metadata$geo_loc_name
  } else if(length(geoTerms)!=0){
    locName<-c()
    for(gt in geoTerms){
      locName <- paste(locName, ":", gt, "=", metadata[,gt], collapse="")
    }
    Metadata$geo_loc_name <- gsub("^:", "", locName)
  }
  
  # Metadata $ all the other terms
  for(varTerm in c("env_package", "env_feature", "env_material", "institutionID", 
                   "nucl_acid_amp", "nucl_acid_ext", "ref_biomaterial", "rel_to_oxygen", 
                   "rightsHolder", "samp_collect_device", "samp_store_dur", "samp_store_loc",
                   "samp_store_temp", "samp_vol_we_dna_ext", "samplingProtocol", 
                   "source_mat_id", "submitted_to_insdc", "investigation_type",
                   "isol_growth_condt", "lib_size")){
    if(varTerm %in% colnames(metadata)){
        Metadata[,varTerm] <- as.character(metadata[,varTerm])
        }
  }
  
  # Metadata $ id
  Metadata$id <- 1:nrow(Metadata)
  # Event $ metadata   !! asuming an event is equivalent to a sequence run (=sample)
  Event$metadata <- Metadata$id
  # Event $ metadata_exists
  Event$metadata_exists<- rep("TRUE", nrow(Event))
  

  # Metadata $ additional_info
  NoteTerms <- intersect(colnames(metadata), c("fieldNotes", "eventRemarks"))
  if(length(NoteTerms)!=0){
    Metadata$additional_info <- paste(metadata[,colnames(metadata) %in% NoteTerms], collapse=" | ")
  }
  
  #------------------------------------------
  #### create Sequences table ####
  #------------------------------------------
  # Sequences $ id [ref: Metadata $ sequence]
  
  SequnencesTerms <- c("id", "sequence_name", "MID", "subspecf_gen_lin", "target_gene",
                       "target_subfragment", "type", "primerName_forward", 
                       "primerName_reverse", "primer_forward", "primer_reverse", "run_type",
                       "seqData_url", "seqData_accessionNumber", "seqData_projectNumber",
                       "seqData_runNumber", "seqData_sampleNumber", "seqData_numberOfBases",
                       "seqData_numberOfSequences")
  
  Sequences <- data.frame(matrix(nrow=num_recs, ncol=length(SequencesTerms), data=""), stringsAsFactors = FALSE)
  colnames(Sequences) <- SequencesTerms
  
  # Sequences $ MID
  if("mid" %in% colnames(metadata)){
    if("additional_mid" %in% colnames(metadata)){
      Sequences$MID <- paste(metadata$mid, metadata$additional_mid, collapse="; ")
    }else{
      Sequences$MID <- metadata$mid
    }
  }
  
  # Sequences $ type
  if("library_source" %in% colnames(metadata)){
    Sequences$type <- metadata$library_source
  }
  
  # Sequences $ primerName_forward
  if("forward_primer" %in% colnames(metadata)){
    Sequences$primerName_forward <- metadata$forward_primer
  }
  # Sequences $ primerName_reverse
  if("reverse_primer" %in% colnames(metadata)){
    Sequences$primerName_reverse <- metadata$reverse_primer
  }
  
  # Sequences $ primer_forward
  if("primer_forward_sequence+forward" %in% colnames(metadata)){
    Sequences$primer_forward <- metadata$primer_forward_sequence
  }
  # Sequences $ primer_reverse
  if("primer_reverse_sequence" %in% colnames(metadata)){
    Sequences$primer_reverse <- metadata$primer_reverse_sequence
  }
  # Sequences $ url
  if("seqData_url" %in% colnames(metadata)){
    Sequences$url <- metadata$seqData_url
  }
  # Sequences $ seqData_accessionNumber | seqData_sampleNumber
  if("seqData_accessionNumber" %in% colnames(metadata)){
    Sequences$seqData_accessionNumber <- metadata$seqData_url
  }
  # Sequences $ seqData_accessionNumber | seqData_sampleNumber
  accTerms <- intersect(colnames(metadata), c("genbank_accession_numbers", "INSDC_SampleID"))
    if(length(accTerms)!=0){
      Sequences$seqData_accessionNumber <- paste(metadata[,colnames(metadata) %in% accTerms], collapse="; ")
      Sequences$seqData_sampleNumber <- Sequences$seqData_accessionNumber
    }
  # Sequences $ seqData_projectNumber
  if("bioproject" %in% colnames(metadata)){
    Sequences$seqData_projectNumber <- metadata$bioproject
  }
  # Sequences $ seqData_runNumber
  if("sra_run_number" %in% colnames(metadata)){
    Sequences$seqData_runNumber <- metadata$sra_run_number
  }
  # Sequences $ seqData_numberOfBases
  if("number_of_bases_predicted" %in% colnames(metadata)){
    Sequences$seqData_numberOfBases <- metadata$number_of_bases_predicted
  }
  # Sequences $ seqData_numberOfSequences 
  if("lib_reads_seqd" %in% colnames(metadata)){
    Sequences$seqData_numberOfSequences <- metadata$lib_reads_seqd
  }
  
  # remaining terms
  for(varTerm in c("subspecf_gen_lin", "target_gene",  "target_subfragment", "run_type")){
    if(varTerm %in% colnames(metadata)){
      Sequences[,varTerm] <- as.character(metadata[,varTerm])
    }
  }
  
  # Sequences $ id
  Sequences$id <- 1:nrow(Sequences)
  # Metadata $ sequence ## assuming every metadata record corresponds to a sequence record
  Metadata$sequence <- Sequences$id
  

  
  #------------------------------------------
  #### create Biome table ####
  #------------------------------------------
  # Biome $ id [ref: Metadata $ env_biome]
  
  BiomeTerms <- c("id", "name", "biome_level", "parent_biome")
  
  # Metadata $ env_biome
  
  
  
  #------------------------------------------
  #### create Environment table ####
  #------------------------------------------
  EnvironmentTerms <- c("id", "env_sample_name", "created_at", "Latitude", "Longitude",
                        "link_climate_info", "env_variable", "env_method", "env_units", 
                        "sequences", "env_numeric_value", "env_text_value")
  
  
  
  
  
  ### writing the data
  baseName <- "testdata"
  write.csv(Event, paste(dest.dir, "/", baseName, "_Event.csv", sep=""), na="", row.names = FALSE)
  write.csv(ParentEvent, paste(dest.dir, "/", baseName, "_ParentEvent.csv", sep=""), na="", row.names = FALSE)
  write.csv(Reference, paste(dest.dir, "/", baseName, "_Reference.csv", sep=""), na="", row.names = FALSE)
  write.csv(event_type, paste(dest.dir, "/", baseName, "_event_type.csv", sep=""), na="", row.names = FALSE)
  write.csv(Metadata, paste(dest.dir, "/", baseName, "_Metadata.csv", sep=""), na="", row.names = FALSE)
  write.csv(Sequences, paste(dest.dir, "/", baseName, "_Sequences.csv", sep=""), na="", row.names = FALSE)

  

  

  

  

  


  
  
  
  ## Event $ parent_sample (name) => add with parent event table
  
  ## Event $ occurence, environment
  
  ## Event $ occurence_exists, environment_exists
  
  # to be worked out with DwC:
  #if("occurenceID" %in% colnames(metadata)){
  #  Event_recs$occurrence<-metadata$occurenceID
  #  Event_recs$occurrence_exists<-rep(TRUE, num_recs)
  #} else{
  #  Event_recs$occurrence_exists<-rep(FALSE, num_recs)
  #}
  

  

  
  

  
  
  ############### done up to here  ###############
  ################################################
  
  
  
  
  
  
  
  
  
  


  # 6. Package table
  PackageTerms <- c("id", "name")
  
  
    if("env_package" %in% colnames(metadata)){
      Package_recs <- metadata$env_package
      ColsDone <- c(ColsDone, "env_package")
    } else if(metapackage != "multiple_packages"){
      Package_recs <- rep(metapackage, num_recs)
    }

  
  
  
  for(trm in setdiff(EventTerms, c("footprintWKT", "sample_name", "occurrence",
                                   "metadata",  "environment", "metadata_exists",
                                   "occurrence_exists", "environment_exists"))){
    if(trm %in% colnames(metadata)){
      Event_recs[,trm]<-metadata[,trm]
      ColsDone <- c(ColsDone, trm)
    }
  }
  
  
  
  # 7. Occurence table
  OccurenceTerms <- c("id", "occurrenceID", "taxon", "occurrence_notes", "occurrence_status",
                      "occurrence_class", "catalog_number", "date_identified", "other_catalog_numbers",
                      "recorded_by", "associated_sequences" )
  
  
  if(metaformat == "DwC"){
    Occurence_recs <- data.frame(matrix(nrow=num_recs, ncol=length(OccurenceTerms), data=NA))
    colnames(Occurence_recs) <- OccurenceTerms
    rownames(Occurence_recs) <- rownames(coredata)
    for(trm in setdiff(OccurenceTerms, "taxon")){
      if(trm %in% colnames(coredata)){
        Occurence_recs[,trm]<-coredata[,trm]
        ColsDone <- c(ColsDone, trm)
      }
    }
    
    ### taxon
    #####################
    ### TO BE RESOLVED:
    ### change taxon to taxonID?
    #####################
    
  }
  
  # 8. Taxa table
  TaxaTerms <- c("id", "name", "TaxonRank", "taxonID", "parent_taxa")
  
  
  if(metaformat == "MIxS"){
    if("subspecf_gen_lin" %in% colnames(metadata)){
      ### DO SOMETHING
    }
  }else if(metaformat == "DwC"){
    allRanks <- c("kingdom","phylum","class","order","family","genus","subgenus", 
                 "specificEpithet","infraspecificEpithet")
    
    allRelations <- data.frame(child=allRanks, parent=c("root", allRanks[-length(allRanks)]))
    
    txTerms <- as.character(polaaarLib[polaaarLib$polaaarDB_table=="Taxa",]$polaaarDB_name)

    
    txdata <- coredata[,intersect(colnames(coredata), allRanks)]
    
    tx_recs <- wideTab.to.hierarchicalTab(txdata, allRelations)
    
    ### before adding to the DB: check if not already in there!!!!
  }
  
 # split taxon into kingdom|phylum|class|order|family|genus|subgenus|specificEpithet|infraspecificEpithet
  
  as.character(polaaarLib[polaaarLib$polaaarDB_table=="Taxa",]$polaaarDB_name)
  Occurence_recs <- data.frame(matrix(nrow=num_recs, ncol=length(OccurenceTerms), data=NA))
  colnames(Occurence_recs) <- OccurenceTerms
  rownames(Occurence_recs) <- rownames(coredata)
  
  

   
  # 12. Units table
  UnitTerms <- c("id", "name", "html_tag")
  
  # 13. Variable table
  VariableTerms <- c("id", "name", "var_units", "method", "var_type")
  #var_type either TXT (text) or NUM (numeric)
  

  # 17. sampling_method table
  # links to variable table
  sampling_methodTerms <- c("id", "shortname", "description")

  

  

  
  
  
  

  