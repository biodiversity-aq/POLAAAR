###library() the polartools package for the QC functions

polaaarLib<-read.csv("/Users/msweetlove/OneDrive_RBINS/mARS_NewSeqData/polaRRTools/Data/PolaaarLibrary.csv", 
                                 header=TRUE)



write.metadata.MIxS.as.mars <- function(metadata.object, name.prefix=NULL, out.dir=getwd(),
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
  if(!dir.exists(file.path(out.dir))){
    stop("Could not find the directory provided in the dest.dir argument.")
  }
  mimarks.name <- paste(out.dir, "/MiMARKS", name.prefix, ".csv", sep="")
  seqset.name <- paste(out.dir, "/SeqSet", name.prefix, ".csv", sep="")
  
  # 3. create MiMARKS tabble
  # 3.1 first row
  Row01 <- paste("section", "Structured Comment Name", "units", 
                 paste(rownames(metadata), collapse=","), sep=",")
  # 3.2 general terms
  mimarksTerms <- as.character(MarsLib[MarsLib$mars_mimarks==TRUE,]$name) #the basic terms that need to be in every dataset
  #associated units
  mimarksUnits <- as.character(TermsLib[match(setdiff(mimarksTerms, names(metaunits)), TermsLib$name),]$expected_unit)
  names(mimarksUnits) <- setdiff(mimarksTerms, names(metaunits))
  mimarksUnits <- c(mimarksUnits, metaunits)
  #assocated MIxS sections
  mimarksSection <- as.character(TermsLib[match(setdiff(mimarksTerms, names(metasection)), TermsLib$name),]$section)
  names(mimarksSection) <- setdiff(mimarksTerms, names(metasection))
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
    Row0x_data <- gsub("NA", "", Row0x_data)
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
      Row0x_data <- gsub("NA", "", Row0x_data)
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
      Row0x_data <- metadata[,colnames(metadata)==mimarksEquiv]
      if(mimarksEquiv == "investigation_type"){ ###special case
        Row0x_data <- gsub("mimarks-survey", "marker gene", Row0x_data)
      }
      if(ncol(Row0x_data)==0){ #case it had no mimarks equivalent or wasn't in the colnames of metadata
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
      }
    }
    Row0x_data[is.na(Row0x_data)]<-""
    Row0x_data <- gsub("NA", "", Row0x_data)
    Row0x_data <- paste(Row0x_data, collapse=",")
    Row0x <- paste(Row0x_name, Row0x_data, sep=",")
    Row02 <-  paste(Row02, Row0x, sep="\n")
  }
  seqset.data <- paste(Row01, Row02, sep="")
  
  # 5. finalize
  cat(paste("The data had been written to ", out.dir, "\n",sep=""))
  
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

get.polaaar.EML.data <- function(EML.url){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description 
  #' @param EML.url
  #' @details 
  #' @return a list of the data needed for polaaar

  require(EML)
  require(emld)
  
  #EML.url <- "https://api.gbif.org/v1/dataset/ebe73e19-11eb-48d2-9263-fb0caa3b7b5a/document"
  
  eml_file <- emld::as_emld(EML.url, from = "xml")
  
  # abstract
  abstract <- eml_get(eml_file,"abstract", from="xml")$para
  
  # start_date | end_date
  dates <- eml_get(eml_file,"temporalCoverage", from="xml")
  if("rangeOfDates" %in% names(dates)){
    start_date  <- gsub('\n| |\t', '', dates$rangeOfDates$beginDate)
    end_date  <- gsub('\n| |\t', '', dates$rangeOfDates$endDate)
  } else if("singleDateTime" %in% names(dates)){
    start_date  <- gsub('\n| |\t', '', dates$singleDateTime)
    end_date <- start_date
  }else{
    print(names(dates))
  }
  
  # bounding_box
  boundingBox <- unlist(EML::eml_get(eml_file, "boundingCoordinates", from = "xml"))
  boundingBox <- paste("SRID=4326;POLYGON ((",
                       boundingBox["northBoundingCoordinate"], ", ",
                       boundingBox["eastBoundingCoordinate"], ", ",
                       boundingBox["southBoundingCoordinate"], ", ",
                       boundingBox["westBoundingCoordinate"], "))", sep="")
  
  # is_public
length(EML::eml_get(eml_file, "pubDate", from = "xml"))
class(EML::eml_get(eml_file, "dataset", from = "xml"))


  # associated_references
  # associated_media

  
  
}




metadata.to.polaaar <- function(metadata.object, EML.url=NA, user_ID=""){
  #' @author Maxime Sweetlove ccBY 4.0 2019
  #' @description metadata.to.polaaar writes a metadata.MIxS or metadata.DwC object to the polaaar database
  #' @param metadata.object a metadata.MIxS class object or metadata.DwC. The object to be formatted to fit the polaaar database schema
  #' @param EML.url a named vector. URL to the EML (IPT) associated with the dataset, written as a named vector, e.g. c(project1="http://url1", project2="http://url1"). The names in the vector must correspond to the different projects listed under "project_name" in the metadata object. If only one link is provided and is unnamed, than this will automatically be the parent project to all the samples in the metadata object.
  #' @param user_ID ID of the user that send the data
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
  # project_creator = user id of person that emailed the spreadsheet
  # boundingbox notation:
  # SRID=4326;POLYGON ((N, E, S, W))
  ProjectMetadataTerms <- c("id", "project_name", "start_date", "end_date", "EML_URL", "abstract", 
                            "bounding_box", "is_public", "associated_references", 
                            "associated_media", "created_on", "updated_on", 
                            "project_creator", "project_qaqc")
  ProjectMetadata <- data.frame(matrix(nrow=length(EML.url), ncol=length(ProjectMetadataTerms), data=""))
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
  ProjectMetadata[,colnames(ProjectMetadata %in% c("abstact"))] <- as.character(ProjectMetadata[,colnames(ProjectMetadata %in% c("abstact"))])
  for(emlurl in EML.url){
    eml_row <- as.numeric(row.names(ProjectMetadata[match(emlurl,ProjectMetadata$EML_URL),]))

    
    }

  
  # ProjectMetadata $ created_on / updated_on



  
  
  
  
  
  ### change ParentEvent $ event_type to ID
  


  ## Event $ parent_sample (name) => add with parent event table
  
  ## Event $ occurence, metadata, environment
  
  ## Event $ occurence_exists, metadata_exists, environment_exists
  
  # to be worked out with DwC:
  #if("occurenceID" %in% colnames(metadata)){
  #  Event_recs$occurrence<-metadata$occurenceID
  #  Event_recs$occurrence_exists<-rep(TRUE, num_recs)
  #} else{
  #  Event_recs$occurrence_exists<-rep(FALSE, num_recs)
  #}
  

  

  
  

  
  
  ############### done up to here  ###############
  ################################################
  
  
  
  
  
  
  
  
  
  
  # 2.2. EventType table
  #*****************************
  # id = [ParentEvent]:event_type
  EventTypeTerms <- c("id", "name")
  if("investigation_type" %in% colnames(metadata)){
    types <- unique(metadata$investigation_type)
    EventType_recs <- data.frame(matrix(nrow=length(types), ncol=length(EventTypeTerms), data=NA))
    colnames(EventType_recs) <- EventTypeTerms
    EventType_recs$id <- 1:length(types)
    EventType_recs$name <- types
  }
  
  

  

  

  
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
  
  
  # 9. Biome table
  BiomeTerms <- c("id", "name", "biome_level", "parent_biome")
  
  # 10. Environment table
  EnvironmentTerms <- c("id", "env_sample_name", "created_at", "Latitude", "Longitude",
                        "link_climate_info", "env_variable", "env_method", "env_units", 
                        "sequences", "env_numeric_value", "env_text_value")
  
  # 11. Metadata table
  
  #### note: make geographic location also recursive table as biome
  MetadataTerms <- c("metadata_tag", "md_created_on", "metadata_creator", "license", 
                     "continent", "country", "state_province", "waterBody", "islandGroup", 
                     "island", "location", "geo_loc_name", "additional_info", "env_biome", 
                     "env_package", "env_feature", "env_material", "institutionID", 
                     "nucl_acid_amp", "nucl_acid_ext", "ref_biomaterial", "rel_to_oxygen", 
                     "rightsHolder", "samp_collect_device", "samp_store_dur", "samp_store_loc",
                     "samp_store_temp", "samp_vol_we_dna_ext", "samplingProtocol", 
                     "source_mat_id", "submitted_to_insdc", "investigation_type","isol_growth_condt",
                     "lib_size", "sequence")

  
  # 12. Units table
  UnitTerms <- c("id", "name", "html_tag")
  
  # 13. Variable table
  VariableTerms <- c("id", "name", "var_units", "method", "var_type")
  #var_type either TXT (text) or NUM (numeric)
  
  # 14. Sequnences table
  SequnencesTerms <- c("sequence_name", "MID", "subspecf_gen_lin", "target_gene",
                       "target_subfragment", "type", "primerName_forward", 
                       "primerName_reverse", "primer_forward", "primer_reverse", "run_type",
                       "seqData_url", "seqData_accessionNumber", "seqData_projectNumber",
                       "seqData_runNumber", "seqData_sampleNumber", "seqData_numberOfBases",
                       "seqData_numberOfSequences")


  
  # 16. Reference table
  ReferenceTerms <- c("id", "authors_list", "doi", "short_authors", "title", "journal", 
                      "year", "occurences")
  #short author list = et al.
  
  # 17. sampling_method table
  # links to variable table
  sampling_methodTerms <- c("id", "shortname", "description")

  

  

  
  
  
  

  
  
  # 3. create MiMARKS tabble
  # 3.1 first row
  Row01 <- paste("section", "Structured Comment Name", "units", 
                 paste(rownames(metadata), collapse=","), sep=",")
  # 3.2 general terms
  mimarksTerms <- as.character(MarsLib[MarsLib$mars_mimarks==TRUE,]$name) #the basic terms that need to be in every dataset
  #associated units
  mimarksUnits <- as.character(TermsLib[match(setdiff(mimarksTerms, names(metaunits)), TermsLib$name),]$expected_unit)
  names(mimarksUnits) <- setdiff(mimarksTerms, names(metaunits))
  mimarksUnits <- c(mimarksUnits, metaunits)
  #assocated MIxS sections
  mimarksSection <- as.character(TermsLib[match(setdiff(mimarksTerms, names(metasection)), TermsLib$name),]$section)
  names(mimarksSection) <- setdiff(mimarksTerms, names(metasection))
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
    Row0x_data <- gsub("NA", "", Row0x_data)
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
      Row0x_data <- gsub("NA", "", Row0x_data)
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
      Row0x_data <- metadata[,colnames(metadata)==mimarksEquiv]
      if(mimarksEquiv == "investigation_type"){ ###special case
        Row0x_data <- gsub("mimarks-survey", "marker gene", Row0x_data)
      }
      if(ncol(Row0x_data)==0){ #case it had no mimarks equivalent or wasn't in the colnames of metadata
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
      }
    }
    Row0x_data[is.na(Row0x_data)]<-""
    Row0x_data <- gsub("NA", "", Row0x_data)
    Row0x_data <- paste(Row0x_data, collapse=",")
    Row0x <- paste(Row0x_name, Row0x_data, sep=",")
    Row02 <-  paste(Row02, Row0x, sep="\n")
  }
  seqset.data <- paste(Row01, Row02, sep="")
  
  # 5. finalize
  cat(paste("The data had been written to ", out.dir, "\n",sep=""))
  
  write.table(mimarks.data, file=mimarks.name, 
              col.names = FALSE, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  write.table(seqset.data, file=seqset.name, 
              col.names = FALSE, row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  }
