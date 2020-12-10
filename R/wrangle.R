#############################################################
#
# Functions for data-wrangling
#
############################################################

#' Wrangle input data for format required by simulation
#'
#' @param inndata data-frame with inndata comming from "get_inndata"
#' @param start_year Year of staring recording as introductions
#' @param end_year Year of ending recording as introductions
#' @param geoselect_native A spatial.polygon.dataframe indicating the native distribution range of the species in epsg:4326
#' @param focal_species The latin (canonical) name of the species
#' @return A data.frame
#' @examples
#'
#' setwd("~/temp") # simulation working directory
#' # Download example input files
#' get_inndata(serveradress="my-server-adress.no",datafolder="./data") # will take som time to load
#' download.file("https://ntnu.box.com/shared/static/kzgmsa898gcks57iacogpet1g2yue8ib.rds",
#'              "./data/geoselect_native_Esox_lucius.rds") # not yet procedures to draw this rigth from DB
#'
#' inndata <- readRDS("./data/view_occurrence_by_event.rds")
#' geoselect_native <- readRDS("./data/geoselect_native_Esox_lucius.rds")
#' focal_species <- "Esox lucius"
#' mynewdata_aslist <- wrangle_and_slice(inndata=inndata,start_year=1970,end_year=2017)
#' mynewdata <- mynewdata_aslist[["data"]]
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import lubridate
#' @import data.table
#' @import FNN
#' @import sp
#'
#' @export
#'
#'

wrangle_and_slice <- function(start_year,end_year,inndata,focal_species,geoselect_native){

  #######################################################
  # Reclassify establishMeans of focal species ----------------------------
  # based upon distr. polygon
  #######################################################
  inndata_species <- inndata[inndata$scientificName==focal_species,]

  # load spatial filter (polygon) and set projection to epsg:4326
  latlong = "+init=epsg:4326"
  geoselect_spdf <- spTransform(geoselect_native,latlong)

  # convert inndata to SpatialPointDataFrame
  #inndata<-outdata_timeslot
  inndata_species <- as.data.frame(inndata_species)
  coordinates(inndata_species) = ~decimalLongitude + decimalLatitude

  # project outdata_event
  latlong = "+init=epsg:4326"
  proj4string(inndata_species) = CRS(latlong)

  # filter out values from outdata_event based upon geoselect_spdf
  # Note: clipping spatial data in R by using [] from http://robinlovelace.net/r/2014/07/29/clipping-with-r.html
  inndata_geoselect <- inndata_species[geoselect_spdf, ]

  inndata_species$establishmentMeans <- ifelse(inndata_species$waterBodyID %in% inndata_geoselect@data$waterBodyID,"native",
                                       "introduced")
  inndata_species <- as.data.frame(inndata_species)
  inndata_species <- inndata_species %>% select(eventID,establishmentMeans_tmp=establishmentMeans)

  inndata <- left_join(inndata,inndata_species,by="eventID")
  inndata$establishmentMeans <- ifelse(inndata$scientificName==focal_species,
                                       as.character(inndata$establishmentMeans_tmp),inndata$establishmentMeans)

  #######################################################
  # spread data to wide format by species -----------------------------
  #######################################################
  inndata <- inndata %>% mutate(scientificName2=str_replace(scientificName, " ", "_"))
  inndata$year <- year(parse_date_time(inndata$dateEnd,"ymd"))
  inndata <- inndata[complete.cases(inndata$year),]

  spread_data <- inndata %>% dplyr::select(eventID,scientificName2)
  spread_data$presence <- 1
  spread_data <- spread_data %>% spread(scientificName2,presence,fill=0)

  # # select variables to be included. NB! No occurrence level variables, only event level and above
  # inndata2 <- inndata %>% dplyr::select(dateEnd,year,eventID,scientificName,
  #                                       recordedBy,samplingTaxaRange,samplingProtocol,
  #                                       datasetID,county,municipality,countryCode,eventID,#decimalLatitude,decimalLongitude,
  #                                       utm_x,utm_y,waterBodyID,area_km2,distance_to_road,
  #                                       perimeter_m,elevation,eurolst_bio01,eurolst_bio02,
  #                                       eurolst_bio03,eurolst_bio05,eurolst_bio06,eurolst_bio07,
  #                                       eurolst_bio10,eurolst_bio11) %>%
  #   distinct()
  inndata2 <- left_join(inndata,spread_data,by="eventID")



  #######################################################
  # Get distance to closest pop at start of time slot----
  #######################################################

  # Create data.frame with all introductions of focal-species
  # within the given time-slot given as introduction=1,
  # distance_to_closest population of focal species
  # given as distance at start of time-slot

  # Select locations of all populations of focal species at start of timeslot i.
  # Useatd for calcualation distance to closest population
  data1 <- inndata2 %>%
    filter(year<=start_year,scientificName==focal_species) %>%
    dplyr::select(utm_x,utm_y,waterBodyID) %>% distinct()

  # Select all locations for which there where observations of fish
  data2 <- inndata2 %>%
    filter(year<=end_year) %>%
    dplyr::select(utm_x,utm_y,waterBodyID) %>% distinct()


  # Calculate distance to closest population of focal species to any lake in dataset
  # at start of time-slot for all locations with fish observations at end
  # of time-slot i. Use get.knnx from the FNN package.
  nn <- get.knnx(data1[c("utm_x","utm_y")],data2[c("utm_x","utm_y")],2)

  dist_to_closest_pop <- ifelse(nn$nn.dist[,1]==0,nn$nn.dist[,2],nn$nn.dist[,1])
  waterBodyID <- data2$waterBodyID
  distance_data <- as.data.frame(cbind(dist_to_closest_pop,waterBodyID))
  # run average distance for each waterBodyID (the problem of duplicated
  # waterBodyIDs here comes from input coords may resemble both outlet
  # or centroid of waterbody)
  distance_data <- distance_data %>%
    group_by(waterBodyID) %>%
    summarize(dist_to_closest_pop=mean(dist_to_closest_pop))

  # match back with inndata
  inndata3 <- left_join(inndata2,distance_data,by="waterBodyID")

  #########################################################################
  # Add 0/1 column of introductions occuring during timeslot -------------
  #########################################################################
  # In order for this to classify as an introduction event,
  # we first need to check out that the species not have been
  # observed in the given waterbody before.

  wb_fish_j <- inndata3 %>% filter(scientificName==focal_species,year<=start_year) %>%
    dplyr::select(waterBodyID)

  # Then, select events with observation of focal species in time-slot,
  # where focal species are classified as introduced.
  outdata_temp2 <- inndata3 %>%
    dplyr::filter(scientificName==focal_species,
                                             year<=end_year,
                                             year>=start_year,
                                             establishmentMeans=="introduced",
                                             !waterBodyID %in% wb_fish_j$waterBodyID) %>%
    dplyr::select(eventID)

  # Then classify events as introduction events of species j.
  inndata3$introduced <- ifelse(inndata$eventID %in% outdata_temp2$eventID,
                                    1,0)

  inndata3$focal_species <- focal_species


  #########################################################################################
  # Group by waterBodyID and species -----
  #######################################################################################

  inndata4 <- inndata3 %>%
    group_by(waterBodyID,county,municipality,countryCode) %>%
    summarise(dist_to_closest_pop=mean(dist_to_closest_pop))

  # get presence/absence aggregated for each time-slot,waterBody and focal species
  # and join to outdata_timeslot

  # select out the numeric variables from innndata to aggregate on by waterbody
  tmp <- inndata3 %>% select_if(is.numeric)

  aggdata2 <- tmp %>% group_by(waterBodyID) %>% summarise_all(funs(max))
  aggdata2$focal_species <- focal_species

  # remove all presence/absence data for species not recorded in selected area
  # (to clean up outdata)
  specieslist2 <- unique(inndata3$scientificName2)
  colsToRemove <- specieslist2[colSums(aggdata2[specieslist2])==0]
  aggdata2 <- aggdata2[!names(aggdata2) %in% colsToRemove]
  aggdata2 <- aggdata2 %>% select(-dist_to_closest_pop)

  # link back to waterbody etc...
  aggdata3 <- left_join(inndata4,aggdata2,by="waterBodyID")

  #######################################
  # Add some derived variables on the fly -------------------------
  #######################################
  #make SCI, and index portraying shore line complexity (see Voje etal. 2013 for similar approach)

  aggdata3$SCI<-(aggdata3$perimeter_m/1000)/(2*sqrt(pi*aggdata3$area_km2))
  aggdata3$area_km2_log<-log(aggdata3$area_km2)
  aggdata3$introduced<-as.integer(aggdata3$introduced)
  aggdata3$dist_to_closest_pop_log<-log(aggdata3$dist_to_closest_pop)
  aggdata3$county<-as.factor(aggdata3$county)
  aggdata3$distance_to_road<-aggdata3$distance_to_road+0.1
  aggdata3$distance_to_road_log<-log(aggdata3$distance_to_road)

  outdata <- list()
  outdata[["data"]] <- aggdata3
  outdata[["start_year"]] <- start_year
  outdata[["end_year"]] <- end_year
  outdata[["focal_species"]] <- focal_species

  return(outdata)

} # end of function

