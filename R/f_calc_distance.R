#############################
# Calc. dist to closest pop #
#############################

f_calc_dist <- function(outdata,species){
  require(dplyr)
  require(FNN)
  species2=str_replace(species, " ", "_")
  # Select locations of all populations of focal species at start of timeslot i.
  # Used for calcualation distance to closest population
  data1 <-as.data.frame(outdata[outdata[[species2]] == 1,] %>%
    dplyr::select(utm_x,utm_y,waterBodyID) %>% distinct())

  # Select all locations for which there where observations of
  # fish at the end of time-slot i.
  data2 <- outdata %>%
    dplyr::select(utm_x,utm_y,waterBodyID) %>%
    distinct() %>%
    as.data.frame()

  # Calculate distance to closest population of focal species -------
  # at start of time-slot for all locations with fish observations at end
  # of time-slot i. Use get.knnx from the FNN package.
  nn <- get.knnx(data1[c("utm_x","utm_y")],data2[c("utm_x","utm_y")],2)

  dist_to_closest_pop <- ifelse(nn$nn.dist[,1]==0,nn$nn.dist[,2],nn$nn.dist[,1])
  waterBodyID <- data2$waterBodyID
  distance_data <- as.data.frame(cbind(dist_to_closest_pop,waterBodyID))
  #outdata <- left_join(outdata,distance_data,by="waterBodyID")
  return(distance_data)
}
