##################################################################################
#
# Build scenarios for species introductions
# based upon new introductions (as predicted from model)
# and secondary spread as predicted using connectivity matrix
#
# NB! Now runs from Rproj home directory (translocation main folder). Setting relative paths accordingly
#
# coexistence pike/browntrout: http://rspb.royalsocietypublishing.org/content/281/1775/20132641
# migration barriers pike: http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2007.01382.x/full
#
#################################################################################
# load libraries 
library(sp)
library(dplyr)
library(data.table)

#...................................
# 1. load input data and funcs ------
#...................................

# load inndata
# outdata_timeslot_Rutilius_rutilus.rds
tmpdata <- tempfile()
download.file(url="https://ntnu.box.com/shared/static/1dgrjwozx4r6zcprta8oj35ocjtdj2nl.rds",
              destfile=tmpdata)
inndata_timeslot <- readRDS(tmpdata)

# load functions
source("./R/f_predict_introduction_events.R")
source("./R/f_calc_distance.R")
source("./R/f_postsim_processing.R")

# load model object used for predictions
# Mort.... 
tmpdata2 <- tempfile()
download.file(url="https://ntnu.box.com/shared/static/ftt7twf3fhdbm4izei0qktj382do3bx8.rds",
              destfile=tmpdata2)
brt_mod <- readRDS(tmpdata2)


# load connectivity matrix for a given area (this defines the geographic scope of simulation)
tmpdata3 <- tempfile()
download.file(url="https://ntnu.box.com/shared/static/2odyiz3nx5fcrkna3ztwbkj2ge01w59y.rds",
              destfile=tmpdata3)
connectivity <- readRDS(tmpdata3)
#connectivity <- readRDS("./Data/connectivity_troendelag.rds") # load model object 
# convinience renaming and selecting
connectivity <- connectivity %>% 
  rename(waterBodyID=lakeid) %>%
  select(waterBodyID,downstream_lakes,downstream_stream_length_km,
         upstream_lakes,upstream_lakes_slope_max,upstream_lakes_slope_max_max,
         upstream_lakes_distance_km)

#.........................................................................................
# Define input variables ----
#.........................................................................................

# Select inndata geographic range based upon connectivity matrix extent
inndata <- inndata_timeslot[inndata_timeslot$waterBodyID %in% connectivity$waterBodyID,]

# species specific stuff
species <- "Rutilus rutilus"
species_var <- stringr::str_replace(species," ","_") # variable describing present/absence of focal species
temp_inc <- 0 # temperature increas

# simulation and time specific stuff
Nsims <- 500 # number of iterations
sim_duration <- 50 # Duration of the scenario, in years (will be corced to multiple of time_slot_length)
time_slot_length <- 5 # Duration of the time-slot, in years (time-span of each time-slot)
gmb_time_slot_length <- 30 # Duration of the time-slot, in years, used to estimate establishment probability
n_time_slots <- as.integer(sim_duration/time_slot_length)
start_year <- 2017

# secondary dispersal stuff
with_secondary <- TRUE # should secondary spread be included in simulations?
slope_barrier <- 0.001 # max slope for migration upstream

# Before each simulation run!!!! 
# Create new dataframe / vars for simulation bookkeeping.
# Use latest time-slot (if multiple) in inndata
inndata_sim1 <- inndata[inndata$t_slot==unique(inndata$t_slot)[1],]

# Exterminate all pressent populations in VFO Trondelag.... 
inndata_sim1[species_var] <- 0
#...................................................................
# Run simulation ----------------
# Output is a data.frame listing waterBodyID of each introduction
# and a true/false if this was a primary introduction or result of 
# secondary dispersal. Also giving simulation run (j) and time-slot (i)
# as index.
#...................................................................

# j simulation runs...
for(j in 1:Nsims){
  
  inndata_sim <- as.data.frame(inndata_sim1) # reset dataset for each simulation run
  
  # simulate across time periods j
  # to run in parallelization, see: https://www.r-bloggers.com/parallel-r-loops-for-windows-and-linux/
  for(i in 1:n_time_slots){
    
    ### i.1 predict translocations and store new introductions in temp object
    tmp_trans <- f_predict_introduction_events_gmb(inndata_sim,brt_mod,species,temp_inc,time_slot_length,gmb_time_slot_length)
    
    # include secondary dispeersal?
    if(with_secondary==TRUE){
      
      # get wbID from introductions in run i
      introduction_lakes <- tmp_trans[tmp_trans$introduced==1,]$waterBodyID
      
      #.............................................................
      # Downstream dispersal
      #.............................................................
      # get vector of wbID for downstream lakes to those with introduction in simrun i
      downstream_lakes <- connectivity$downstream_lakes[connectivity$waterBodyID %in% introduction_lakes]
      downstream_lakes <- paste(downstream_lakes,sep=",",collapse=",")
      downstream_lakes <- na.omit(as.numeric(unlist(strsplit(downstream_lakes,split=","))))
      downstream_lakes <- unique(downstream_lakes) # introductions only listed once (remove duplicated lakeIDs)
      
      # select out downstream lakes that does not have pike at start of time-slot
      downstream_lakes <- downstream_lakes[!(downstream_lakes %in% inndata_sim$waterBodyID[inndata_sim[species_var]==1])]
      
      # finally assign introduction to downstream lakes (without previous obs/intro)
      tmp_trans$introduced <- ifelse(tmp_trans$waterBodyID %in% downstream_lakes,1,tmp_trans$introduced)
      
      #.............................................................
      # Upstream dispersal - NB! Check this part.... unequal length of upstream_lakes and upstream_slopes vector!!!!
      #.............................................................
      # get vector of wbID for upstream lakes 
      upstream_lakes <- connectivity$upstream_lakes[connectivity$waterBodyID %in% introduction_lakes]      
      upstream_lakes <- paste(upstream_lakes,sep=",",collapse=",")
      upstream_lakes <- as.numeric(unlist(strsplit(upstream_lakes,split=",")))
      # get vector of upstream slopes 
      upstream_slopes <- connectivity$upstream_lakes_slope_max_max[connectivity$waterBodyID %in% introduction_lakes]
      upstream_slopes <- paste(upstream_slopes,sep=",",collapse=",")
      upstream_slopes <- unlist(strsplit(upstream_slopes,split=","))
      upstream_slopes <- gsub(" ","",upstream_slopes)
      upstream_slopes <- as.numeric(upstream_slopes)
      
      # finally reachable lakes and assign introduction 
      upstream_lakes_reachable <- na.omit(upstream_lakes[upstream_slopes<slope_barrier])
      upstream_lakes_reachable <- unique(upstream_lakes_reachable) # introductions only listed once (remove duplicated lakeIDs)
      
      # add upstream_lake intros to introduced vector
      tmp_trans$introduced <- ifelse(tmp_trans$waterBodyID %in% upstream_lakes_reachable,1,tmp_trans$introduced)
      
    } # end of secondary==TRUE if statement
    
    ### Store output from time-period i, simulationrun j
    
    # create variables to store
    intro <- tmp_trans$waterBodyID[tmp_trans$introduced==1]
    intro_is_secondary <- ifelse(intro %in% downstream_lakes,TRUE,FALSE)
    time_slot_i <- rep(i,length(intro))
    sim_j <- rep(j,length(intro))
    start_year_i <- rep(( start_year+((i-1)*time_slot_length) ),
                        length(intro))
    end_year_i <- rep(( start_year+(i*time_slot_length)-1 ),
                      length(intro))
    tmp_output <- data.frame(intro=intro,
                             intro_is_secondary=intro_is_secondary,
                             time_slot_i=time_slot_i,
                             sim_j=sim_j,
                             start_year_i=start_year_i,
                             end_year_i=end_year_i)
    
    # on first sim-run, create new object to store output
    if(i==1 & j==1){
      sim_output <- data.frame()
    }
    # append sim_output with data from time-slot i, sim j 
    sim_output <- bind_rows(sim_output,tmp_output)
    
    ### modify inndata_sim with new introductions to provide innput to time_slot i+1
    # Add new introductions to "species_var" column in inndata_sim
    tmp1 <- inndata_sim[species_var]
    inndata_sim[species_var] <- ifelse(tmp_trans$introduced==1,1,tmp1[,1])
    
    # recalculate distance to closest population and replace values
    # in inndata_sim where distance is smaller than previous; 
    # i.e. accounting for situations where closest population is outside
    # the geographic area beein investigated
    # NB! need a try function here to make run in cases where there are non intros in simrun
    # and f_calc_dist fail
    tmp <- try(f_calc_dist(outdata=inndata_sim,species=species),silent = TRUE)
    if(!is.null(dim(tmp))){
      inndata_sim$dist_to_closest_pop <- ifelse(tmp$dist_to_closest_pop<inndata_sim$dist_to_closest_pop,
                                                tmp$dist_to_closest_pop,
                                                inndata_sim$dist_to_closest_pop)
    } else {
      inndata_sim$dist_to_closest_pop <-  inndata_sim$dist_to_closest_pop
    }
    
    tmp <- data.frame()

    # and log variable for model predictions
    inndata_sim$dist_to_closest_pop_log <- log(inndata_sim$dist_to_closest_pop)
    
  } # end of i loop
  
  # display progress
  print(paste("sim-run: ",j,", out of: ",Nsims))
  
} # end of j loop



#............................................................................
# sum up sim_output pr lake and save / write to database  
#............................................................................
# store inndata1, raw and lake aggregated sim_output as list in rds object 
sim_output_lake <- f_sim_output_lake(sim_output,inndata_sim1,Nsims,n_time_slots)

tmpout <- list()
tmpout[["sim_output"]] <- sim_output
tmpout[["inndata_sim1"]] <- inndata_sim1
tmpout[["sim_output_lake"]] <- sim_output_lake
tmpout[["focal_species"]] <- species
tmpout[["time_slot_length"]] <- time_slot_length
tmpout[["start_year"]] <- start_year

# Write output to local disk
url <- paste("./Data/sim_out_",species_var,"_extermination_scenario.rds",sep="")
saveRDS(tmpout,url)

# write output to BOX
library(boxr) 
box_auth() # When run on server, create OAuth2.0 token locally and copy to server - see boxr vignett (remember .gitignore if needed)
box_ul(dir_id = 29103527607, file = paste("./Data/sim_out_",species_var,".rds",sep=""))

# Write lake-specific summary to database 
dataToWrite <- tmpout
nameOfTable <- tolower(paste("sim_output_",species_var,"_extermination_scenario",sep=""))
f_write_simresult_to_db(dataToWrite=sim_output_lake,nameOfTable)

