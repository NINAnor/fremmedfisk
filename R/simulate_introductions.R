#!/usr/bin/env Rscript

#########################
#      brt-model        #
#########################

#' Genneral information about BRT-models (Boosted Regression Trees)
#' Brt-models are a type of regression model technique which is very flexible. It thus has some benefits when fitting ecological data, which is usually non-linear and so on.
#' Given care in the model-fittingg, brt can give predictive advantages over methods as e.g. glm or gam.
#' The following script uses the dismo and gmb packages to optimise a brt model for the given species.
#' Analytically, BRT regularization involves jointly optimizing the number of  trees, learning rate and tree complexity.
#' Here optimizing numbers of trees is done through the gbm.step function, and the script includes a function that tries to aid in the optimazation prosess for learning rate and tree complexity (the get.train.diganostic.func).
#' The next step is to fit the actual model used for predictions (brt_mod)
#'  Its is recomended to read through "A working guide to boosted regression trees" (Elith, et al 2009), before atempting your first go.


if (!require("optparse", character.only=T, quietly=T)) {
  install.packages("optparse")
  require("optparse")
}
library("optparse")

option_list = list(
  make_option(c("-s", "--species_group"), type="character", default=NULL, 
              help="Species group to process", metavar="character"),
  make_option(c("-w", "--workdir"), type="character", default="/data/R/Prosjekter/13845000_invafish/", 
              help="Name of the database to use", metavar="character"),
  make_option(c("-S", "--slope_measure"), type="character", default="slope_7_maximum", 
              help="Slope measure to use in connectivity analysis", metavar="character"),
  make_option(c("-t", "--slope_thresholds"), type="character", default="600,700,800", 
              help="Slope threshods used as estimate for dispersal barriers", metavar="character"),
  make_option(c("-n", "--nsims"), type="integer", default=500, 
              help="Number of simulations", metavar="integer"),
  make_option(c("-N", "--scenario_name"), type="character", default="scenario_0", 
              help="Name of the scenario", metavar="character"),
  make_option(c("-p", "--percentage_extermination"), type="double", default=0.0, 
              help="Percentage of populations to simulate extermination for", metavar="double"),
  make_option(c("-e", "--exterminate_waterbodyid"), type="character", default=NULL, 
              help="List of waterbodies to simulate extermination for", metavar="character"),
  make_option(c("-f", "--fish_barrier"), type="character", default=NULL, 
              help="List of waterbodies to simulate a fish barrier for", metavar="character"),
  make_option(c("-b", "--below"), type="logical", default=FALSE, 
              help="Simulate fish barrier below listed waterbodies [default is above, TRUE means below]", metavar="logical"),
  make_option(c("-T", "--temperature_increase"), type="double", default=0.0, 
              help="Temperature increase in simulation in degree celsius", metavar="double")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# To run in IDE/RStudio skip option parsing and use this opt list
# opt <- list("species_group"="Agn", "workdir"="/data/R/Prosjekter/13845000_invafish/", "nsims"=100, "slope_measure"="slope_7_maximum", "slope_thresholds"="600,700,800")

# To clarify
# - "start populations" kun introduserte eller alle, hva med ukjente
# - Hvordan håndtere forekomst av flere arter fra en gruppe i et vann
#   Relevant i f_predict_introduction_events_gmb
# - Hvordan håndtere sampling bias i rom og tid (begge er ofte skev fordelt)
#   Relevant i f_predict_introduction_events_gmb og kanskje også bygging av gbm
# Parallelize simulations (they are slow (mainly predictions with many lakes)
# Split script into three parts (1, data collection/wrangling, 2. model building,
#   3. simulation and writing of data) so it can be run on the command line to avoid timeouts

# Run with arguments:
# Rscript --vanilla build_brt.R -s "Agn"



if (!require('devtools', character.only=T, quietly=T)) {
  # Requires OpenSSL on the system
  install.packages('devtools')
  require('devtools')
}

if (!require('easypackages', character.only=T, quietly=T)) {
  devtools::install_github("jakesherman/easypackages")
  require('easypackages')
  }

library(easypackages)

req_packages <- c('dplyr', 'getPass', 'dbplyr', 'lubridate', 'FNN', 'stringr', 'gbm', 'dismo', 'doParallel', 'RPostgreSQL', 'pool', 'postGIStools')

for (p in req_packages) {
  if(!require(p, character.only=T, quietly=T)){
    install.packages(p)
    print(p)
  }
}
libraries(req_packages)

source('./R/f_calc_distance.R')
source('./R/git_predict_introduction_events.R')
# Run git_access_connectivity_matrix2.R
source('./R/git_access_connectivity_matrix2.R')

workdir <- opt$workdir
focal_group <- opt$species_group

Nsims <- opt$nsims # number of iterations
slope_measure <- opt$slope_measure
focal_species_slope <- as.integer(str_split(opt$slope_thresholds, ",")[[1]])
scenario_name <- opt$scenario_name
if("exterminate_waterbodyid" %in% names(opt)) {
  exwaterbodyID <- as.integer(str_split(opt$exterminate_waterbodyid, ",")[[1]])
} else {
  exwaterbodyID <- NA
}
if("fish_barrier" %in% names(opt)) {
  fish_barrier  <- as.integer(str_split(opt$fish_barrier, ",")[[1]])
} else {
  fish_barrier <- NA
}


if (slope_measure & focal_species_slope) {
  sec_disp <- TRUE
  use_slope_barrier <- TRUE
} else{
  sec_disp <- FALSE
  use_slope_barrier <- FALSE
}

# secondary dispersal stuff
with_secondary <- sec_disp # should secondary spread be included in simulations?
percentage_exter <- opt$percentage_extermination # Give the percentage of focal species populations one wants to exterminate before simulation

temp_inc <- opt$temperature_increase # temperature increase

# Set on or the other to true or false (for upstream dispersal). Probability is based on analyses from Sam and Stefan
use_disp_probability <- FALSE # Currently not supported

# Load wrangled data
load(paste0(workdir,"_",focal_group,".RData"))


time_slot_length <- time_slot_length # Duration of the time-slot, in years (time-span of each time-slot), used to estimate establishment probability
n_time_slots <- 1 #as.integer(sim_duration/time_slot_length)
start_year_sim <- end_year
end_year_sim <- start_year_sim + time_slot_length

pg_drv <- RPostgreSQL::PostgreSQL()
pg_db <- 'nofa'
### Setup database connection
if (file.exists("~/.pgpass") & "host" %in% names(opt) & "user" %in% names(opt)) {
  pgpass <- readLines("~/.pgpass")
  pg_host <-opt$host
  pg_db <- opt$dbname
  pg_user <- opt$user
  pg_password <- strsplit(pgpass[grepl(paste(opt$host, opt$port, opt$dbname, opt$user, sep=":"), pgpass)], ":")[[1]][5]
} else {
#Set connection parameters
host_msg <- 'Please enter host name:'
user_msg <- 'Please enter your user name:'
pw_msg <- "Please enter your password:"
if (rstudioapi::isAvailable()) {
 pg_host <- rstudioapi::showPrompt(title='Hostname', message=host_msg, default='')
 pg_user <- rstudioapi::showPrompt(title='Username', message=user_msg, default='')
 pg_password <- rstudioapi::askForPassword(prompt=pw_msg)
} else {
 pg_host <- getPass(msg=host_msg)
 pg_user <- getPass(msg=user_msg)
 pg_password <- getPass(msg=pw_msg)
}
}
pool <- dbPool(
  drv = pg_drv,
  dbname = pg_db,
  host = pg_host,
  user = pg_user,
  password = pg_password,
  idleTimeout = 36000000
)

con <- poolCheckout(pool)


# git_wrangle.R needs a code review with regards to NULL handling in dates
source('./R/git_wrangle.R')
#add species to new loc_env dataframe based on outdata_data_gjedde
source('./R/get_lake_environment.R')
#add recalculated closest distance based on new data
source('./R/f_calc_distance.R')
source('./R/git_predict_introduction_events.R')
# Run git_access_connectivity_matrix2.R
source('./R/git_access_connectivity_matrix2.R')


brt_model <- paste0("brt_mod_norway_", project_schema, "_", focal_group ,".rds")

# Get model from git_make_brt.R
# source('./R/git_make_brt.R')
brt_mod <- readRDS(paste0(simdir, "/brt_mod_norway_",focal_species_str,".rds"))
covariates <- brt_mod$var.names

outdata <- inndata_timeslot[inndata_timeslot$countryCode =="NO",]  #lake_env[lake_env$countryCode =="NO",]
outdata$countryCode <- factor(outdata$countryCode)
#outdata <- outdata %>% filter(!(county %in% c("Finnmark","Troms","Nordland")))
outdata$county <- factor(outdata$county)
analyse.df <- as.data.frame(outdata)
if(focal_speciesid==26436){
  analyse.df <- subset(analyse.df, select = -c(dist_to_closest_pop_log))
  analyse.df <- merge(analyse.df, lake_env[,c("waterBodyID", "dist_to_closest_pop_log")])
}

#.........................................................................................
# Define input variables ----
#.........................................................................................

# Select inndata geographic range based upon connectivity matrix extent
#inndata <- loc_env[loc_env$waterBodyID %in% wbid_wrid$waterBodyID,]

tmpout <- list()

# Before each simulation run!!!!
print(paste(scenario_name, sec_disp))

# Create new dataframe / vars for simulation bookkeeping.
# Use latest time-slot (if multiple) in inndata
inndata_sim1 <- lake_env#[inndata$t_slot==unique(inndata$t_slot)[1],]

#inndata_sim1<-as.data.frame(inndata_sim1)

# Exterminate percentage of present populations of focal species.
if(with_secondary) {
  barriers <- focal_species_slope
} else {
  barriers <- FALSE
}

if(!is.na(exwaterbodyID)) {
  if(exwaterbodyID == "all") {
    exwaterbodyID <- inndata_sim1$waterBodyID[inndata_sim1[focal_species_str]==1]
  } else if(exwaterbodyID == "random") {
    # random extermination
    # if (percentage_exter > 0 & percentage_exter < 1.0) {
    # inndata_sim1[ sample( which(inndata_sim1[focal_species_str]==1), round(percentage_exter*length(which(inndata_sim1[focal_species_str]==1)))), ][focal_species_str] <- 0
    # } else if (percentage_exter == 1.0) {
    # # Exterminate all pressent populations in VFO Trondelag....
    # inndata_sim1[focal_species_str] <- 0
    #}
  }

  inndata_sim1[,focal_species_str][inndata_sim1$waterBodyID %in% exwaterbodyID] <- 0

  if(sum(inndata_sim1[focal_species_str])==0) { # e.g. solabbor
    inndata_sim1$n_pop <- 0
    inndata_sim1$dist_to_closest_pop_log <- Inf
  } else {
    con <- poolCheckout(pool)
    # recalculate n_pop
    geoselect_no_occ_pop_5000 <-  dbGetQuery(con, paste0('SELECT al.id AS "waterBodyID", count(ol.geom) FROM
                                              (SELECT id, geom FROM nofa.lake WHERE ARRAY[\'NO\'::varchar] <@ "countryCodes") AS al,
                                       (SELECT geom FROM nofa.lake WHERE id IN (SELECT "waterBodyID" FROM ', occurrence_view, ') AND id NOT IN (',toString(unique(exwaterbodyID)),')) AS ol
                                       WHERE ST_DWithin(al.geom, ol.geom, 5000)
                                       GROUP BY al.id'))
    
    inndata_sim1 <- subset(inndata_sim1, select = -c(n_pop))
    inndata_sim1 <- merge(inndata_sim1, geoselect_no_occ_pop_5000, all.x=TRUE)
    inndata_sim1["n_pop"][is.na(inndata_sim1["n_pop"])] <- 0

    # recalculate dist_to_closest_pop_log
    a <- f_calc_dist(outdata=inndata_sim1,species=focal_species)
    inndata_sim1$dist_to_closest_pop_log <- log(a$dist_to_closest_pop)
    poolReturn(con)
  }
}

  if(!is.na(fish_barrier)){
    con <- poolCheckout(pool)
    fish_barrier_loc <- fish_barrier
    fish_barrier <- as.integer(strsplit(dbGetQuery(con, paste0('SELECT upstream_lakes FROM ', conmat_schema, '.', conmat_summary_table, ' WHERE "lakeID" = ', fish_barrier_loc, ';'))$upstream_lakes, ',')[[1]])
    if(scenario[["below"]]){
      fish_barrier <- append(fish_barrier, fish_barrier_loc)
    } else {

    }
  }
for(slope_barrier in barriers) { # max slope for migration upstream

  # j simulation runs...
  for(j in 1:Nsims){

    inndata_sim <- as.data.frame(inndata_sim1) # reset dataset for each simulation run

    # simulate across time periods j
    # to run in parallelization, see: https://www.r-bloggers.com/parallel-r-loops-for-windows-and-linux/
    for(i in 1:n_time_slots){

      ### i.1 predict translocations and store new introductions in temp object
      tmp_trans <- f_predict_introduction_events_gmb(inndata_sim,brt_mod,analyse.df,focal_species_str,temp_inc,start_year_sim, covariates)
      tmp_trans <- tmp_trans[!is.na(tmp_trans[focal_species_str]),]
      # include secondary dispeersal?
      if(with_secondary==TRUE && length(tmp_trans[tmp_trans$introduced==1,]$waterBodyID)>0){

        con <- poolCheckout(pool)

        # get wbID from introductions in run i
        introduction_lakes <- tmp_trans[tmp_trans$introduced==1,]$waterBodyID
        species_lakes <- tmp_trans$waterBodyID[tmp_trans[focal_species_str]==1]
        introduction_wrid <- wbid_wrid$wrid[wbid_wrid$waterBodyID %in% species_lakes]
        species_wrid_wbid <- wbid_wrid[wbid_wrid$waterBodyID %in% species_lakes,]
        disperse_input <- species_wrid_wbid %>%
             group_by(wrid) %>%
             summarise(waterBodyIDs = toString(waterBodyID))

        introduction_lakes <- wbid_wrid$waterBodyID[wbid_wrid$waterBodyID %in% species_lakes]
        #assign new introductions
        tmp_trans$introduced <- ifelse(tmp_trans$waterBodyID %in% introduction_lakes,1,tmp_trans$introduced)

        # select out downstream lakes that does not have species at start of time-slot
        if(use_slope_barrier==TRUE){
          # wbid_wrid_array <- get_wbid_wrid_array(con, species_lakes)
          reachable_lakes_species <- get_reachable_lakes(con, disperse_input$wrid, disperse_input$waterBodyIDs, "slope_7_maximum", slope_barrier, conmat_schema, conmat_table)
          if(length(reachable_lakes_species)>0){
            reachable_lakes_species <- reachable_lakes_species[,1]
            if(length(fish_barrier) > 0) {
            if(length(intersect(fish_barrier, tmp_trans[tmp_trans$introduced==1,]$waterBodyID))==0) {
                reachable_lakes_species <- setdiff(reachable_lakes_species, fish_barrier)
                }
            }
            #downstream_lakes <- get_downstream_lakes(con, unique(introduction_lakes), unique(introduction_wrid))
            # select out downstream lakes that does not have species at start of time-slot
            #downstream_lakes <- downstream_lakes[!(downstream_lakes$downstream_lakes %in% inndata_sim$waterBodyID[inndata_sim[focal_species_str]==1]),]
            reachable_lakes_species <- setdiff(reachable_lakes_species, inndata_sim$waterBodyID[inndata_sim[focal_species_str]==1])
            # finally assign introduction to downstream lakes (without previous obs/intro)

            tmp_trans$introduced <- ifelse(tmp_trans$waterBodyID %in% reachable_lakes_species,1,tmp_trans$introduced)
          }
        }
        ##..or dispersal probability based on analyses from Sam and Stefan
        if(use_disp_probability==TRUE){
          lakes_reachable <- get_reachable_lakes_species(con,unique(introduction_lakes))
          # select out upstream lakes that does not have species at start of time-slot
          lakes_reachable <- lakes_reachable[!(lakes_reachable$accessible_lake %in% inndata_sim$waterBodyID[inndata_sim[focal_species_str]==1]),]
          # add lake intros to introduced vector, based on probability
          tmp_trans$introduced <- ifelse(tmp_trans$waterBodyID %in% lakes_reachable$accessible_lake,rbinom(length(tmp_trans$waterBodyID), size = 1, prob=lakes_reachable$likelihood),tmp_trans$introduced)
        }

      poolReturn(con)

      } # end of secondary==TRUE if statement

      ### Store output from time-period i, simulationrun j

      # create variables to store
      intro <- tmp_trans$waterBodyID[tmp_trans$introduced==1]
      intro_is_secondary <- 0

      # on first sim-run, create new object to store output
      if(i==1 & j==1){
        sim_output <- data.frame()
      }

      if(with_secondary==TRUE){
      intro_is_secondary <- ifelse(intro %in% reachable_lakes_species,TRUE,FALSE)
      }
      if(length(tmp_trans[tmp_trans$introduced==1,]$waterBodyID)>0){

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

      # append sim_output with data from time-slot i, sim j
      sim_output <- bind_rows(sim_output,tmp_output)

      ### modify inndata_sim with new introductions to provide innput to time_slot i+1
      # Add new introductions to "focal_species_str" column in inndata_sim
      tmp1 <- inndata_sim[focal_species_str]
      inndata_sim[focal_species_str][inndata_sim$waterBodyID %in% tmp_trans$waterBodyID[tmp_trans$introduced==1],] <- 1


      # n_pop should be recalculated as well!!!
      #system.time(geoselect_no_species_pop_5000 <- dbGetQuery(con, paste0('SELECT al.id AS "waterBodyID", count(ol.geom) FROM
      #                                       (SELECT id, geom FROM nofa.lake WHERE ecco_biwa_wr IN (',toString(unique(wbid_wrid$wrid)),')) AS al,
      #                                       (SELECT geom FROM nofa.lake WHERE id IN (', toString(unique(intro)), ')) AS ol
      #                                       WHERE ST_DWithin(al.geom, ol.geom, 5000) GROUP BY al.id'))
      #)

      # recalculate distance to closest population and replace values
      # in inndata_sim where distance is smaller than previous;
      # i.e. accounting for situations where closest population is outside
      # the geographic area being investigated
      # NB! need a try function here to make run in cases where there are non intros in simrun
      # and f_calc_dist fail
      tmp <- try(f_calc_dist(outdata=inndata_sim,species=focal_species),silent = TRUE)
      if(!is.null(dim(tmp))){
        inndata_sim$dist_to_closest_pop_log <- ifelse(log(tmp$dist_to_closest_pop)<inndata_sim$dist_to_closest_pop,
                                                  log(tmp$dist_to_closest_pop),
                                                  inndata_sim$dist_to_closest_pop_log)
      } else {
        inndata_sim$dist_to_closest_pop_log <-  inndata_sim$dist_to_closest_pop_log
      }

      tmp <- data.frame()

      # and log variable for model predictions
      # inndata_sim$dist_to_closest_pop_log <- log(inndata_sim$dist_to_closest_pop)
      }
    } # end of i loop

    # display progress
    print(paste("sim-run: ",j,", out of: ",Nsims))

  } # end of j loop

  #............................................................................
  # sum up sim_output pr lake and save / write to database
  #............................................................................
  # store inndata1, raw and lake aggregated sim_output as list in rds object
  source('./R/f_postsim_processing.R')
  sim_output_lake <- f_sim_output_lake(sim_output,inndata_sim1,Nsims)
  sim_output_lake$slope_threshold <- slope_barrier


  tmpout[[scenario_name]] <- list()
  tmpout[[scenario_name]][["sim_output"]] <- sim_output
  tmpout[[scenario_name]][["inndata_sim1"]] <- inndata_sim1
  tmpout[[scenario_name]][["sim_output_lake"]] <- sim_output_lake
  tmpout[[scenario_name]][["focal_species"]] <- focal_species
  tmpout[[scenario_name]][["time_slot_length"]] <- time_slot_length
  tmpout[[scenario_name]][["start_year"]] <- start_year
  tmpout[[scenario_name]][["secondary_dispersal"]] <- sec_disp

  # Write lake-specific summary to database
  dataToWrite <- sim_output_lake

  slopeout <- list()
  slopeout[[paste0("sim_output_", slope_barrier)]] <- sim_output_lake

  con <- poolCheckout(pool)
  dbWriteTable(con, c(project_schema,
                      paste("sim",
                             tolower(focal_species_str),
                             scenario_name,
                             Nsims,
                             "simu",
                             ifelse(sec_disp, paste0(slope_barrier, "_deg"),
                                    "no_second")
                             , sep="_")
                      ), value=dataToWrite,overwrite=TRUE)
  poolReturn(con)
  
}


# Write output to local disk
url <- paste(workdir, "sim_", tolower(focal_species_str), "_",
             format(Sys.time(), "%Y_%m_%d"), ".rds", sep="")
saveRDS(tmpout, url)

# Rydde opp i forbindelse mot server
poolReturn(con)
poolClose(pool)


# Check effect of distance
# plot(brt_mod, i.var=2)
# brt_mod_o <- readRDS(paste0(workdir, brt_models[[as.character(26138)]]))
# brt_mod_s <- readRDS(paste0(workdir, brt_models[[as.character(26436)]]))
# plot(brt_mod_s, i.var=2)
# plot(brt_mod_o, i.var=2)
