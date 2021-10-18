#!/usr/bin/env Rscript


# To clarify
# - "start populations" kun introduserte eller alle, hva med ukjente
# - Hvordan håndtere forekomst av flere arter fra en gruppe i et vann
# - Hvordan håndtere forekomst av flere arter fra en gruppe i et vann
#   Relevant i f_predict_introduction_events_gmb
# - Hvordan håndtere sampling bias i rom og tid (begge er ofte skev fordelt)
#   Relevant i f_predict_introduction_events_gmb og kanskje også bygging av gbm
# Parallelize simulations (they are slow (both predictions and with many lakes)


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

scriptdir <- '~/INVAFISH-sim/'
setwd(scriptdir)

slope_thresholds <- c(600, 700, 800)

brt_dir <- "/data/R/Prosjekter/13845000_invafish/"
simdir <- paste0(brt_dir, 'temp/')
try(system(paste0('mkdir ', simdir)))

focal_group <- "Agn"

start_year <- 1967
end_year <- 2017

conmat_schema <- "invafish"
conmat_table <- "invafish_lake_connectivity"
conmat_summary_table <- "invafish_lake_connectivity_summary"

project_schema <- "fremmedfisk"

slope_measure <- "slope_7_maximum"

### Setup database connection
#Set connection parameters
pg_drv <- RPostgreSQL::PostgreSQL()
pg_db <- 'nofa'
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

pool <- dbPool(
  drv = pg_drv,
  dbname = pg_db,
  host = pg_host,
  user = pg_user,
  password = pg_password,
  idleTimeout = 36000000
)

con <- poolCheckout(pool)


# Run git_access_connectivity_matrix2.R
source('./R/git_access_connectivity_matrix2.R')
# git_wrangle.R needs a code review with regards to NULL handling in dates
source('./R/git_wrangle.R')
#add species to new loc_env dataframe based on outdata_data_gjedde
source('./R/get_lake_environment.R')
#add recalculated closest distance based on new data
source('./R/f_calc_distance.R')
source('./R/git_predict_introduction_events.R')


brt_models <- list("Agn"="brt_mod_norway_fremmedfisk_Agn.rds", # Agn
                   "26436"="brt_mod_norway_fremmedfisk_Lepomis gibbosus.rds", # solabbor
                   "26138"="brt_mod_norway_fremmedfisk_Phoxinus phoxinus.rds"  # ørekyt
)

spec_scenarios <- list("Agn"=list("scenario_0"=list("exwaterbodyID"=NA, fish_barrier=NA, below=NA),
                                    "scenario_1"=list("exwaterbodyID"=4118447, fish_barrier=NA, below=NA),
                                    "scenario_2"=list("exwaterbodyID"=c(4118446, 4118447), fish_barrier=NA, below=NA)
), # Agn
"26436"=list("scenario_0"=list("exwaterbodyID"=NA, fish_barrier=NA, below=TRUE),
             "scenario_1"=list("exwaterbodyID"="all", fish_barrier=NA, below=TRUE)
), # solabbor
"26138"=list("scenario_0"=list("exwaterbodyID"=NA, fish_barrier=NA, below=NA),
             "scenario_1"=list("exwaterbodyID"=4110252, fish_barrier=4110252, below=TRUE),
             "scenario_2"=list("exwaterbodyID"=NA, fish_barrier=4110252, below=FALSE)
) # ørekyt
)



focal_species_group <- dbGetQuery(con, paste0('SELECT "species_group", string_agg(CAST(taxonid AS varchar), \',\') AS taxonid,  string_agg(CAST("scientificName" AS varchar), \',\') AS scientificname FROM (SELECT * FROM fremmedfisk.species_groups NATURAL INNER JOIN (SELECT "taxonID" AS taxonid, "scientificName" FROM nofa.l_taxon) AS t WHERE species_group = \'',focal_group,'\') AS x GROUP BY "species_group";')) # WHERE "taxonID" = ', focal_group, ';'))[,1]
  #dbGetQuery(con, paste0('SELECT "taxonID", scientificName" FROM nofa.l_taxon WHERE "taxonID" = ', focal_group, ';'))[,1]
focal_species_str <- tolower(focal_group)
focal_species_slope <- slope_thresholds
focal_species_scenarios <- spec_scenarios[[as.character(focal_group)]]
focal_species <- strsplit(focal_species_group[,3], ",")[[1]]

project_schema <- "fremmedfisk"
occurrence_view_name <- paste0(focal_species_str, "_occurrence")
occurrence_view <- paste0('"', project_schema, '"."', occurrence_view_name, '"')


dbGetQuery(con, paste0('DROP VIEW IF EXISTS ', occurrence_view, ';
CREATE VIEW ', occurrence_view, ' AS SELECT 
  l.*,
  e."dateStart",
  e."dateEnd",
  e.modified,
  o.*
    FROM (SELECT 
        "taxonID",
        "occurrenceStatus", 
        "establishmentMeans", 
        "eventID"
        FROM nofa.occurrence
        WHERE "taxonID" IN (', gsub(",", ", ", focal_species_group[,2]), ') AND "occurrenceStatus" = \'present\') AS o
NATURAL LEFT JOIN nofa.event AS e
LEFT JOIN (SELECT "waterBodyID", "locationID" FROM nofa.location WHERE "countryCode" = \'NO\') AS l USING ("locationID")
WHERE "locationID" IS NOT NULL AND "waterBodyID" IS NOT NULL;'))


### Hent ut alle vannregioner fra innsjødatasettet
wbid_wrid <- dbGetQuery(con, 'SELECT ecco_biwa_wr AS wrid, id AS "waterBodyID" FROM nofa.lake WHERE ecco_biwa_wr IS NOT NULL;') # get_wbid_wrid(con)

# Get occurrences
inndata <- tbl(con, sql("SELECT * FROM nofa.view_occurrence_by_event")) %>% collect()

# From get_geoselect_native.R - to be replaced by quality check from Trygve
# source('./R/get_geoselect_native.R')
#if(focal_speciesid != 26436) {
#geoselect_native <- get_historic_distribution(con, focal_speciesid)
#} else {
  geoselect_native <- SpatialPolygonsDataFrame(SpatialPolygons(list(), proj4string = CRS(as.character("+proj=longlat +datum=WGS84 +no_defs"))), data=data.frame())
#}

# Should we keep county specific data handling? maybe better to drop it...
counties <- dbGetQuery(con, 'SELECT DISTINCT (county) county FROM "AdministrativeUnits"."Fenoscandia_Municipality_polygon" WHERE "countryCode" = \'NO\';')  # AND county NOT IN (\'Finnmark\',\'Troms\',\'Nordland\');')

geoselect_no_species_pop_5000 <- dbGetQuery(con, paste0('SELECT al.id AS "waterBodyID", count(ol.geom) FROM
                                                  (SELECT id, geom FROM nofa.lake WHERE ARRAY[\'NO\'::varchar] <@ "countryCodes") AS al,
                                           (SELECT geom FROM nofa.lake WHERE id IN (SELECT "waterBodyID" FROM ', occurrence_view, ')) AS ol
                                           WHERE ST_DWithin(al.geom, ol.geom, 5000)
                                           GROUP BY al.id'))
names(geoselect_no_species_pop_5000)[2] = "n_pop"

inndata <- merge(inndata, geoselect_no_species_pop_5000, all.x=TRUE)
inndata["n_pop"][is.na(inndata["n_pop"])] <- 0

# From git_wrangle.R
outdata_list <- wrangle_and_slice(start_year,end_year,inndata,focal_species,focal_species_str,geoselect_native)

inndata_timeslot <- outdata_list$data

# Get lake environment data from get_lake_environment.R
#add species to new loc_env dataframe based on outdata_data_gjedde
lake_env <- get_lake_environment(con, wbid_wrid$waterBodyID)

for(s in focal_species) {
  lake_env[focal_species_str] <- inndata_timeslot[gsub(" ", "_", s)][match(as.numeric(lake_env$waterBodyID), inndata_timeslot$waterBodyID),]
}
lake_env[focal_species_str][is.na(lake_env[focal_species_str])] <- 0

lake_env$introduced <- inndata_timeslot$introduced[match(as.numeric(lake_env$waterBodyID), inndata_timeslot$waterBodyID)]
lake_env$introduced[is.na(lake_env$introduced)] <- 0

lake_env <- merge(lake_env, geoselect_no_species_pop_5000, by="waterBodyID", all.x=TRUE)
lake_env$n_pop <- ifelse(is.na(lake_env$n_pop), 0, lake_env$n_pop)

#add recalculated closest distance based on new data
a <- f_calc_dist(outdata=lake_env,species=focal_species_str)
lake_env$dist_to_closest_pop_log <- log(a$dist_to_closest_pop)






# Get model from git_make_brt.R
# source('./R/git_make_brt.R')
brt_mod <- readRDS(paste0(brt_dir, brt_models[[as.character(focal_speciesid)]]))
covariates <- brt_mod$var.names

outdata <-inndata_timeslot[inndata_timeslot$countryCode =="NO",]  #lake_env[lake_env$countryCode =="NO",]
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

# simulation and time specific stuff
Nsims <- 100 # number of iterations
sim_duration <- 1 # Duration of the scenario, in years (will be corced to multiple of time_slot_length)
time_slot_length <- 50 # Duration of the time-slot, in years (time-span of each time-slot)
gmb_time_slot_length <- 50 # Duration of the time-slot, in years, used to estimate establishment probability
n_time_slots <- 1#as.integer(sim_duration/time_slot_length)
start_year_sim <- 2017
end_year_sim <- start_year_sim + time_slot_length

tmpout <- list()


for(sec_disp in TRUE) { #c(FALSE, TRUE)
  temp_inc <- 0 # temperature increase

  # secondary dispersal stuff
  with_secondary <- sec_disp # should secondary spread be included in simulations?
  percentage_exter <- 0.0 # Give the percentage of focal species populations one wants to exterminate before simulation

  # Set on or the other to true or false (for upstream dispersal). Probability is based on analyses from Sam and Stefan
  use_slope_barrier <- TRUE
  use_disp_probability <- FALSE

  # Before each simulation run!!!!
  for(scenario_name in names(focal_species_scenarios)) {

    print(paste(scenario_name, sec_disp))

    scenario <- focal_species_scenarios[[scenario_name]]
    exwaterbodyID <- scenario[["exwaterbodyID"]]
    fish_barrier  <- scenario[["fish_barrier"]]

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
  }
}

# Write output to local disk
url <- paste(brt_dir, "sim_", tolower(focal_species_str), "_",
             format(Sys.time(), "%Y_%m_%d"), ".rds", sep="")
saveRDS(tmpout, url)

# Rydde opp i forbindelse mot server
poolReturn(con)
poolClose(pool)


# Check effect of distance
# plot(brt_mod, i.var=2)
# brt_mod_o <- readRDS(paste0(brt_dir, brt_models[[as.character(26138)]]))
# brt_mod_s <- readRDS(paste0(brt_dir, brt_models[[as.character(26436)]]))
# plot(brt_mod_s, i.var=2)
# plot(brt_mod_o, i.var=2)
