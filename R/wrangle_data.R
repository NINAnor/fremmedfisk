#!/usr/bin/env Rscript


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
  make_option(c("-D", "--scriptdir"), type="character", default="./", 
              help="Name of the database to use", metavar="character"),
  make_option(c("-d", "--dbname"), type="character", default="nofa", 
              help="Name of the database to use", metavar="character"),
  make_option(c("-p", "--port"), type="integer", default=5432, 
              help="Name of the database to use", metavar="integer"),
  make_option(c("-u", "--user"), type="character", default=NULL, 
              help="Name of the user to connect to the database", metavar="character"),
  make_option(c("-H", "--host"), type="character", default=NULL, 
              help="Name of the host where the database runs", metavar="character"),
  make_option(c("-c", "--connectivity_schema"), type="character", default="invafish", 
              help="Name of the database schema where the connectivity matrix lives", metavar="character"),
  make_option(c("-P", "--project_schema"), type="character", default="fremmedfisk", 
              help="Name of the database schema where the project data is supposed to be stored", metavar="character"),
  make_option(c("-y", "--start_year"), type="integer", default=1967, 
              help="Start year for introduction modeling", metavar="integer"),
  make_option(c("-l", "--time_slot_length"), type="integer", default=50, 
              help="Length of the time slot for introduction modeling", metavar="integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Run with arguments:
# Rscript --vanilla ~/fremmedfisk/R/wrangle_data.R -s Predator -H "vm-srv-wallace.vm.ntnu.no" -u stefan -D ~/fremmedfisk/ -c invafish -P fremmedfisk -w /data/R/Prosjekter/13845000_invafish/fremmedfisk_opsjon/simulering/ -y 1967 -l 50
#
# Run script on the command line in parallel for all species groups:
# species_groups=$(psql -h "vm-srv-wallace.vm.ntnu.no" -d nofa -U "stefan" -At -c "SELECT DISTINCT ON (species_group) species_group FROM fremmedfisk.species_groups;")
# echo "$species_groups" | awk '{print("Rscript --vanilla ~/fremmedfisk/R/wrangle_data.R -s " $1 " -H \"vm-srv-wallace.vm.ntnu.no\" -u stefan -D ~/fremmedfisk/ -c invafish -P fremmedfisk -w /data/R/Prosjekter/13845000_invafish/fremmedfisk_opsjon/simulering/ -y 1967 -l 50 \0")}' | xargs -0 -P 5 -I{} bash -c "{}"

# 
# To run in IDE/RStudio skip option parsing and use this opt list
# opt <- list("species_group"="Predator", "dbname"="nofa", "user"="stefan", "port"=5432, "host"="vm-srv-wallace.vm.ntnu.no", "scriptdir"='~/fremmedfisk/', "connectivity_schema"="invafish", "project_schema"="fremmedfisk", "workdir"="/data/R/Prosjekter/13845000_invafish/fremmedfisk_opsjon/simulering/", "start_year"=1967, "time_slot_length"=50)

# To clarify
# - "start populations" kun introduserte eller alle, hva med ukjente
# - Hvordan håndtere forekomst av flere arter fra en gruppe i et vann
#   Relevant i f_predict_introduction_events_gmb
# - Hvordan håndtere sampling bias i rom og tid (begge er ofte skev fordelt)
#   Relevant i f_predict_introduction_events_gmb og kanskje også bygging av gbm
# - Vurdere effekt som utenfor naturlig utbredelse (vanskelig 7 for å ikke si umulig
#   når artsgrupper med ulik naturlig utbredelse analyseres)

# Neste steg:
# Aggregering på vannregion / bi-kuber nivå
# Parallelize simulations (they are slow (mainly predictions with many lakes)
# Split script into three parts (1, data collection/wrangling, 2. model building,
#   3. simulation and writing of data) so it can be run on the command line to avoid timeouts


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

# Potential arguments
focal_group <- opt$species_group
scriptdir <- opt$scriptdir
setwd(scriptdir)
workdir <- opt$workdir
start_year <- opt$start_year
time_slot_length <- opt$time_slot_length
end_year <- start_year + time_slot_length
conmat_schema <- opt$connectivity_schema
conmat_table <- paste0(conmat_schema, "_lake_connectivity")
conmat_summary_table <- paste0(conmat_schema, "_lake_connectivity_summary")
project_schema <- opt$project_schema


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

focal_species_group <- dbGetQuery(con, paste0('SELECT "species_group", string_agg(CAST(taxonid AS varchar), \',\') AS taxonid,  string_agg(CAST("scientificName" AS varchar), \',\') AS scientificname FROM (SELECT * FROM fremmedfisk.species_groups NATURAL INNER JOIN (SELECT "taxonID" AS taxonid, "scientificName" FROM nofa.l_taxon) AS t WHERE species_group = \'',focal_group,'\') AS x GROUP BY "species_group";')) # WHERE "taxonID" = ', focal_group, ';'))[,1]
  #dbGetQuery(con, paste0('SELECT "taxonID", scientificName" FROM nofa.l_taxon WHERE "taxonID" = ', focal_group, ';'))[,1]
focal_species_str <- tolower(focal_group)
focal_species <- strsplit(focal_species_group[,3], ",")[[1]]

occurrence_view_name <- paste0(focal_species_str, "_occurrence")
occurrence_view <- paste0('"', project_schema, '"."', occurrence_view_name, '"')

view_exists <- NA #tryCatch(dbGetQuery(con, paste0('SELECT "waterBodyID" FROM ', occurrence_view, ' LIMIT 1;')), error=function(cond) {return(NA)}, warning=function(cond) {return(NA)})

if(is.na(view_exists)) {
# year 500 is before any occurrences are registered in nofa
dbGetQuery(con, paste0('DROP VIEW IF EXISTS', occurrence_view, ';\n CREATE VIEW ', occurrence_view, ' AS SELECT 
"waterBodyID" FROM (SELECT (nofa.get_occurrence_status_timeslot(
	\'', focal_species_group$taxonid,'\'::text, \'500-01-01\'::timestamp, \'', start_year,'-01-01\'::timestamp
)).*) AS a
WHERE "occurrenceStatus" = \'present\'
ORDER BY "waterBodyID" ASC'))
}

### Hent ut alle vannregioner fra innsjødatasettet
wbid_wrid <- dbGetQuery(con, 'SELECT ecco_biwa_wr AS wrid, id AS "waterBodyID" FROM nofa.lake WHERE ecco_biwa_wr IS NOT NULL AND ARRAY[\'NO\'::varchar] <@ "countryCodes";') # get_wbid_wrid(con)

# Get occurrences
observation_lakes <- tbl(con, sql("SELECT DISTINCT ON (\"waterBodyID\") \"waterBodyID\" FROM nofa.view_occurrence_by_event WHERE \"countryCode\" = 'NO'")) %>% collect()

# From get_geoselect_native.R - to be replaced by quality check from Trygve
# source('./R/get_geoselect_native.R')
#if(focal_speciesid != 26436) {
#geoselect_native <- get_historic_distribution(con, focal_speciesid)
#} else {
#  geoselect_native <- SpatialPolygonsDataFrame(SpatialPolygons(list(), proj4string = CRS(as.character("+proj=longlat +datum=WGS84 +no_defs"))), data=data.frame())
#}

# Should we keep county specific data handling? maybe better to drop it...
#counties <- dbGetQuery(con, 'SELECT DISTINCT (county) county FROM "AdministrativeUnits"."Fenoscandia_Municipality_polygon" WHERE "countryCode" = \'NO\';')  # AND county NOT IN (\'Finnmark\',\'Troms\',\'Nordland\');')

geoselect_no_species_pop_5000 <- dbGetQuery(con, paste0('SELECT al.id AS "waterBodyID", count(ol.geom) FROM
                                                  (SELECT id, geom FROM nofa.lake WHERE ARRAY[\'NO\'::varchar] <@ "countryCodes") AS al,
                                           (SELECT nl.id, nl.geom FROM (SELECT "waterBodyID" AS id FROM ', occurrence_view, ') AS ov NATURAL LEFT JOIN nofa.lake AS nl) AS ol
                                           WHERE ST_DWithin(al.geom, ol.geom, 5000) AND al.id != ol.id
                                           GROUP BY al.id'))
names(geoselect_no_species_pop_5000)[2] = "n_pop"


# Get lake environment data from get_lake_environment.R
# add species to new loc_env dataframe based on outdata_data_gjedde
lake_env <- get_lake_environment(con, wbid_wrid$waterBodyID)
lake_env$decimalLatitude_R<-round(lake_env$decimalLatitude,1)

lake_env <- merge(lake_env, geoselect_no_species_pop_5000, all.x=TRUE)
lake_env["n_pop"][is.na(lake_env["n_pop"])] <- 0


lake_env <- merge(lake_env, data.frame(waterBodyID=observation_lakes, observation_lake=1), all.x=TRUE)
lake_env["observation_lake"][is.na(lake_env["observation_lake"])] <- 0

# From git_wrangle.R
#outdata_list <- wrangle_and_slice(start_year,end_year,inndata,focal_species,focal_species_str,geoselect_native)

#inndata_timeslot <- outdata_list$data

occurrences_start <- dbGetQuery(con, paste0('SELECT  "waterBodyID", "taxonID", "occurrenceStatus", "establishmentMeans" 
FROM (SELECT (nofa.get_occurrence_status_timeslot(
  \'', focal_species_group$taxonid,'\'::text, \'500-01-01\'::timestamp, \'', start_year,'-01-01\'::timestamp, NULL, NULL, \'NO\'
)).*) AS a
-- WHERE "occurrenceStatus" = \'present\''))
names(occurrences_start) <- c("waterBodyID", "taxonID", paste("occurrenceStatus", start_year, sep='_'), paste("establishmentMeans", start_year, sep='_'))

occurrences_end <- dbGetQuery(con, paste0('SELECT "waterBodyID", "taxonID", "occurrenceStatus", "establishmentMeans", "presentWithinTimeslot" 
FROM (SELECT (nofa.get_occurrence_status_timeslot(
  \'', focal_species_group$taxonid,'\'::text, \'', start_year,'-01-01\'::timestamp, \'', end_year-1,'-12-31\'::timestamp, NULL, NULL, \'NO\'
)).*) AS a -- WHERE "occurrenceStatus" = \'present\''))
names(occurrences_end) <- c("waterBodyID", "taxonID", paste("occurrenceStatus", end_year, sep='_'), paste("establishmentMeans", end_year, sep='_'), "presentWithinTimeslot")

occurrences_timeslot <- merge(occurrences_start, occurrences_end, by=c("waterBodyID", "taxonID"), all.x=TRUE, all.y=TRUE)

occurrences_timeslot[paste("occurrenceStatus", start_year, sep='_')] <- as.integer(ifelse(occurrences_timeslot[paste("occurrenceStatus", start_year, sep='_')]=='absent' | is.na(occurrences_timeslot[paste("occurrenceStatus", start_year, sep='_')]), 0, 1))
occurrences_timeslot[paste("occurrenceStatus", end_year, sep='_')] <- as.integer(ifelse(occurrences_timeslot[paste("occurrenceStatus", end_year, sep='_')]=='absent' | is.na(occurrences_timeslot[paste("occurrenceStatus", end_year, sep='_')]), 0, 1))
occurrences_timeslot["presentWithinTimeslot"] <- as.integer(ifelse(is.na(occurrences_timeslot["presentWithinTimeslot"]), 0, 1))
occurrences_timeslot$introduced <- as.integer(ifelse(occurrences_timeslot[paste("occurrenceStatus", start_year, sep='_')] != 1 & (occurrences_timeslot[paste("occurrenceStatus", end_year, sep='_')] == 1 | occurrences_timeslot$presentWithinTimeslot) & occurrences_timeslot[paste("establishmentMeans", end_year, sep='_')] == 'introduced', 1, 0))

lake_env_species <- merge(lake_env, occurrences_timeslot[c("waterBodyID", "taxonID", paste("occurrenceStatus", start_year, sep='_'), paste("occurrenceStatus", end_year, sep='_'), "introduced")], by="waterBodyID", all.x=TRUE, all.y=TRUE)
lake_env_species[paste("occurrenceStatus", start_year, sep='_')][is.na(lake_env_species[paste("occurrenceStatus", start_year, sep='_')])] <- as.integer(0)
lake_env_species[paste("occurrenceStatus", end_year, sep='_')][is.na(lake_env_species[paste("occurrenceStatus", end_year, sep='_')])] <- as.integer(0)
lake_env_species["introduced"][is.na(lake_env_species["introduced"])] <- as.integer(0)

lake_env_species <- lake_env_species[lake_env_species$observation_lake == 1,]

summary_columns <- c(paste("occurrenceStatus", start_year, sep='_'), paste("occurrenceStatus", end_year, sep='_'), "introduced")
occurrences_timeslot_wbid <- as.data.frame(occurrences_timeslot[append("waterBodyID",summary_columns)] %>%
  group_by(waterBodyID) %>%
  summarize_at(.vars = summary_columns, .funs = max))


lake_env <- merge(lake_env, occurrences_timeslot_wbid, by="waterBodyID", all.x=TRUE)
lake_env[paste("occurrenceStatus", start_year, sep='_')][is.na(lake_env[paste("occurrenceStatus", start_year, sep='_')])] <- as.integer(0)
lake_env[paste("occurrenceStatus", end_year, sep='_')][is.na(lake_env[paste("occurrenceStatus", end_year, sep='_')])] <- as.integer(0)
lake_env["introduced"][is.na(lake_env["introduced"])] <- as.integer(0)

#for(s in focal_species) {
#  lake_env[focal_species_str] <- inndata_timeslot[gsub(" ", "_", s)][match(as.numeric(lake_env$waterBodyID), inndata_timeslot$waterBodyID),]
#}
#lake_env[focal_species_str][is.na(lake_env[focal_species_str])] <- 0

#lake_env$introduced <- as.integer(
#  ifelse(lake_env[paste("occurrenceStatus", start_year, sep='_')] == 0 & lake_env[paste("occurrenceStatus", end_year, sep='_')] == 1 & lake_env[paste("establishmentMeans", end_year, sep='_')] == 'introduced', 1, 0))
# inndata_timeslot$introduced[match(as.numeric(lake_env$waterBodyID), inndata_timeslot$waterBodyID)]
# lake_env$introduced[is.na(lake_env$introduced)] <- 0

#lake_env <- merge(lake_env, geoselect_no_species_pop_5000, by="waterBodyID", all.x=TRUE)
#lake_env$n_pop <- ifelse(is.na(lake_env$n_pop), 0, lake_env$n_pop)

#add recalculated closest distance based on new data
a <- f_calc_dist(outdata=lake_env,occurrences=paste("occurrenceStatus", start_year, sep='_'))
a$dist_to_closest_pop_log <- log(a$dist_to_closest_pop)
lake_env <- merge(lake_env, a, by="waterBodyID", all.x=TRUE)
lake_env_species <- merge(lake_env_species, a, by="waterBodyID", all.x=TRUE)
rm(a)
# Remove objects that should not be kept for the following operations
# rm(opt)

lake_env_species$countryCode <- factor(lake_env_species$countryCode)
lake_env_species$county <- factor(lake_env_species$county)

# Rydde opp i forbindelse mot server
poolReturn(con)
poolClose(pool)
rm(pool)

save.image(paste0(workdir,"_",focal_group,".RData"))

