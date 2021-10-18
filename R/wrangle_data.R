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
# Rscript --vanilla ~/fremmedfisk/R/wrangle_data.R -s Agn -H "vm-srv-wallace.vm.ntnu.no" -u stefan -D ~/fremmedfisk/ -c invafish -P fremmedfisk -w /data/R/Prosjekter/13845000_invafish/fremmedfisk_opsjon/simulering/ -y 1967 -l 50

# Run script on the command lin in parallel for all species groups:
# species_groups=$(psql -h "vm-srv-wallace.vm.ntnu.no" -d nofa -U "stefan" -At -c "SELECT DISTINCT ON (species_group) species_group FROM fremmedfisk.species_groups;")
# echo "$species_groups" | awk '{print("Rscript --vanilla ~/fremmedfisk/R/wrangle_data.R -s " $1 " -H \"vm-srv-wallace.vm.ntnu.no\" -u stefan -D ~/fremmedfisk/ -c invafish -P fremmedfisk -w /data/R/Prosjekter/13845000_invafish/fremmedfisk_opsjon/simulering/ -y 1967 -l 50 \0")}' | xargs -0 -P 5 -I{} bash -c "{}"

# 
# To run in IDE/RStudio skip option parsing and use this opt list
# opt <- list("species_group"="Eksotisk", "dbname"="nofa", "user"="stefan", "port"=5432, "host"="vm-srv-wallace.vm.ntnu.no", "scriptdir"='~/fremmedfisk/', "connectivity_schema"="invafish", "project_schema"="fremmedfisk", "workdir"="/data/R/Prosjekter/13845000_invafish/fremmedfisk_opsjon/simulering/", "start_year"=1967, "time_slot_length"=50)

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

save.image(paste0(workdir,"_",focal_group,".RData"))

# Rydde opp i forbindelse mot server
poolReturn(con)
poolClose(pool)
