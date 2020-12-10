library(postGIStools)
library(RPostgreSQL)
#Set connection parameters
#pg_drv <- dbDriver("PostgreSQL")
#user_msg <- 'Please enter your user:'
#pw_msg <- "Please enter your password:"
#if (any(grepl("RStudio", .libPaths()))) {
#  pg_user <- rstudioapi::askForPassword(user_msg)
#  pg_password <- rstudioapi::askForPassword(pw_msg)
#} else {
#  pg_user <- readline(prompt=user_msg)
#  pg_password <- readline(prompt=pw_msg)
#}
#pg_host <- "vm-srv-wallace.vm.ntnu.no"
#pg_db<-"nofa"

#Initialise connection
#con<-dbConnect(pg_drv,dbname=pg_db,user=pg_user, password=pg_password,host=pg_host)

get_historic_distribution <- function(connection, taxonID) {
  species_native <- get_postgis_query(con, paste0('SELECT
                                                 ST_Transform(geom, 4326) AS geom, occurrence
                                                 FROM nofa."l_taxon_historicDistribution"
                                                 WHERE occurrence > 1 AND "taxonID" = ', taxonID),  geom_name = "geom")
  species_native
}
