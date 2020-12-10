library(RPostgreSQL)
library(pool)

#Set connection parameters
#pg_drv <- RPostgreSQL::PostgreSQL()
#pg_host <- "vm-srv-wallace.vm.ntnu.no"
#pg_db <- 'nofa'
#user_msg <- 'Please enter your user:'
#pw_msg <- "Please enter your password:"
#if (any(grepl("RStudio", .libPaths()))) {
#  pg_user <- rstudioapi::askForPassword(user_msg)
#  pg_password <- rstudioapi::askForPassword(pw_msg)
#} else {
#  pg_user <- readline(prompt=user_msg)
#  pg_password <- readline(prompt=pw_msg)
#}

#pool <- dbPool(
#  drv = pg_drv,
#  dbname = pg_db,
#  host = pg_host,
#  user = pg_user,
#  password = pg_password,
#  idleTimeout = 36000000
#)

#con <- poolCheckout(pool)

####################################################################################
# N?

### Hent ut milj?data for alle innsj?er i Agder
get_lake_environment <- function(db_conection, waterBodyIDs) {
  sql_string <- paste("SELECT * FROM nofa.view_lake_environment AS a WHERE \"waterBodyID\" IN (", toString(waterBodyIDs, sep=','), ");", sep='')
  res <- dbGetQuery(db_conection, sql_string)
  res
}
# Eksempel
# lake_env <- get_lake_environment(con, wbid_wrid$waterBodyID)
