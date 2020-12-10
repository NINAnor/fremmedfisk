# library(RPostgreSQL)
# library(pool)

# # Set connection parameters
# pg_drv <- RPostgreSQL::PostgreSQL()
# pg_host <- "vm-srv-finstad.vm.ntnu.no"
# pg_db <- 'nofa'
# pg_user <- rstudioapi::askForPassword("enter username")
# pg_password <- rstudioapi::askForPassword("enter psw")

# pool <- dbPool(
#   drv = pg_drv,
#   dbname = pg_db,
#   host = pg_host,
#   user = pg_user,
#   password = pg_password,
#   idleTimeout = 36000000
# )

# con <- poolCheckout(pool)

####################################################################################
# N?

### Hent ut milj?data for alle innsj?er i Agder
get_location_environment <- function(db_conection) {
  sql_string <- "SELECT * FROM nofa.view_location_environment AS a WHERE county IN ('Vest-Agder', 'Aust-Agder');"
  res <- dbGetQuery(db_conection, sql_string)
  res
}
# Eksempel
#loc_env <- get_location_environment(con)
