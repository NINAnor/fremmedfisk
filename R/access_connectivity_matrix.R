# Før
# # get vector of wbID for upstream lakes 
# upstream_lakes <- connectivity$upstream_lakes[connectivity$waterBodyID %in% introduction_lakes]      
# upstream_lakes <- paste(upstream_lakes,sep=",",collapse=",")
# upstream_lakes <- as.numeric(unlist(strsplit(upstream_lakes,split=",")))
# # get vector of upstream slopes 
# upstream_slopes <- connectivity$upstream_lakes_slope_max_max[connectivity$waterBodyID %in% introduction_lakes]
# upstream_slopes <- paste(upstream_slopes,sep=",",collapse=",")
# upstream_slopes <- unlist(strsplit(upstream_slopes,split=","))
# upstream_slopes <- gsub(" ","",upstream_slopes)
# upstream_slopes <- as.numeric(upstream_slopes)

# # finally reachable lakes and assign introduction 
# upstream_lakes_reachable <- na.omit(upstream_lakes[upstream_slopes<slope_barrier])
# upstream_lakes_reachable <- unique(upstream_lakes_reachable) # introductions only listed once (remove duplicated lakeIDs)

# SELECT
# *
# FROM 
# temporary_agder_connectivity.lake_connectivity_summary WHERE 


####################################################################################
# Sette opp en test case

slope_barrirer <- 0.01

library(RPostgreSQL)
library(pool)

#Set connection parameters
pg_drv <- RPostgreSQL::PostgreSQL()
pg_host <- "vm-srv-wallace.vm.ntnu.no"
pg_db <- 'nofa'
pg_user <- rstudioapi::askForPassword("enter username")
pg_password <- rstudioapi::askForPassword("enter psw")


pool <- dbPool(
  drv = pg_drv,
  dbname = pg_db,
  host = pg_host,
  user = pg_user,
  password = pg_password,
  idleTimeout = 36000000
)

con <- poolCheckout(pool)

####################################################################################
# Nå

### Hent ut alle vannregioner (for Agder (mange inneholder ingen innsjøer)
get_wrids <- function(db_conection) {
  sql_string <- "SELECT DISTINCT ON (a.gid) a.gid AS wrid FROM \"Hydrography\".\"waterregions_dem_10m_nosefi\" AS a, (SELECT geom FROM \"AdministrativeUnits\".\"Fenoscandia_Municipality_polygon\" WHERE county IN ('Vest-Agder', 'Aust-Agder')) AS b WHERE ST_Intersects(a.geom, b.geom);"
  res <- dbGetQuery(db_conection, sql_string)[,1]
  res
}
# Eksempel
wrids <- get_wrids(con)

### Hent ut data frame med kombinasjon av waterbodyID for alle innsjøer (kolonne 1) og id for vannregioner (kolonne 2) (her er det kun vannregioner som inneholder innsjøer)
get_wbid_wrid <- function(db_conection, eb_waterregionID) {
  sql_string <- paste("SELECT id AS \"waterBodyID\", ecco_biwa_wr AS wrid FROM nofa.lake WHERE ecco_biwa_wr IN (", toString(eb_waterregionID, sep=','), ")", sep='')
  res <- dbGetQuery(db_conection, sql_string)
  res
}
# Eksempel
wbid_wrid <- get_wbid_wrid(con, wrids)

### Hent ut data frame med unike kombinasjon av waterbodyID (kolonne 1) og id for innsjøer som ligger nedstrøms (kolonne 2)
get_downstream_lakes <- function(db_conection, waterbodyID, eb_waterregionID) {
  sql_string <- paste("SELECT \"lakeID\" AS  \"waterBodyID\", CAST(unnest(string_to_array(downstream_lakes, ',')) AS integer) AS downstream_lakes FROM
                      agder.lake_connectivity_summary WHERE 
                      wrid IN (", toString(eb_waterregionID, sep=','), ") AND
                      \"lakeID\" IN (", toString(waterbodyID, sep=','),");", sep='')
  res <- dbGetQuery(db_conection, sql_string)
  res
}
# Eksempel
downstream_lakes <- get_downstream_lakes(con, unique(wbid_wrid[,1][1:100]), unique(wbid_wrid[,2][1:100]))


### Hent ut data frame med kombinasjon av waterbodyID (kolonne 1) og id for innsjøer som ligger oppstøms og der skråning i forbindelsen er lavere enn slope_barrier (kolonne 2)
get_reachable_upstream_lakes <- function(db_conection, waterbodyID, eb_waterregionID, slope_barrier) {
  sql_string <- paste("SET constraint_exclusion = on;
                      SELECT
                      from_lake AS source_lake, to_lake AS upstream_lake
                      FROM
                     agder.lake_connectivity
                      WHERE
                      wrid in (", toString(eb_waterregionID, sep=','), ") AND
                      from_lake in (", toString(waterbodyID, sep=','), ") AND
                      upstream_slope_max_max <= ", slope_barrier, "
                      UNION ALL SELECT
                      from_lake AS source_lake, to_lake AS upstream_lake
                      FROM
                      agder.lake_connectivity
                      WHERE
                      wrid in (", toString(eb_waterregionID, sep=','),") AND
                      to_lake in (", toString(waterbodyID, sep=','),") AND
                      downstream_slope_max_max <= ", slope_barrier, ";", sep='')
  res <- dbGetQuery(db_conection, sql_string)
  res
}
# Eksempel
upstream_lakes <- get_reachable_upstream_lakes(con, unique(wbid_wrid[,1][1:100]), unique(wbid_wrid[,2][1:100]), slope_barrirer)

### Hent ut data frame med kombinasjon av waterbodyID (kolonne 1) og sannsynlighet for tilgjengelighet for gjedde (likelihood)
get_reachable_lakes_liklihood_pike <- function(db_conection, waterbodyID) {
  sql_string <- paste("SELECT DISTINCT ON (accessible_lake) * FROM (
    SELECT (accessible_lakes_pike(wbid)).*
        FROM
            (SELECT unnest(ARRAY[", toString(waterbodyID, sep=','),"]) AS wbid) AS l
    ) AS al
WHERE likelihood > 0.000001
ORDER BY accessible_lake ASC, likelihood DESC;", sep='')
  res <- dbGetQuery(db_conection, sql_string)
  res
}
# Eksempel
accessible_lakes_likelihood_pike <- get_reachable_upstream_lakes(con, unique(wbid_wrid[,2][1:100]))

get_reachable_lakes <- function(db_conection, waterbodyID, slope_barrier) {
    sql_string <- paste("SELECT DISTINCT ON (accessible_lakes)
accessible_lakes_threshold(lake, ", slope_barrier, ") AS accessible_lakes FROM (SELECT unnest(ARRAY[", toString(waterbodyID, sep=','),"]) AS lake) AS x;", sep='')
    res <- dbGetQuery(db_conection, sql_string)
    res
}
# Eksempel (get lakes with up to 7 deg. upstream slope)
reachable_lakes <- get_reachable_upstream_lakes(con, unique(wbid_wrid[,2][1:100]), 700)
# Eksempel (get only downstream or adjacent lakes)
reachable_lakes <- get_reachable_upstream_lakes(con, unique(wbid_wrid[,2][1:100]), 0)

# Funksjon for å hente ut spårednins-sannsynlighet for gjedde
# Returns table with two columns "accessible_lake (waterBodyID) og sansynlighet for tilgjengelighet for gjedde (0.0 -1.0)
get_reachable_lakes_pike <- function(db_conection, waterbodyID) {
    sql_string <- paste("SELECT * FROM (SELECT DISTINCT ON (accessible_lake) (accessible_lakes_pike_test(lake)).* FROM 
      (SELECT unnest(ARRAY[", toString(waterbodyID, sep=','),"]) AS lake) AS l) AS y WHERE likelihood > 0.00001;", sep='')
    res <- dbGetQuery(db_conection, sql_string)
    res
}
# Eksempel (get lakes with up to 7 deg. upstream slope)
reachable_lakes_pike <- get_reachable_lakes_pike(con, unique(wbid_wrid[,2][1:100]))
# Eksempel (get only downstream or adjacent lakes)
reachable_lakes_pike <- get_reachable_lakes_pike(con, unique(wbid_wrid[,2][1:100]))


####################################################################################
# Rydde opp
poolReturn(con)
poolClose(pool)
