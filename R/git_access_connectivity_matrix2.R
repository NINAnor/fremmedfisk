#!/usr/bin/env Rscript
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

# slope_barrirer <- 700

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

### Hent ut alle vannregioner (for Agder (mange inneholder ingen innsj?er)
get_wrids <- function(db_conection, connectivity_schema, connectivity_table) {
  sql_string <- paste0('SELECT DISTINCT ON (wrid) wrid FROM ', connectivity_schema, '.', connectivity_table, ';') 
  res <- dbGetQuery(db_conection, sql_string)[,1]
  res
}
# Eksempel
# wrids <- get_wrids(con, "fremmedfisk", "fremmedfisk_lake_connectivity_summary")

### Hent ut data frame med kombinasjon av waterbodyID for alle innsj?er (kolonne 1) og id for vannregioner (kolonne 2) (her er det kun vannregioner som inneholder innsjøer)
get_wbid_wrid <- function(db_conection, connectivity_schema) {
  sql_string <- paste0('SELECT ecco_biwa_wr AS wrid, id AS "waterBodyID" FROM nofa.lake WHERE ecco_biwa_wr IN (SELECT wrid FROM ', connectivity_schema, '.wrids);')
  res <- dbGetQuery(db_conection, sql_string)
  res
}
# Eksempel
# wbid_wrid <- get_wbid_wrid(con, "fremmedfisk", "fremmedfisk_lake_connectivity_summary")

get_wbid_wrid_array <- function(db_conection, waterBodyID) {
  sql_string <- paste("SELECT array_to_string(array_agg(id), ',') AS \"waterBodyID\", ecco_biwa_wr AS wrid FROM nofa.lake WHERE id IN (", toString(waterBodyID, sep=','), ") GROUP BY ecco_biwa_wr;", sep='')
  res <- dbGetQuery(db_conection, sql_string)
  res
}

### Hent ut data frame med unike kombinasjon av waterbodyID (kolonne 1) og id for innsj?er som ligger nedstr?ms (kolonne 2)
get_downstream_lakes <- function(db_conection, waterbodyID, eb_waterregionID) {
  sql_string <- paste("SET constraint_exclusion = on;
                      SELECT \"lakeID\" AS  \"waterBodyID\", CAST(unnest(string_to_array(downstream_lakes, ',')) AS integer) AS downstream_lakes FROM
                      agder.lake_connectivity_summary WHERE
                      wrid IN (", toString(eb_waterregionID, sep=','), ") AND
                      \"lakeID\" IN (", toString(waterbodyID, sep=','),");", sep='')
  res <- dbGetQuery(db_conection, sql_string)
  res
}
# Eksempel
#downstream_lakes <- get_downstream_lakes(con, unique(wbid_wrid[,1][1:100]), unique(wbid_wrid[,2][1:100]))


### Hent ut data frame med kombinasjon av waterbodyID (kolonne 1) og id for innsj?er som ligger oppst?ms og der skr?ning i forbindelsen er lavere enn slope_barrier (kolonne 2)
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
#upstream_lakes <- get_reachable_upstream_lakes(con, unique(wbid_wrid[,1][1:100]), unique(wbid_wrid[,2][1:100]), slope_barrirer)

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
#accessible_lakes_likelihood_pike <- get_reachable_upstream_lakes(con, unique(wbid_wrid[,2][1:100]))

get_reachable_lakes <- function(db_conection, wrid, waterBodyID, slope_measure, slope_threshold, conmat_schema, conmat_table) {
  sql_string <- paste0("SELECT DISTINCT ON (accessible_lakes_threshold_simple) nofa.accessible_lakes_threshold_simple(wrid, wbid, CAST('", slope_measure, "' AS text), ", slope_threshold, ", CAST('", conmat_schema, "' AS text), CAST('", conmat_table, "' AS text)) FROM (SELECT unnest(ARRAY[", toString(wrid), ']) AS wrid, CAST(unnest(ARRAY[', toString(paste0("'", waterBodyID, "'", sep="")), "]) AS text) AS wbid) AS x")
  res <- dbGetQuery(db_conection, sql_string)
  res
}

# Eksempel (get lakes with up to 7 deg. upstream slope)
#reachable_lakes <- get_reachable_upstream_lakes(con, unique(wbid_wrid[,2][1:100]), 700)
# Eksempel (get only downstream or adjacent lakes)
#reachable_lakes <- get_reachable_upstream_lakes(con, unique(wbid_wrid[,2][1:100]), 0)
get_reachable_lakes_wbid <- function(db_conection, waterbodyID, slope_barrier, conmat_schema, conmat_table) {
  sql_string <- paste0('SELECT c."waterBodyID" AS source, l.to_lake AS acclake
      FROM agder.lake_connectivity AS l,
(
  SELECT ecco_biwa_wr AS wrid, id AS "waterBodyID" FROM nofa.lake WHERE id IN(', toString(waterbodyID, sep=','),')
) AS c
  WHERE l.wrid = c.wrid AND
  l.from_lake = c."waterBodyID" AND l.upstream_slope_max_max <= 700
  UNION ALL
  SELECT c."waterBodyID" AS source, l.from_lake AS acclake
  FROM agder.lake_connectivity AS l,
  (
  SELECT ecco_biwa_wr AS wrid, id AS "waterBodyID" FROM nofa.lake WHERE id IN(', toString(waterbodyID, sep=','), ')
) AS c

      WHERE l.wrid = c.wrid AND
      l.to_lake = c."waterBodyID" AND l.downstream_slope_max_max <= 700')
                       res <- dbGetQuery(db_conection, sql_string)
                       res
}


get_reachable_lakes_wrid <- function(db_conection, wrid, waterbodyID, slope_barrier) {
  sql_string <- paste(
    "--SELECT count(acclake), acclake FROM (
    SELECT l.to_lake AS acclake
    -- , CASE WHEN l.upstream_slope_max_max <= 0 THEN CAST(0 AS smallint)
    -- ELSE l.upstream_slope_max_max END AS slope
    -- , CASE
    -- WHEN l.downstream_slope_max_max <= 0 AND l.upstream_slope_max_max > 0 THEN CAST(''upstreams'' AS character varying(25))
    -- WHEN l.upstream_slope_max_max <= 0 AND l.downstream_slope_max_max > 0 THEN CAST(''downstreams'' AS character varying(25))
    -- ELSE CAST(''up/-donwstreams'' AS character varying(25)) END AS type
    FROM agder.lake_connectivity AS l
    WHERE l.wrid = ", wrid, " AND
    l.from_lake IN (", toString(waterbodyID, sep=','), ") AND l.upstream_slope_max_max <= ", slope_barrier, "
    UNION ALL
    SELECT l.from_lake AS acclake
    -- , CASE WHEN l.upstream_slope_max_max <= 0 THEN CAST(0 AS smallint)
    -- ELSE l.upstream_slope_max_max END AS slope
    -- , CASE
    -- WHEN l.downstream_slope_max_max <= 0 AND l.upstream_slope_max_max > 0 THEN CAST(''upstreams'' AS character varying(25))
    -- WHEN l.upstream_slope_max_max <= 0 AND l.downstream_slope_max_max > 0 THEN CAST(''downstreams'' AS character varying(25))
    -- ELSE CAST(''up/-donwstreams'' AS character varying(25)) END AS type
    FROM agder.lake_connectivity AS l
    WHERE l.wrid = ", wrid, " AND
    l.to_lake IN (", toString(waterbodyID, sep=','), ") AND l.downstream_slope_max_max <= ", slope_barrier, "
    --) AS y
    --GROUP BY acclake
    ;"
    , sep='')
  res <- dbGetQuery(db_conection, sql_string)
  res
}

# Example
# rl <- get_reachable_lakes_wrid(con, 246992, wbid_wrid[,1][wbid_wrid[,2] == 246992][1:100], 700)

# Funksjon for ? hente ut sp?rednins-sannsynlighet for gjedde
# Returns table with two columns "accessible_lake (waterBodyID) og sansynlighet for tilgjengelighet for gjedde (0.0 -1.0)
get_reachable_lakes_pike <- function(db_conection, waterbodyID) {
  sql_string <- paste("SELECT * FROM (SELECT DISTINCT ON (accessible_lake) (accessible_lakes_pike_test(lake)).* FROM
                      (SELECT unnest(ARRAY[", toString(waterbodyID, sep=','),"]) AS lake) AS l) AS y WHERE likelihood > 0.00001;", sep='')
  res <- dbGetQuery(db_conection, sql_string)
  res
}
# Eksempel (get lakes with up to 7 deg. upstream slope)
#reachable_lakes_pike <- get_reachable_lakes_pike(con, unique(wbid_wrid[,2][1:100]))
# Eksempel (get only downstream or adjacent lakes)
#reachable_lakes_pike <- get_reachable_lakes_pike(con, unique(wbid_wrid[,2][1:100]))

####################################################################################
# Rydde opp
#poolReturn(con)
#poolClose(pool)
