library(postGIStools)
library(RPostgreSQL)
#Set connection parameters
pg_drv <- dbDriver("PostgreSQL")
pg_user <- rstudioapi::askForPassword("enter username")
pg_password <- rstudioapi::askForPassword("enter psw")
pg_host <- "vm-srv-finstad.vm.ntnu.no"
pg_db<-"nofa"

#Initialise connection
con<-dbConnect(pg_drv,dbname=pg_db,user=pg_user, password=pg_password,host=pg_host)

geoselect_native_gjedde <- get_postgis_query(con, 'SELECT ST_Transform(geom, 4326) AS geom, ogc_fid, occurrence  FROM temporary.gjedde ',  geom_name = "geom")
geoselect_native_mort <- get_postgis_query(con, 'SELECT ST_Transform(geom, 4326) AS geom, ogc_fid, occurrence  FROM temporary.mort ',  geom_name = "geom")
geoselect_native_soerv <- get_postgis_query(con, 'SELECT ST_Transform(geom, 4326) AS geom, ogc_fid, occurrence  FROM temporary.soerv ',  geom_name = "geom")
geoselect_native_oerekyte <- get_postgis_query(con, 'SELECT ST_Transform(geom, 4326) AS geom, ogc_fid, occurrence  FROM temporary.oerekyte ',  geom_name = "geom")
