####################################################
#
# Postprosessing of sim-output in order to get some
# usefull format
#
####################################################

# output pr lake
f_sim_output_lake <- function(sim_output,inndata_sim1,Nsims){
  require(dplyr)
  sim_output_lake <- sim_output %>%
    group_by(intro) %>%
    summarize(n_intro=n(),intro_is_secondary=tail(names(sort(table(intro_is_secondary))),1)) %>%
    rename(waterBodyID=intro)
  sim_output_lake$p_intro <- sim_output_lake$n_intro / (Nsims)
  sim_output_lake <- left_join(inndata_sim1[c("waterBodyID","utm_x","utm_y")],
                               sim_output_lake,by="waterBodyID")# %>%
    #rename(lat=decimalLatitude,long=decimalLongitude)
  sim_output_lake$p_intro[is.na(sim_output_lake$p_intro)] <- 0
  sim_output_lake$n_intro[is.na(sim_output_lake$n_intro)] <- 0
  return(sim_output_lake)
}



# raw output spatially enabled
f_sim_output_raw_spatial <- function(sim_output,inndata_sim1) {
  require(dplyr)
  sim_output_spatial <- left_join(sim_output,inndata_sim1[c("waterBodyID","decimalLongitude","decimalLatitude")],
                                  by=c("intro"="waterBodyID")) %>%
    rename(lat=decimalLatitude,long=decimalLongitude)
  return(sim_output_spatial)
}


#................................................
# N antall bestander med focal_species
#
# Function aggregates sim_output to cummulative
# number of populations present of focal species
# at the time of each sim-run.
#------------------------------------------------

f_simout_cummulative_N <- function(sim_output,inndata_sim1,species_var){

  # Number of populations at start of simulation
  n_pop_start <- sum(inndata_sim1[[species_var]])

  # create data frame with N number of intro's pr sim-run and time-slot
  sim_output_n_j <- sim_output %>%
    group_by(sim_j,time_slot_i) %>%
    summarize(N=(n()+n_pop_start))

  # calculate cummulative sums of pop's along time_slot i within j
  # order sim_output_n_j first just to be sure...
  sim_output_n_j <- sim_output_n_j[order(sim_output_n_j$sim_j,sim_output_n_j$time_slot_i),]
  sim_runs <- unique(sim_output_n_j$sim_j) # index for for loop
  for(j in 1:length(sim_runs)){
    # create vector to store output at first iteration
    if(j==1){
      N_cum <- numeric()
    }
    cumsum_tmp <- cumsum(sim_output_n_j$N[sim_output_n_j$sim_j==sim_runs[j]])
    N_cum <- append(N_cum,cumsum_tmp)
  }
  sim_output_n_j$N_cum <- N_cum # attach cumsums to data.frame

  return(sim_output_n_j)
} # end of function

# example
# species_var <- "Esox_lucius"
# f_simout_cummulative_N(sim_output,inndata_sim1,species_var)

#..................................................................
#
# Write output to DB
#.............................................................

f_write_simresult_to_db <- function(dataToWrite,nameOfTable){
  require(RPostgreSQL)

  # connect to db
  con <- dbConnect(drv=dbDriver("PostgreSQL"),
                   dbname = "nofa",
                   host =rstudioapi::askForPassword("Please enter server adress"),
                   port = 5432,
                   user =rstudioapi::askForPassword("Please enter your user"),
                   password =rstudioapi::askForPassword("Please enter your psw"),
                   options="-c search_path=temporary"
                   )
  # write table to db
  dbWriteTable(con, name=paste(nameOfTable), value=dataToWrite,overwrite=TRUE)

  # grant priviliges on new table
  dbSendStatement(conn=con,paste("GRANT SELECT ON", noquote(nameOfTable) ," TO nofa_reader;"))
  dbSendStatement(conn=con,paste("GRANT ALL PRIVILEGES ON", noquote(nameOfTable) ," TO nofa_admin;"))
  dbDisconnect(con)
}

# # example (require user with write permission to temporary schema)
# dataToWrite <- sim_output_pr_lake
# nameOfTable <- "simulation_out_test"
# f_write_simresult_to_db(dataToWrite=sim_output_pr_lake,nameOfTable="simulation_out_test")

