##############################################
# Function for predicting introduction-events#
##############################################

#' Creates function which produces a dataframe from outdata, including introduction probabilities and introduction events based on these

#' @param outdata data-frame comming from "Wrangle.R"
#' @param brt_mod A boosted regression tree model from make_boosted_glm
#' @param species The common name of the species
#' @param temp_inc An optional yearly temp increment if nothing temp inc = 0. Only included as an aditive 
#' @param start_year Year of staring recording as introductions, same as in wrangle.R
#' @return A data.frame with indtroduction probabilities on a yearly basis, including introduction events based on these probabilities (rbinom)

# use boosted regression tree
f_predict_introduction_events_gmb <-function(outdata,brt_mod,species,temp_inc,start_year, covariates){
  require(gbm)
  require(stringr)
  require(dplyr)
  species2=stringr::str_replace(species, " ", "_") 
  #outdata$distance_to_road<-outdata$distance_to_road+0.01 
  #outdata$distance_to_road_log<-log(outdata$distance_to_road)
  #outdata$dist_to_closest_pop_log<-log(outdata$dist_to_closest_pop)

  data_no_species<-outdata[outdata[[species2]] == 0,]
  data_w_species<-outdata[outdata[[species2]] == 1,]
  data_no_species$eurolst_bio01 <- data_no_species$eurolst_bio10+temp_inc #increase annual mean summer temperature for the period (scale to eurolist_bio10)
  data_no_species$prob_introduction<-predict.gbm(brt_mod,data_no_species,n.trees=brt_mod$gbm.cal$best.trees, type="response")#make probabilties over the full period
  #data_no_species$prob_introduction<-data_no_species$prob_introduction #annual estimates based on total introductions in time period
  data_no_species$sim_introduced<-rbinom(length(data_no_species$prob_introduction), size = 1, prob=data_no_species$prob_introduction)
  #data_no_species[[species2]]<-ifelse(data_no_species$introduced==1,1,0)
  data_w_species$sim_introduced<-0
  outdata<-bind_rows(data_no_species,data_w_species)
  
  return(outdata)
}

# use boosted regression tree
f_predict_introduction_events_gmb_group <-function(outdata,brt_mod,species,temp_inc,start_year){
  require(gbm)
  require(stringr)
  require(dplyr)
  species2=stringr::str_replace(species, " ", "_") 
  #outdata$distance_to_road<-outdata$distance_to_road+0.01 
  #outdata$distance_to_road_log<-log(outdata$distance_to_road)
  #outdata$dist_to_closest_pop_log<-log(outdata$dist_to_closest_pop)
  
  data_no_species<-outdata[outdata[[species2]] == 0,]
  data_w_species<-outdata[outdata[[species2]] == 1,]
  data_no_species$eurolst_bio01<- data_no_species$eurolst_bio10+temp_inc #increase annual mean summer temperature for the period (scale to eurolist_bio10)
  data_no_species$prob_introduction<-predict.gbm(brt_mod,data_no_species,n.trees=brt_mod$gbm.cal$best.trees, type="response")#make probabilties over the full period
  #data_no_species$prob_introduction<-data_no_species$prob_introduction #annual estimates based on total introductions in time period
  data_no_species$introduced<-rbinom(length(data_no_species$prob_introduction), size = 1, prob=data_no_species$prob_introduction)
  data_no_species[[species2]]<-ifelse(data_no_species$introduced==1,1,0)
  outdata<-bind_rows(data_no_species,data_w_species)
  
  return(outdata)
}
