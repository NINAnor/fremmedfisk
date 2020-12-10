##############################################
# Function for predicting introduction-events#
##############################################

#' Creates function which produces a dataframe from outdata, including introduction probabilities and introduction events based on these

#' @param outdata data-frame comming from "Wrangle.R"
#' @param start_year Year of staring recording as introductions, same as in wrangle.R
#' @param end_year Year of ending recording as introductions, same as in wrangle.R
#' @param brt_mod A boosted regression tree model from make_boosted_glm
#' @param species The common name of the species
#' @param temp_inc An optional yearly temp increment if nothing temp inc = 0. Only included as an aditive
#' @return A data.frame with indtroduction probabilities on a yearly basis, including introduction events based on these probabilities (rbinom)









# use boosted regression tree
f_predict_introduction_events_gmb <-function(outdata,brt_mod,species,summer_temp_inc,winter_temp_inc,start_year, end_year){
  require(gbm)
  require(stringr)
  require(dplyr)
  species2=stringr::str_replace(species, " ", "_")
  outdata$distance_to_road<-outdata$distance_to_road+0.01
  outdata$distance_to_road_log<-log(outdata$distance_to_road)
  outdata$dist_to_closest_pop_log<-log(outdata$dist_to_closest_pop)
  data_no_pike<-outdata[outdata[[species2]] == 0,]
  data_w_pike<-outdata[outdata[[species2]] == 1,]
  data_no_pike$eurolst_bio10<- data_no_pike$eurolst_bio10+summer_temp_inc #increase annual mean summer temperature for the period (scale to eurolist_bio10)
  data_no_pike$eurolst_bio11<- data_no_pike$eurolst_bio11+winter_temp_inc
  data_no_pike$prob_introduction<-predict.gbm(brt_mod,data_no_pike,n.trees=brt_mod$gbm.cal$best.trees, type="response")#make probabilties over the full period
  data_no_pike$prob_introduction<-data_no_pike$prob_introduction/(end_year-start_year) #annual estimates based on total introductions in time period
  data_no_pike$introduced<-rbinom(length(data_no_pike$prob_introduction), size = 1, prob=data_no_pike$prob_introduction)
  data_no_pike[[species2]]<-ifelse(data_no_pike$introduced==1, 1,0)
  # Add back introductions from before simulation
  outdata<-bind_rows(data_no_pike,data_w_pike)
  return(outdata)
}


