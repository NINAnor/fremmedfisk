#########################
#      brt-model        #
#########################

#' Genneral information about BRT-models (Boosted Regression Trees)
#' Brt-models are a type of regression model technique which is very flexible. It thus has some benefits when fitting ecological data, which is usually non-linear and so on.
#' Given care in the model-fittingg, brt can give predictive advantages over methods as e.g. glm or gam.
#' The following script uses the dismo and gmb packages to optimise a brt model for the given species.
#' Analytically, BRT regularization involves jointly optimizing the number of  trees, learning rate and tree complexity.
#' Here optimizing numbers of trees is done through the gbm.step function, and the script includes a function that tries to aid in the optimazation prosess for learning rate and tree complexity (the get.train.diganostic.func). 
#' The next step is to fit the actual model used for predictions (brt_mod)
#'  Its is recomended to read through "A working guide to boosted regression trees" (Elith, et al 2009), before atempting your first go.

#' Parameters to include
#' @param 



###############################
# part 1: load and filter data
###############################
library(stringr)
library(dismo)
library(dplyr)


source("./R/f_geoselect.R")
outdata <- f_geoselect_inverse_spdf(geoselect="./Data/geoselect_native_Rutilus_rutilus.rds",inndata=outdata) #needs to be adressed
# make spatial selection for model estimation - Norway minus Finnmark, Troms and Nordland. 
# The distribution and native area for finnamark would create a lot of missery
outdata <- outdata[outdata$countryCode =="NO",]
outdata$countryCode<-factor(outdata$countryCode)
outdata <- outdata %>% filter(!(county %in% c("Finnmark","Troms","Nordland")))
outdata$county<-factor(outdata$county)


# remove all populations of focal species where focal species is present at start of time-slot
# i.e. focal_specie. No idea how to do this in data.tables, but it's straith-forward with dplyr
# using the programmable version of functions identified by underscore at the end (in this case filter_)
focal_species_var <- stringr::str_replace(string=focal_species_vec[j], pattern=" ", replacement="_")
select_focal <- paste("!(",focal_species_var,"==1 & introduced==0)")
analyse.df <- outdata %>% dplyr::filter_(select_focal)
analyse.df<-as.data.frame(analyse.df) # convert to data.frame - needed for gbm.step input


###############################
# part 2: Make the brt model
###############################


#It is encuraged to do this with paralell computing speeds the prosess up to some extent.
#Identify cores on current system
cores<-detectCores(all.tests = FALSE, logical = FALSE)
cores

#Create training function for gbm.step
get.train.diganostic.func=function(tree.com,learn){
  #set seed for reproducibility
  k1<-gbm.step(data=analyse.df, 
               gbm.x = c( "distance_to_road_log", "dist_to_closest_pop_log","county","SCI","eurolst_bio10", "area_km2_log"), #Include variables at will here
               gbm.y = "introduced",
               family = "bernoulli", 
               tree.complexity = tree.com,
               learning.rate = learn,
               bag.fraction = 0.6,
               prev.stratify=TRUE,
               n.folds=10,
               n.trees=500,
               step.size=100,
               silent=TRUE,
               plot.main = FALSE,
               n.cores=cores)
  
  k.out=list(interaction.depth=k1$interaction.depth,
             shrinkage=k1$shrinkage,
             n.trees=k1$n.trees,
             AUC=k1$self.statistics$discrimination,
             cv.AUC=k1$cv.statistics$discrimination.mean,
             deviance=k1$self.statistics$mean.resid,
             cv.deviance=k1$cv.statistics$deviance.mean)  
  return(k.out)
}

#define complexity and learning rate
tree.complexity<-c(1:9)
learning.rate<-c(0.01,0.025,0.005,0.0025,0.001)

#setup parallel backend to use n processors
cl<-makeCluster(cores)
registerDoParallel(cl)

#Run the actual function
foreach(i = tree.complexity) %do% {
  foreach(j = learning.rate) %do% {
    nam=paste0("gbm_tc",i,"lr",j)
    assign(nam,get.train.diganostic.func(tree.com=i,learn=j))
    
  }
}

#Stop parallel
stopCluster(cl)
registerDoSEQ()

#Find all item in workspace that contain "gbm_tc"
train.all<-ls(pattern="gbm_tc")

#cbind each list that contains "gbm_tc"
train.results<-list(do.call(cbind,mget(train.all)))

#Place in a data frame
train.results<- do.call(rbind, lapply(train.results, rbind))
train.results <- data.frame(matrix(unlist(train.results),ncol=7 , byrow=T))

#Change column names
colnames(train.results)<-c("TC","LR","n.trees", "AUC", "cv.AUC", "dev", "cv.dev")

#Round 4:7 down to 3 digits
train.results[,4:7]<-round(train.results[,4:7],digits=3)

#Sort by cv.dev, cv.AUC, AUC
train.results<-train.results[order(train.results$cv.dev,-train.results$cv.AUC, -train.results$AUC),]

train.results #Includes a dataframe with ordered (numbered) choice based on AUC cv.dev and cv.AUC, be aware that there are mutiple ways of judging the models...

# Use best parametrization from train.results 

brt_mod<-gbm.fixed(data=analyse.df, gbm.x = c( "distance_to_road_log", "dist_to_closest_pop_log","county","SCI","eurolst_bio10", "area_km2_log"), gbm.y = "introduced",family = "bernoulli",tree.complexity = 10, learning.rate = 0.001,bag.fraction = 0.60,n.trees=3700)
names(brt_mod$gbm.call)[1] <- "dataframe"

# save modell object as .rds
saveRDS(brt_mod,paste("./Data/brt_mod_",focal_species_var,".rds",sep=""))




