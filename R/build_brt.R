#!/usr/bin/env Rscript

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


if (!require("optparse", character.only=T, quietly=T)) {
  install.packages("optparse")
  require("optparse")
}
library("optparse")

option_list = list(
  make_option(c("-s", "--species_group"), type="character", default=NULL,
              help="Species group to process", metavar="character"),
  make_option(c("-w", "--workdir"), type="character", default="/data/R/Prosjekter/13845000_invafish/",
              help="Name of the database to use", metavar="character"),
  make_option(c("-c", "--cores"), type="integer", default=4,
              help="Number of cores to use for parallel processing", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Run with arguments:
# Rscript --vanilla build_brt.R -s "Agn" -workdir "/data/R/Prosjekter/13845000_invafish/fremmedfisk_opsjon/simulering/" -cores 20
#
# Run script on the command line in parallel for all species groups:
# species_groups=$(psql -h "vm-srv-wallace.vm.ntnu.no" -d nofa -U "stefan" -At -c "SELECT DISTINCT ON (species_group) species_group FROM fremmedfisk.species_groups;")
# echo "$species_groups" | awk '{print("Rscript --vanilla ~/fremmedfisk/R/build_brt.R -s " $1 " -w /data/R/Prosjekter/13845000_invafish/fremmedfisk_opsjon/simulering/ -c 9 \0")}' | xargs -0 -P 5 -I{} bash -c "{}"

# To run in IDE/RStudio skip option parsing and use this opt list
# opt <- list("species_group"="Predator", "workdir"="/data/R/Prosjekter/13845000_invafish/fremmedfisk_opsjon/simulering/", "cores"=40)

# To clarify
# - "start populations" kun introduserte eller alle, hva med ukjente
# - Hvordan håndtere forekomst av flere arter fra en gruppe i et vann
#   Relevant i f_predict_introduction_events_gmb
# - Hvordan håndtere sampling bias i rom og tid (begge er ofte skev fordelt)
#   Relevant i f_predict_introduction_events_gmb og kanskje også bygging av gbm
# Parallelize simulations (they are slow (mainly predictions with many lakes)
# Split script into three parts (1, data collection/wrangling, 2. model building,
#   3. simulation and writing of data) so it can be run on the command line to avoid timeouts



if (!require('devtools', character.only=T, quietly=T)) {
  # Requires OpenSSL on the system
  install.packages('devtools')
  require('devtools')
}

if (!require('easypackages', character.only=T, quietly=T)) {
  devtools::install_github("jakesherman/easypackages")
  require('easypackages')
  }

library(easypackages)

req_packages <- c('stringr',
                  'gbm',
                  'dismo',
                  'dplyr',
                  'doParallel')

for (p in req_packages) {
  if(!require(p, character.only=T, quietly=T)){
    install.packages(p)
    print(p)
  }
}
libraries(req_packages)

# Potential arguments
# test if there is at least one argument: if not, return an error

workdir <- opt$workdir
focal_group <- opt$species_group
#It is encuraged to do this with paralell computing speeds the prosess up to some extent.
#Identify cores on current system
cores <- opt$cores #detectCores(all.tests = FALSE, logical = FALSE) - 8

# Load wrangled data
load(paste0(workdir,"_",focal_group,".RData"))


###############################
# part 1: load and filter data
###############################
#focal_species_var<-"gjedde"
#source("./R/f_geoselect.R")
#lake_env_species <- f_geoselect_inverse_spdf(geoselect="./Data/geoselect_native_Rutilus_rutilus.rds",inndata=lake_env_species) #needs to be adressed
# make spatial selection for model estimation - Norway minus Finnmark, Troms and Nordland.
# The distribution and native area for finnamark would create a lot of missery

# e.g.
# lake_env_species <- lake_env_species_data_gjedde[lake_env_species_data_gjedde$countryCode =="NO",]
# or
# lake_env_species <- lake_env[lake_env$observation_lake == 1,] # inndata_timeslot[inndata_timeslot$countryCode =="NO",]  #lake_env[lake_env$countryCode =="NO",]

#lake_env_species <- merge(inndata_timeslot[,c('waterBodyID')], lake_env, by='waterBodyID', all.y=FALSE)

#lake_env_species$countryCode <- factor(lake_env_species$countryCode)
#lake_env_species <- lake_env_species %>% filter(!(county %in% c("Finnmark","Troms","Nordland")))
#lake_env_species <- lake_env_species %>% filter((county %in% c("Aust-Agder","Vest-Agder","Telemark","Rogaland")))
#lake_env_species$county <- factor(lake_env_species$county)
# lake_env_species <- lake_env_species[lake_env_species$minimumElevationInMeters>0,]
# lake_env_species$n_pop <- NA
# lake_env_species$n_pop <- ifelse(lake_env_species$waterBodyID %in% geoselect_no_gjedde_pop_5000$waterBodyID,geoselect_no_gjedde_pop_5000$count,lake_env_species$n_pop)
#making a rougher latitude variable to make it les specific
#lake_env_species$decimalLatitude_R<-round(lake_env_species$decimalLatitude,1)
covariates <- c("distance_to_road_log","elevation" ,"decimalLatitude_R","dist_to_closest_pop_log", "SCI", "eurolst_bio10", "buffer_5000m_population_2011" ,"area_km2_log", "n_pop")

#focal_species_vec <- unique(lake_env_species$focal_species)
# remove all populations of focal species where focal species is present at start of time-slot
# i.e. focal_specie. No idea how to do this in data.tables, but it's straith-forward with dplyr
# using the programmable version of functions identified by underscore at the end (in this case filter_)
#focal_species_var <- stringr::str_replace(string=focal_species_vec, pattern=" ", replacement="_")
#select_focal <- paste("!(",focal_species_var,"==1 & introduced==0)")
#analyse.df <- lake_env_species %>% dplyr::filter_(select_focal)
analyse.df <- as.data.frame(lake_env_species) # convert to data.frame - needed for gbm.step input

# Remove NoData
for (n in append(covariates, "introduced")) {
  analyse.df <- analyse.df[!is.na(analyse.df[n]),]
}
analyse.df$introduced <- as.integer(analyse.df$introduced)

###############################
# part 2: Make the brt model
###############################


# Outer loop has 9 items, the inner 5
#Create training function for gbm.step
get.train.diganostic.func=function(tree.com,learn,indf){

  k.out=list(interaction.depth=NA,
             shrinkage=NA,
             n.trees=NA,
             AUC=NA,
             cv.AUC=NA,
             deviance=NA,
             cv.deviance=NA)

  #set seed for reproducibility
  k1<-try(gbm.step(data=indf,
                   gbm.x = covariates, #  Include variables at will here,"county"
                   gbm.y = "introduced",
                   family = "bernoulli",
                   tree.complexity = tree.com,
                   learning.rate = learn,
                   bag.fraction = 0.8,
                   prev.stratify=TRUE,
                   n.folds=10,
                   n.trees=500,
                   step.size=100,
                   silent=TRUE,
                   plot.main = FALSE,
                   n.cores=1))

  if(exists("k1")) {
    if(! is.vector(k1) & !is.null(k1)) {
      k.out=try(list(interaction.depth=k1$interaction.depth,
                     shrinkage=k1$shrinkage,
                     n.trees=k1$n.trees,
                     AUC=k1$self.statistics$discrimination,
                     cv.AUC=k1$cv.statistics$discrimination.mean,
                     deviance=k1$self.statistics$mean.resid,
                     cv.deviance=k1$cv.statistics$deviance.mean))
    }
  }
  return(k.out)
}

#define complexity and learning rate
tree.complexity<-c(1:9)
learning.rate<-c(0.01, 0.025, 0.005, 0.0025, 0.001)
stepsize <- 100

#setup parallel backend to use n processors
cl<-parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)

start.time <- Sys.time()
#Run the actual function
gbms <- foreach(i=tree.complexity, .packages = c('gbm', 'dismo', 'doParallel')) %:%
  foreach(j=learning.rate, .packages = c('gbm', 'dismo', 'doParallel')) %dopar% {
    try(get.train.diganostic.func(tree.com=i,learn=j,indf=analyse.df))
  }
end.time <- Sys.time()

#exec_time <- end.time - start.time
end.time - start.time

#Stop parallel
stopCluster(cl)
registerDoSEQ()

# Create data frame for collecting training results
train.results.par <- data.frame(tc=numeric(),
                                lr=numeric(),
                                interaction.depth=numeric(),
                                shrinkage=numeric(),
                                n.trees=numeric(),
                                AUC=numeric(),
                                cv.AUC=numeric(),
                                deviance=numeric(),
                                cv.deviance=numeric()
)

# Collect training results
for (i in 1:length(tree.complexity)) {
  for (j in 1:length(learning.rate)) {
    train.results.par[nrow(train.results.par) + 1,] <- list(
      tc=ifelse(is.null(tree.complexity[i]),NA,tree.complexity[i]),
      lr=ifelse(is.null(learning.rate[j]),NA,learning.rate[j]),
      interaction.depth=ifelse(is.null(gbms[[i]][[j]]$interaction.depth),NA,gbms[[i]][[j]]$interaction.depth),
      shrinkage=ifelse(is.null(gbms[[i]][[j]]$shrinkage),NA,gbms[[i]][[j]]$shrinkage),
      n.trees=ifelse(is.null(gbms[[i]][[j]]$n.trees),NA,gbms[[i]][[j]]$n.trees),
      AUC=ifelse(is.null(gbms[[i]][[j]]$AUC),NA,gbms[[i]][[j]]$AUC),
      cv.AUC=ifelse(is.null(gbms[[i]][[j]]$cv.AUC),NA,gbms[[i]][[j]]$cv.AUC),
      deviance=ifelse(is.null(gbms[[i]][[j]]$deviance),NA,gbms[[i]][[j]]$deviance),
      cv.deviance=ifelse(is.null(gbms[[i]][[j]]$cv.deviance),NA,gbms[[i]][[j]]$cv.deviance)
    )
  }
}



#Find all item in workspace that contain "gbm_tc"
#train.all<-ls(pattern="gbm_tc")

#cbind each list that contains "gbm_tc"
#train.results<-list(do.call(cbind,mget(train.all)))

#Place in a data frame
#train.results<- do.call(rbind, lapply(train.results, rbind))
#train.results <- data.frame(matrix(unlist(train.results),ncol=7 , byrow=T))

#Change column names
#colnames(train.results)<-c("TC","LR","n.trees", "AUC", "cv.AUC", "dev", "cv.dev")

#Round 4:7 down to 3 digits
train.results.par[,6:9] <- round(train.results.par[,6:9],digits=3)

#Sort by cv.dev, cv.AUC, AUC
train.results.par <- train.results.par[order(train.results.par$cv.deviance,-train.results.par$cv.AUC, -train.results.par$AUC),]

# Results deviate when using %do% and %dopar%
train.results.par #Includes a dataframe with ordered (numbered) choice based on AUC cv.dev and cv.AUC, be aware that there are mutiple ways of judging the models...

# save training results as .rds for manual inspection
saveRDS(train.results.par, paste0(workdir, "/brt_mod_norway_training_",focal_species_str,".rds"))


best_params <- train.results.par[1,]


# summed_likelihood_for_non_introduction <- prod(1-ifelse(is.na(tmp_trans$prob_introduction),0,tmp_trans$prob_introduction))
# summed_likelihood_for_introduction <- 1 - summed_likelihood_for_non_introduction

# For Agder with 378 rows in data.table
#tc     lr interaction.depth shrinkage n.trees   AUC cv.AUC deviance
#2 0.0050                 2    0.0050     900 0.956  0.950    0.042
#3 0.0025                 3    0.0025    1500 0.963  0.949    0.038
#2 0.0010                 2    0.0010    5200 0.962  0.947    0.041
#6 0.0025                 6    0.0025    1000 0.977  0.941    0.034
#6 0.0010                 6    0.0010    2900 0.976  0.936    0.031
# Use best parametrization from train.results
#
# brt_mod<-gbm.fixed(data=analyse.df, gbm.x = c( "distance_to_road_log", "dist_to_closest_pop_log","SCI","minimumElevationInMeters","buffer_5000m_population_2006" ,"area_km2_log","n_pop"), gbm.y = "introduced",family = "bernoulli",tree.complexity = 9, learning.rate = 0.001,bag.fraction = 1,n.trees=4000)
brt_mod <- gbm.step(data=analyse.df, gbm.x=covariates, gbm.y="introduced", family="bernoulli", tree.complexity=best_params$tc, step.size=100, learning.rate=best_params$lr, n.trees=ifelse(best_params$n.trees < (10 * stepsize),stepsize, best_params$n.trees - (10 * stepsize)), max.trees=best_params$n.trees + (10 * stepsize))
names(brt_mod$gbm.call)[1] <- "dataframe"

predictors<-gbm.simplify(brt_mod, n.folds = 10, n.drops = "auto", alpha = 1, prev.stratify = TRUE,
                         eval.data = NULL, plot = TRUE)
# Plot suggests to possibly also drop a second predictor (buffer_5000m_population_2006)
brt_mod_simp<-gbm.step(data=analyse.df, keep.data = TRUE, gbm.x=predictors$pred.list[[1]], gbm.y="introduced", family="bernoulli", tree.complexity=best_params$tc, step.size=100, learning.rate=best_params$lr, n.trees=best_params$n.trees - (5 * stepsize), max.trees=best_params$n.trees - (15 * stepsize))
#gbm.step(data=analyse.df, gbm.x = predictors$pred.list[[1]], gbm.y = "introduced",family = "bernoulli",tree.complexity = 8,step.size=100 ,learning.rate = 0.01,n.trees=1500)

# save model object as .rds
saveRDS(brt_mod_simp,paste0(workdir, "/brt_mod_norway_",focal_species_str,".rds"))

