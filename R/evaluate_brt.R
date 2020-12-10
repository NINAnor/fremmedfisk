cols <- c("n_pop", "dist_to_closest_pop_log", "area_km2_log", "distance_to_road_log", "eurolst_bio10", "buffer_5000m_population_2006", "dataset")

predict <- transform(lake_env, dataset=rep('prediction.df', length(lake_env[,1])))[cols]
train <- transform(analyse.df, dataset=rep('analyse.df', length(analyse.df[,1])))[cols]
compare <- rbind(train, predict)

library(lattice)
library(gtools)
# Open a pdf file
par(mfrow=c(length(cols)-1), 2)

pdf("~/brt_mod_simp_plots.pdf")

par(mfrow=c(1, 1))

for (c in cols[1:(length(cols)-1)]) {
print(c)
hist(analyse.df[,c],
     ylab=c,
     xlab="analyse.df",
     main=NULL)
hist(lake_env[,c],
     ylab=c,
     xlab="predict.df",
     main=NULL)
}

dfs_train <- as.data.frame(summary(train))
dfs_predict <- as.data.frame(summary(predict[complete.cases(predict),]))
compare <- cbind(dfs_predict, dfs_train)[,c(2,3,6)]
names(compare) <- c('variable', 'prediction_space', 'training_space')

grid::grid.newpage()
tt <- gridExtra::ttheme_default(base_size=8)
gridExtra::grid.table(compare, theme=tt)

lapply(1:(length(cols)-1), function(c) plot(brt_mod_simp, i=c))


interact <- combinations((length(cols)-1),2,cols[1:(length(cols)-1)])
#for (i in 1:length(interact[,1])) {
#  plot(brt_mod_simp, i=c(interact[,1][i], interact[,2][i]))
#  }

lapply(1:length(interact[,1]), function(i) plot(brt_mod_simp, i=c(interact[,1][i], interact[,2][i])))
# Close the pdf file
dev.off()
