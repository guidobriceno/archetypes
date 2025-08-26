##%######################################################%##
#                                                          #
####                    Prepare data                    ####
#                                                          #
##%######################################################%##

library(kamila)
library(terra)
library(data.table)
library(stringr)
library(ggplot2)
library(mclust)

# Open metric files for each site
# ---------------------------------------------------------------- - - -
# Cotriguacu (7533)
Cotriguacu_metrics <- readRDS("results/Cotriguacu/Cotriguacu_metrics.rds")
# Guaviare (27797)
Guaviare_metrics <- readRDS("results/Guaviare/Guaviare_metrics.rds")
# Paragominas (55740 cells)
Paragominas_metrics <- readRDS("results/Paragominas/Paragominas_metrics.rds")

# All sites on one file
metrics_all <- rbindlist(list(Cotriguacu_metrics, Guaviare_metrics, Paragominas_metrics))
metrics_all<-na.omit(metrics_all)
metrics_all[, `:=`(
  def_speed = 100*def_speed/max(def_speed),
  deg_speed = 100*deg_speed/max(deg_speed)
)]
# Area interest
# metrics_all <- metrics_all[def_baseline>5 | def_forest_left<100] # sup 2% de deg ou def ou reg
metrics_all <- metrics_all[def_forest_loss>5 | deg_forest_loss>5 | reg_forest_gain>5] # sup 2% de deg ou def ou reg

# qualitative metrics as factor
cat_df <- metrics_all[, .SD, .SDcols = patterns("activeness")]
cat_df <- cat_df[, lapply(.SD, factor)]

# quantitative metrics 
con_df <- metrics_all[, .SD, .SDcols = !patterns("activeness")][, -c("cell","site")] 

# réduire ou centrer-réduire ???
#con_df <- data.table(base::scale(con_df)) # centré réduit
#con_df <- con_df[, lapply(.SD, function(x) x/max(x))] # réduit


##%######################################################%##
#                                                          #
####                    Classify                        ####
#                                                          #
##%######################################################%##

# Selecting the number of clusters (Kamila) 
# ---------------------------------------------------------------- - - -
#ClassSelectionKamila <- kamila(con_df,cat_df,numClust=3:6, numInit=10,calcNumClust = "ps", numPredStrCvRun = 15, predStrThresh = 0.5)
#png(file="results/nbe_class.png")
#plot(3:6, ClassSelectionKamila$nClust$psValues)
#dev.off()
# Apply the Kamila classification on all sites
# ---------------------------------------------------------------- - - -
ncluster = 4
class_kamila <-  kamila(con_df, cat_df, ncluster, 10, maxIter = 100)
metrics_all$cluster <- class_kamila$finalMemb

#sensitivity test
ncluster = 4
ntest <- 100
clTest <- matrix(0,nrow(con_df),ntest)
for( i in 1:ntest){
  clTest[,i] <- kamila(con_df, cat_df, ncluster, numInit=10, maxIter = 100)$finalMemb
}

concord <- rep(0,ntest*(ntest-1)/2)
k <- 1
for(i in 1:(ntest-1)){
  for(j in (i+1):ntest){
    print(k)
    concord[k] <- mclust::adjustedRandIndex(clTest[,i],clTest[,j])
    k <- k+1
  }
}
plot(density(concord),t="l")
summary(concord)

# Convert to images
sites <- c("Paragominas", "Guaviare", "Cotriguacu")

for (i in 1:length(sites)){
  Site <- sites[i]
  img <- rast(stringr::str_glue("results/{Site}/baseline_def.tif"))
  empty_img <- rast(img, nlyrs=1, vals=NA)
  Clust_site <- metrics_all[site == Site,]
  set.values(empty_img, Clust_site$cell, Clust_site$cluster)
  plot(empty_img, colNA="blue")
  writeRaster(empty_img, stringr::str_glue("results/{Site}/Archetype_{Site}.tif"), overwrite=T)
}

saveRDS(metrics_all,'results/metrics_all.rds')
# buffer -2km

# Boxplot by archetypes
# ---------------------------------------------------------------- - - -

# names(metrics_Pgm)<-c("Activeness of deforestation","Activeness of degradation","Kamila (4cl)",
#                       'Baseline forest',"Percentage of regrowth","Forest left" ,"Percentage of forest loss",
#                       "Percentage of degradation","Speed of forest loss", "Speed of degradation")

# to put in order according to the label of the classes
Titles = c("Past gradual frontiers", "Consolidated Frontiers", "Rampant frontiers", "Vulnerable frontiers")

for (i in 1:ncluster){
  # class1
  metrics_class <- metrics_all[which(metrics_all$cluster == i),]
  metrics_class <- metrics_class[,-c(1,2,7,12)]
  metrics_class <- melt(metrics_class, id.vars = "site")
  
  image <-ggplot(metrics_class, aes(x=variable, y=value, fill=site))+
    geom_boxplot()+
    labs(x=" ", y=" Percent",title=Titles[i])+theme_bw()+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 2))
  
  Site = site[i]
  ggsave(stringr::str_glue("results/{Site}"), image, path= stringr::str_glue("results/{Site}"))  
  
}


