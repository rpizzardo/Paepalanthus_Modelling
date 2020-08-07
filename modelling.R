######## Species Distribution Modelling: present and future predictions #######
#### R. C. Pizzardo, MAY 2020 ####

rm(list=ls()) #uncomment to clean your environment

wd <- "C:/Users/raque/Paepalanthus/"
setwd(wd)
library(sdm)
library(dismo)
library(raster)
library(shiny)
library(CoordinateCleaner)
library(maptools)
library(mapview)
library(magrittr)
library(dplyr)
library(usdm)
library(OpenImageR)
library(rJava)
data(wrld_simpl)

### Organizing the data ###
##Getting the distribution points
sp <- read.csv("distribution_points.csv")
dat_cl<-sp

##Flaging the data
sp <- data.frame(sp)
sp[is.na(sp)] = 0
flags <- clean_coordinates(x = sp, 
                           lon = "Lon", 
                           lat = "Lat",
                           species = "species",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros", "seas", "duplicates", "urban"),
                           seas_ref = buffland) # most test are on by default
suma <- summary(flags)

##Cleaning the data
coord_no <- !is.na(sp$lon) | !is.na(sp$lat) # flag records without coordinates
coord_prec <- sp$coordinateUncertaintyInMeters/1000 <= 100 | 
  is.na(sp$coordinateUncertaintyInMeters) #remove records with low coordinate precision
rec_age <- sp$year > 1945 # Record age

#create clean dataset
rm <- flags$.summary
if(length(coord_no) > 0) {rm <- rm & coord_no}
if(length(coord_prec) > 0) {rm <- rm & coord_prec}
if(length(rec_age) > 0) {rm <- rm & rec_age}
dat_cl <-sp%>%filter(rm)

#create output and write to disk
flag_sum <- suma %>%
  t()%>%
  data.frame()%>%
  mutate(
    coordinate_missing = sum(!coord_no),
    coordinate_precision = sum(!coord_prec),
    record_age = sum(!rec_age),
    raw_records = nrow(sp),
    raw_taxa = length(unique(sp$species)),
    clean_records = nrow(dat_cl),
    clean_taxa = length(unique(dat_cl$species)))

i <- dat_cl$lon > 0 & dat_cl$lat > 0
dat_cl$lon[i] <- -1 * dat_cl$lon[i]
dat_cl$lat[i] <- -1 * dat_cl$lat[i]
#dat_cl <- dat_cl[dat_cl$lat < -5, ]

write.csv(dat_cl, "clean_distribution_points.csv")

#_________________________________________________

### Modelling ###
#Geting the variable
bio <- raster::getData('worldclim', var='bio', res=2.5) # for the present
biof_50 <- raster::getData('CMIP5', var='bio',res=2.5,rcp=85,year=50, model='AC') #for 2050
biof_70 <- raster::getData('CMIP5', var='bio',res=2.5,rcp=85,year=70, model='AC') #for 2070

names(biof_50)<-names(bio)
names(biof_70)<- names(bio)
biof_50 <- crop(biof_50, extent(-60, -20, -35, 5))
biof_70 <- crop(biof_70, extent(-60, -20, -35, 5))

#Collinearity test
layersX <- crop(bio, extent(-60, -20, -35, 5)) #Brazil east
v1<-vifstep(layersX)
v2<-vifcor(layersX,th=0.8) #stablished threshold
biom <- exclude(layersX,v2) # excludes variables
current <-stack(biom)

## Creating models 
geo <- dat_cl[,c("species","Lon","Lat")]
geo <- geo[geo$species %in% names(which(table(geo$species) > 3)), ]
species <- names(which(table(geo$species)!=0))

Sys.time() -> start_time_total

for(i in 1:length(species)){
  Sys.time() -> start_time
  sp0 = species[i]
  sp1 <- geo[geo$species==sp0,]
  sp1$species <- 1
  coordinates(sp1) <- ~ Lon + Lat
  
  #Present
  setwd(paste(wd, "distribution_models/1_current/2", sep=""))
  d <- sdmData(species~., sp1, predictors = current, bg = list(n=1000)) # 1000 background points 
  m <- sdm(species~., d, methods =c("glm","brt","rf", "maxent"),        # choosing 4 different methods for modelling
           replication=c('boot'), n=5)                                  # 5 replications for each method = 20 models per species
  
  Stat <- getEvaluation(m,stat=c('TSS', 'AUC'),opt=1, file=paste(sp0, "current.csv", sp="")) # getting the evaluation (TSS and AUC) for each method and model of each species
  write.csv(Stat, file=paste(sp0, "evaluation.csv", sp=""))
  
  
  x <- getVarImp(m, id=1, wtest="test.dep")                         # getting the relative importance of each variable per species
  write.table(data.frame(), file=paste0(sp0, "_varimp.txt"))
  sink(paste0(sp0, "_varimp.txt"))
  print(x)
  sink()
  plot(x,'auc')
  
  
  en <- ensemble(m, current, paste(sp0, "_ensemble.tif", sep="") ,
                 setting=list(method='weighted',stat=c("AUC")))      # ensemble ("fusion model")
  
  #2070
  setwd(paste(wd, "distribution_models/1_future_70/", sep=""))
  enf_70 <- ensemble(m, biof_70, paste(sp0, "_ensemble_70.tif", sep=""),
                     setting = list(method='weighted', stat=c("AUC")))
  
  #2050
  setwd(paste(wd, "distribution_models/1_future_50/", sep=""))
  enf_50 <- ensemble(m, biof_50, paste(sp0, "_ensemble_50.tif", sep=""),
                     setting = list(method='weighted', stat=c("AUC")))
  
  #ploting pdfs for present, 2050 and 2070:  roc curve (AUC), ensemble models and variable importance (graph) 
  setwd(paste(wd, "distribution_models/1_pdfs_current/", sep=""))
  pdf(file=paste(sp0, "_roc", ".pdf", sep=""), width = 8, height = 5)
  roc(m, 1)
  title("glm model", adj=0, add=T)
  roc(m, 2)
  title("brt model", adj=0, add=T)
  roc(m, 3)
  title("rf model", adj=0, add=T)
  roc(m, 4)
  title("maxent model", adj=0, add=T)
  dev.off()
  
  pdf(file=paste(sp0, "_sdm", ".pdf", sep=""), width = 10, height = 5)
  plot(en, zlim=c(0,1))
  plot(sp1, pch=19, cex=0.1, col="tomato1", add=T)
  title(main=paste(sp0, "current distribution", sep=" "))
  dev.off()
  
  pdf(file=paste0(sp0, "_getvarimp", ".pdf"), width = 3, height = 3)
  plot(x, 'cor')
  title("cor")
  plot(x, 'auc')
  title("auc")
  dev.off()
  
  setwd(paste(wd, "distribution_models/1_pdfs_70/", sep=""))
  pdf(file=paste(sp0, "_roc", ".pdf", sep=""), width = 8, height = 5)
  roc(m, 1)
  title("glm model", adj=0, add=T)
  roc(m, 2)
  title("brt model", adj=0, add=T)
  roc(m, 3)
  title("rf model", adj=0, add=T)
  roc(m, 4)
  title("maxent model", adj=0, add=T)
  dev.off()
  
  pdf(file=paste(sp0, "_sdm70", ".pdf", sep=""), width = 10, height = 5)
  plot(enf_70, zlim=c(0,1))
  plot(sp1, pch=19, cex=0.1, col="tomato1", add=T)
  title(main=paste(sp0, "future 2070", sep=" "))
  dev.off()
  
  setwd(paste(wd, "distribution_models/1_pdfs_50/", sep=""))
  pdf(file=paste(sp0, "_roc", ".pdf", sep=""), width = 8, height = 5)
  roc(m, 1)
  title("glm model", adj=0, add=T)
  roc(m, 2)
  title("brt model", adj=0, add=T)
  roc(m, 3)
  title("rf model", adj=0, add=T)
  roc(m, 4)
  title("maxent model", adj=0, add=T)
  dev.off()
  
  pdf(file=paste(sp0, "_sdm50", ".pdf", sep=""), width = 10, height = 5)
  plot(enf_50, zlim=c(0,1))
  plot(sp1, pch=19, cex=0.1, col="tomato1", add=T)
  title(main=paste(sp0, "future 2050", sep=" "))
  dev.off()
  
  Sys.time() -> end_time
  print(c(species[i], "done!"))
  print(end_time-start_time)
}
Sys.time() -> end_time_total
print(end_time_total-start_time_total)


#____________________________________________

### Stacking the models according to the bioregions of campo rupestre ###

##For the present
setwd(paste(wd, "distribution_models/1_current/", sep=""))

#getting the models (south distribution)
south <- c( "Paepalanthus bryoides_ensemble.tif", "Paepalanthus comans_ensemble.tif",
            "Paepalanthus eriophaeus_ensemble.tif", "Paepalanthus macrocephalus_ensemble.tif",
            "Paepalanthus scirpeus_ensemble.tif", "Paepalanthus nigrescens_ensemble.tif") 

south <- stack(south)
south_final <- mean(south) #mean
south_sd <- calc(south, fun=sd) #standad deviation
writeRaster(south_final, "south_final.tif")
writeRaster(south_sd, "south_sd.tif")

#getting the models (north distribution)
north <- c("Paepalanthus barbulatus_ensemble.tif", "Paepalanthus cinereus_ensemble.tif",
           "Paepalanthus harleyi_ensemble.tif", "Paepalanthus macrocaulon var. contasensis_ensemble.tif",
           "Paepalanthus pulvinatus_ensemble.tif", "Paepalanthus spathulatus_ensemble.tif")

north <- stack(north)
north_final <- mean(north) # mean
north_sd <- calc(north, fun=sd) #standad deviation
writeRaster(north_final, "north_final.tif")
writeRaster(north_sd, "north_sd.tif")

#ploting the results for the present
par(mfrow=c(2,2))
plot(north_final)
plot(north_sd)
plot(south_final)
plot(south_sd)

##For 2070
setwd(paste(wd, "distribution_models/1_future_70/", sep=""))

south_70 <- c("Paepalanthus bryoides_ensemble_70.tif", "Paepalanthus comans_ensemble_70.tif",
              "Paepalanthus eriophaeus_ensemble_70.tif", "Paepalanthus macrocephalus_ensemble_70.tif",
              "Paepalanthus scirpeus_ensemble_70.tif", "Paepalanthus nigrescens_ensemble_70.tif") 

south_70 <- stack(south_70)
south_final_70 <- mean(south_70)
south_70_sd <- calc(south_70, fun=sd)
writeRaster(south_final_70, "south_70.tif")
writeRaster(south_70_sd, "south_70_sd.tif")

north_70 <- c("Paepalanthus barbulatus_ensemble_70.tif", "Paepalanthus cinereus_ensemble_70.tif",
              "Paepalanthus harleyi_ensemble_70.tif", "Paepalanthus macrocaulon var. contasensis_ensemble_70.tif",
              "Paepalanthus pulvinatus_ensemble_70.tif", "Paepalanthus spathulatus_ensemble_70.tif")

north_70 <- stack(north_70)
north_final_70 <- mean(north_70)
north_70_sd <- calc(north_70, fun=sd)
writeRaster(north_final_70, "north_70.tif")
writeRaster(north_70_sd, "north_70_sd.tif")

#ploting the results for 2070
par(mfrow=c(2,2))
plot(north_final_70)
plot(north_70_sd)
plot(south_final_70)
plot(south_70_sd)

##For 2050
setwd(paste(wd, "distribution_models/1_future_50/", sep=""))

south_50 <- c("Paepalanthus bryoides_ensemble_50.tif", "Paepalanthus comans_ensemble_50.tif",
              "Paepalanthus eriophaeus_ensemble_50.tif", "Paepalanthus macrocephalus_ensemble_50.tif",
              "Paepalanthus scirpeus_ensemble_50.tif", "Paepalanthus nigrescens_ensemble_50.tif") 

south_50 <- stack(south_50)
south_final_50 <- mean(south_50)
south_50_sd <- calc(south_50, fun=sd)
writeRaster(south_final_50, "south_50.tif")
writeRaster(south_50_sd, "south_50_sd.tif")

north_50 <- c("Paepalanthus barbulatus_ensemble_50.tif", "Paepalanthus cinereus_ensemble_50.tif",
              "Paepalanthus harleyi_ensemble_50.tif", "Paepalanthus macrocaulon var. contasensis_ensemble_50.tif",
              "Paepalanthus pulvinatus_ensemble_50.tif", "Paepalanthus spathulatus_ensemble_50.tif")

north_50 <- stack(north_50)
north_final_50 <- mean(north_50)
north_50_sd <- calc(north_50, fun=sd)
writeRaster(north_final_50, "north_50.tif")
writeRaster(north_50_sd, "north_50_sd.tif")

#ploting the results for 2050
par(mfrow=c(2,2))
plot(north_final_50)
plot(north_50_sd)
plot(south_final_50)
plot(south_50_sd)


### Ploting all the models ###
plot(stack(north_final, north_final_50, north_final_70,
           south_final, south_final_50, south_final_70), zlim=c(0,0.8))
plot(stack(north_sd, north_50_sd, north_70_sd,
           south_sd, south_50_sd, south_70_sd), zlim=c(0,1))



#### End ####
