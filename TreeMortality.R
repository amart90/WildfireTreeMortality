##########  Wildfire tree mortality detection and  ##########
##########   prediction using LiDAR and Landsat    ##########
##########         Senior capstone project         ##########
##########             Anthony Martinez            ##########
##########                June 2017                ##########


setwd("/media/antho/HDD/KingFire/CapstoneU")
#set up to run LiDAR data using parallel processing. Must be using Linux system.

######  Post-fire burn mortality classification  ######

#load packages
library(sp)
library(raster)
library(rgdal)
library(foreach)
library(iterators)
library(parallel)
library(doMC)
library(randomForest)
library(rgeos)

#identify number of cores for parallel processing
registerDoMC(cores = 6)

#load Landsat rasters

preB5i = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/Initial/LC80430332013291LGN00-PRE/LC80430332013291LGN00_B5.TIF')
preB7i = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/Initial/LC80430332013291LGN00-PRE/LC80430332013291LGN00_B7.TIF')

postB1i = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/Initial/LC80430332014294LGN00-POST/LC80430332014294LGN00_B1.TIF')
postB2i = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/Initial/LC80430332014294LGN00-POST/LC80430332014294LGN00_B2.TIF')
postB3i = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/Initial/LC80430332014294LGN00-POST/LC80430332014294LGN00_B3.TIF')
postB4i = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/Initial/LC80430332014294LGN00-POST/LC80430332014294LGN00_B4.TIF')
postB5i = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/Initial/LC80430332014294LGN00-POST/LC80430332014294LGN00_B5.TIF')
postB6i = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/Initial/LC80430332014294LGN00-POST/LC80430332014294LGN00_B6.TIF')
postB7i = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/Initial/LC80430332014294LGN00-POST/LC80430332014294LGN00_B7.TIF')

preB5y = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/1YrPost/LC80430332014246LGN00-PRE/LC80430332014246LGN00_B5.TIF')
preB7y = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/1YrPost/LC80430332014246LGN00-PRE/LC80430332014246LGN00_B7.TIF')

postB1y = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/1YrPost/LC80430332015249LGN00-POST/LC80430332015249LGN00_B1.TIF')
postB2y = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/1YrPost/LC80430332015249LGN00-POST/LC80430332015249LGN00_B2.TIF')
postB3y = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/1YrPost/LC80430332015249LGN00-POST/LC80430332015249LGN00_B3.TIF')
postB4y = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/1YrPost/LC80430332015249LGN00-POST/LC80430332015249LGN00_B4.TIF')
postB5y = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/1YrPost/LC80430332015249LGN00-POST/LC80430332015249LGN00_B5.TIF')
postB6y = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/1YrPost/LC80430332015249LGN00-POST/LC80430332015249LGN00_B6.TIF')
postB7y = raster('/media/antho/HDD/KingFire/Capstone/Landsat8/1YrPost/LC80430332015249LGN00-POST/LC80430332015249LGN00_B7.TIF')
files <- ls()

#clip rasters
KingBuffer <- readOGR('/media/antho/HDD/KingFire/Capstone/R/Capstone/KingBuffer.shp')
x <- foreach(i=1:length(files), .export = files) %dopar% {
  assign(files[[i]], crop(get(files[[i]]), extent(KingBuffer)))
  assign(files[[i]], mask(get(files[[i]]), KingBuffer))}
for(i in 1:length(files)) {assign(files[[i]], x[[i]])}

#NBR, initial
preNBRi <- (1000 * ((preB5i - preB7i) / (preB5i + preB7i)))
postNBRi <- (1000 * ((postB5i - postB7i) / (postB5i + postB7i)))

#NBR, 1 year later
preNBRy <- (1000 * ((preB5y - preB7y) / (preB5y + preB7y)))
postNBRy <- (1000 * ((postB5y - postB7y) / (postB5y + postB7y)))

##dNBR (corr)
##mode calculation i
#modei1 <- (preNBRi - postNBRi)
#ddi <- density(modei1)
#which.max(ddi$y)
#modei <- ddi$x[which.max(ddi$y)]
##mode calculation y
#modey1 <- (preNBRy - postNBRy)
#ddy <- density(modey1)
#which.max(ddy$y)
#modey <- ddy$x[which.max(ddy$y)]

##raster calculations

#dNBRicorr <- ((preNBRi - postNBRi) - modei)
#dNBRycorr <- ((preNBRy - postNBRy) - modey)

#dNBR
dNBRi <- preNBRi - postNBRi
dNBRy <- preNBRy - postNBRy

#RbNBR
RdNBRi <- (dNBRi / (sqrt(abs(preNBRi / 1000))))
RdNBRy <- (dNBRy / (sqrt(abs(preNBRi / 1000))))

#NDMI
NDMIi <- (postB5i - postB6i) / (postB5i + postB6i)
NDMIy <- (postB5y - postB6y) / (postB5y + postB6y)

#B5B4 (B6/B5 in Landsat 8)
B5B4i <- postB6i / postB5i
B5B4y <- postB6y / postB5yvaluetable <- na.omit(getValues(trainingbrick))

#NDVI
NDVIi <- (postB5i - postB4i) / (postB5i + postB4i)
NDVIy <- (postB5y - postB4y) / (postB5y + postB4y)

#TCBRI
TCBRIi <- 0.2049 * postB2i + 0.4158 * postB3i + 0.5524 * postB4i + 0.5741 * postB5i + 0.3124 * postB6i + 0.2303 * postB7i
TCBRIy <- 0.2049 * postB2y + 0.4158 * postB3y + 0.5524 * postB4y + 0.5741 * postB5y + 0.3124 * postB6y + 0.2303 * postB7y

#TCGRE
TCGREi <- -0.1603 * postB2i - 0.2819 * postB3i - 0.4934 * postB4i + 0.79401 * postB5i - 0.0002 * postB6i - 0.1446 * postB7i
TCGREy <- -0.1603 * postB2y - 0.2819 * postB3y - 0.4934 * postB4y + 0.79401 * postB5y - 0.0002 * postB6y - 0.1446 * postB7y

#TCWET
TCWETi <-0.0315 * postB2i + 0.2021 * postB3i + 0.3102 * postB4i + 0.1594 * postB5i - 0.6806 * postB6i - 0.6109 * postB7i
TCWETy <-0.0315 * postB2y + 0.2021 * postB3y + 0.3102 * postB4y + 0.1594* postB5y - 0.6806 * postB6y - 0.6109 * postB7y

#write rasters
rasternames <- c("RdNBRi", "RdNBRy", "dNBRi", "dNBRy", "NDMIi", "NDMIy", "B5B4i", "B5B4y", "NDVIi", "NDVIy", "TCBRIi", "TCBRIy", "TCGREi", "TCGREy", "TCWETi", "TCWETy")
rasterlist <- mget(rasternames)
covs <- stack(rasterlist)

#remove intermediate files
remove <- c(files, rasternames, "postNBRi", "postNBRy", "preNBRi", "preNBRy")
rm(list = remove)
gc()

#Create Random points within fire boundary (specify number of points)
library(splancs)
rndmpts <- spsample(KingOutline, n = 300, type = "random")
projection(rndmpts)  <- projection(KingOutline)
rndmptsdf <- data.frame(FID = row.names(rndmpts), X100pc = rep(NA, times = length(rndmpts)), row.names = row.names(rndmpts))
rndmptsdf <- SpatialPointsDataFrame(rndmpts, data = df, proj4string = project(KingOutline))
writeOGR(obj = rndmptsdf, dsn = getwd(), layer = "rndmpts", driver = "ESRI Shapefile", overwrite_layer = T)

#Stop to Identify with "Ground truthing" methods
stop("Visually identify areas of 100% mortality at point locations at this time")

#import points
#Spoints1 = readOGR('rndmpts.shp')
#Spoints2 = readOGR('/media/antho/HDD/KingFire/Capstone/R/Capstone/Burned100pc.shp')
#Spoints <- readOGR('Newrndmpts.shp')
proj4string(Spoints) <- proj4string(covs)
plot(Spoints1)
plot(Spoints2, add = T)


#rasterize points
SRaster <- rasterize(Spoints, covs, field = 'X100pc')

#mask data to points
covmask <- mask(covs, SRaster)

#combine covariates with sample data to make training dataset
names(SRaster) <- "Burned100"
trainingbrick <- addLayer(covmask, SRaster)

#Extract values to matrix
valuetable <- na.omit(getValues(trainingbrick))
valuetable <- data.frame(valuetable)
valuetable$Burned100 <- as.factor(valuetable$Burned100)

#random forest model
rfm1 <- randomForest(Burned100 ~., data = valuetable, importance = T)
rfm1.error <- mean(rfm1$err.rate)
rfm1i <- rfm1
rfm1i$importance[,"MeanDecreaseAccuracy"] <- rfm1i$importance[,"MeanDecreaseAccuracy"]*100
varImpPlot(rfm1i, main = "Burn Severity Classification Model", type = 1, scale=F)
plot(rfm2)
#use rfm1 to predict
predBurn <- predict(covs, model = rfm1, filename = "predBurn.asc", format = "ascii", na.rm = TRUE, inf.rm = TRUE, overwrite = TRUE)
proj4string(predBurn) <- proj4string(covs)
predBurn <- round(predBurn, digits = 0)
plot(predBurn)

#plot proportion of classified points
h <- hist(predBurn, breaks=2)
c.t <- sum(h$counts) #total cells
c0 <- h$counts[1] #count of 0s
c1 <- h$counts[2] #count of 1s
p0 <- c0/c.t #proportion of 0s
p1 <- c1/c.t #proportion of 1s
hist(predBurn, breaks = 2, col = c("#5C5C76", "#CD6D3A"), freq=T, 
     xaxp = c(0,1,1), main = "Frequency of Classified Burn Status", xlab = "Predicted Value")
legend("topright", legend = c("< 100% Tree Mortality", "100% Tree Mortality"), 
       fill = c("#5C5C76", "#CD6D3A"), title = "Classifications")
text(.25, c1 / 2, labels = paste(round(p0*100, digits = 1), "%"), col = "white", cex=1.5)
text(.75, c1 / 2, labels = paste(round(p1*100, digits = 1), "%"), col = "white", cex=1.5)

######  Prediction of Burn Severity Using LiDAR  ######

##Setup
#load libraries
library(sp)
library(raster)
library(rgdal)
library(lidR)
library(foreach)
library(doParallel)
library(randomForest)
library(GISTools)

#set working directory
setwd("D:/KingFire/Capstone/R/Capstone")

##Input spatial files
#fire perimeter
KingOutline <- readOGR('D:/KingFire/Capstone/R/Capstone/KingOutline.shp')

#bare earth raster
ground <- raster("D:/KingFire/CoverStrata/BareEarthKing.asc")

#100% burned raster (output from random forest #1)
predBurn <- raster('predBurn.asc')

#LAS delivery tiles
DelivTiles <- readOGR('D:/KingFire/Capstone/R/Capstone/LAS_DeliveryTiles.shp')

##List intersecting LAS files
KingOutlineProj <- spTransform(KingOutline, projection(DelivTiles))
plot(KingOutlineProj)
plot(DelivTiles, add = TRUE)
DelivInter <- raster::intersect(DelivTiles, KingOutlineProj)
plot(DelivInter)
Ipf <- as.character(DelivInter@data$Filename)
InputFiles <- paste0("D:/KingFire/PreBurn/Raw/Points/FullCloud/", Ipf)
InputFiles <- sort(InputFiles)
InputFileNames <- gsub(".las", "", gsub("D:/KingFire/PreBurn/Raw/Points/FullCloud/", "", InputFiles, fixed = TRUE))

##Import LAS files, normalize them, perform grid statistics, and remove (large) intermediate files
#abridge list
InputFiles <- InputFiles[1:2]
InputFileNames <- InputFileNames[1:2]

ext <- c(xmax(predBurn), ymax(predBurn))
LAS_Grid <- data.frame()

no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)
registerDoParallel(cl)
stime2 <- system.time(x <- foreach(i = 1:length(InputFileNames)) %dopar% {
  library(lidR)
  library(raster)
  fl <- paste0("D:/KingFire/PreBurn/Raw/Points/FullCloud/", InputFileNames[[i]], ".las")
  LASfile <- paste0("LAS_", InputFileNames[[i]])
  assign(LASfile, readLAS(fl))
  LASnorm <- paste0("LAS_", InputFileNames[[i]], "_norm")
  assign(LASnorm, lasnormalize(get(LASfile), ground))
  get(LASnorm) %>% lasfilter(Z > 2) %>% grid_metrics(.stdmetrics_z, res = 30, start = ext)
  remove <- c(as.character(LASfile), as.character(LASnorm))
  rm(list = remove)})

stopCluster(cl)

for(i in 1:length(x)) {
  LAS_Grid <- rbind(LAS_Grid, x[[i]])}

##Produce rasters
StdMetrics <- colnames(LAS_Grid[,-c(1:2)])
LASrNames <-paste0("LAS_Raster_", StdMetrics)
LAS_Grid_sp <- LAS_Grid
coordinates(LAS_Grid_sp) <- ~ X + Y
gridded(LAS_Grid_sp) <- TRUE
proj4string(LAS_Grid_sp) <- proj4string(DelivTiles)
for(i in StdMetrics) {
  name <- paste0("LAS_Raster_", i)
  assign(name, raster(LAS_Grid_sp, layer = i))
  plot(get(name), main = i)}

covs2 <- list()
for(i in 1:length(LASrNames)) {
  name <- LASrNames[[i]]
  covs2[[which(LASrNames %in% name)]] <- get(name)}
covs2 <- stack(covs2)
predBurn <- raster('predBurn.asc')
predBurn <- projectRaster(predBurn, crs = proj4string(KingOutline))
covs2 <- projectRaster(covs2, crs = proj4string(KingOutline))
predmask <- mask(predBurn, KingOutline)
covmask2 <- mask(covs2, KingOutline)
predmask <- crop(predmask, covmask2)
trainingbrick2 <- addLayer(covmask2, predmask)

#Extract values to matrix
valuetable2 <- na.omit(getValues(trainingbrick2))
valuetable2 <- data.frame(valuetable2)
valuetable2$predBurn <- as.factor(round(valuetable2$predBurn))

#random forest model
rfm2 <- randomForest(predBurn ~., data = valuetable2, importance = T)
varImpPlot(rfm2)
rfm2.error <- mean(rfm2$err.rate)
rfm2$err.rate
str(rfm2)

#use rfm2 to predict
predStrctr <- predict(covmask2, model = rfm2, filename = "predStrctr.asc", format = "ascii", na.rm = TRUE, overwrite = TRUE)
predStrctr.prob <- predict(covmask2, model = rfm2, type = "prob")
plot(predStrctr.prob)
plot(predStrctr)
plot(predmask)

#
predBurn.str.mask <- mask(crop(predBurn, predStrctr), predStrctr)
plot(predBurn.str.mask)
plot(predStrctr.prob, add = T)
mean(predStrctr.prob[predBurn.str.mask == 1], na.rm = T)
mean(predStrctr.prob[predBurn.str.mask == 0], na.rm = T)
predStrctr.error <- predStrctr != predBurn.str.mask
predStrctr.error2 <- predStrctr == predBurn.str.mask
plot(predStrctr.error)

cellStats(predStrctr.error, stat = 'sum', na.rm = T)
cellStats(predStrctr.error2, stat = 'sum', na.rm = T)
rfm2$confusion

#rfm3 most important
imp.list <- c("predBurn","zpcum1","zpcum2","zpcum3","zpcum4","zpcum5","zpcum6","zpcum7","zpcum8","zpcum9","zkurt","zskew")
valuetable3 <- valuetable2[imp.list]
rfm3 <- randomForest(predBurn ~., data= valuetable3, importance = T)
varImpPlot(rfm3)
rfm3$confusion
str3 <- predict(covmask2, rfm3)
str3p <- predict(covmask2, rfm3, type = "prob")
plot(str3, main = "predStruc")
srt.9 <- str3p
srt.9[srt.9 > 0.1] <- 1
plot(srt.9, main = "str.1")
plot(predBurn.str.mask, main = "predBurn")

######  Plot maps for prediction model  ######
library(GISTools)

# Distance raster
f1 <- round(predBurn)
f2 <- mask(gridDistance(f1, 0), KingOutline)
breakpoints <- c(0, 1, 200, 400, 800, maxValue(f2))
colors <- c("#5C5C76", "#F9EB38", "#F5BF1A", "#CD6D3A", "#9F4837")
plot(f2, breaks = breakpoints, col = colors, axes = FALSE, box = OFF, legend = FALSE)
writeRaster(f2, "KingStandReplacement.asc", "ascii", overwrite = TRUE)


#distance PDF
brk <- c("Surviving Forest Canopy", "< 200 m", "200 - 400 m", "400 - 800 m", "> 800 m")
par(mar=c(1,1,1,1))
pdf("KingStandReplacement.pdf", width = 7, height = 10, paper = "letter")
plot(f2, breaks = breakpoints, col = colors, axes = FALSE, box = FALSE, legend = FALSE, main = "King Fire - Stand Replacement Area", add = FALSE)
legend("topleft", legend = brk, fill = colors, title = "Distance to Nearest Live-tree Edge", cex = 0.9)
plot(KingOutline, add = TRUE)
north.arrow(xmin(f2) + 7000, ymin(f2) + 37000, 500, lab = "N")
map.scale(xmin(f2) + 7000, ymin(f2) + 35000, 10 * 1000, "Km", 5, 2)
dev.off()

#distance PDF
f3 <- round(predBurn)
f3 <- mask(f3, KingOutline)
brk2 <- c("100% Tree Mortality", "< 100% Tree Mortality")
par(mar=c(1,1,1,1))
pdf("BurnPrediction.pdf", width = 7, height = 10, paper = "letter")
plot(f3, breaks = c(0,.9,1.5), col = c("#5C5C76", "#CD6D3A"), axes = FALSE, box = FALSE, legend = FALSE, main = "King Fire - Burn Mortality Status", add = FALSE)
legend("topleft", legend = brk2, fill = c("#CD6D3A", "#5C5C76"), title = "Predicted Burn Status", cex = 0.9)
plot(KingOutline, add = TRUE)
north.arrow(xmin(f2) + 7000, ymin(f2) + 37000, 500, lab = "N")
map.scale(xmin(f2) + 7000, ymin(f2) + 35000, 10 * 1000, "Km", 5, 2)
dev.off()

#distance jpeg
par(mar=c(1,1,1,1))
jpeg("KingStandReplacement.jpeg", width = 7, height = 10, units = "in", res = 300, family = "", type = "windows")
plot(f2, breaks = breakpoints, col = colors, axes = FALSE, box = FALSE, legend = FALSE, main = "King Fire - Stand replacement area", add = FALSE)
legend("topleft", legend = brk, fill = colors, title = "Distance to Nearest Live-tree Edge", cex = 0.9)
plot(KingOutline, add = TRUE)
north.arrow(xmin(f2) + 7000, ymin(f2) + 37000, 500, lab = "N")
map.scale(xmin(f2) + 7000, ymin(f2) + 35000, 10 * 1000, "Km", 5, 2)
dev.off()

#burn prediction jpeg
par(mar=c(1,1,1,1))
jpeg("BurnPrediction.jpeg", width = 7, height = 10, units = "in", res = 300, family = "", type = "windows")
plot(f3, breaks = c(0,.9,1.5), col = c("#5C5C76", "#CD6D3A"), axes = FALSE, box = FALSE, legend = FALSE, main = "King Fire - Burn Mortality Status", add = FALSE)
legend("topleft", legend = brk2, fill = c("#CD6D3A", "#5C5C76"), title = "Predicted Burn Status", cex = 0.9)
plot(KingOutline, add = TRUE)
north.arrow(xmin(f2) + 7000, ymin(f2) + 37000, 500, lab = "N")
map.scale(xmin(f2) + 7000, ymin(f2) + 35000, 10 * 1000, "Km", 5, 2)
dev.off()

#LiDAR Prediction jpeg
par(mar=c(1,1,1,1))
jpeg("LiDARPrediction.jpeg", width = 7, height = 10, units = "in", res = 300, family = "", type = "windows")
plot(predBurn.str.mask, breaks = c(0,.9,1.5), col = c("#5C5C76", "#CD6D3A"), axes = FALSE, box = FALSE, legend = FALSE, main = "King Fire - Prediction of Burn severity using LiDAR", add = FALSE)
legend("topleft", legend = brk2, fill = c("#CD6D3A", "#5C5C76"), title = "Predicted Burn Status", cex = 0.9)
plot(KingOutline, add = TRUE)
north.arrow(xmin(predBurn.str.mask) + 1600, ymin(predBurn.str.mask) + 16000, 500, lab = "N")
map.scale(xmin(predBurn.str.mask) + 1600, ymin(predBurn.str.mask) + 15000, 6 * 1000, "Km", 3, 2)
dev.off()
