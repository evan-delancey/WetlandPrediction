###################################################################
#------------------------------------------------------------------
#------------------------------------------------------------------
#Boreal wetland probability
#Filename: "Boreal_WetlandProbability.R"
#Written and developed by Evan DeLancey - GIS Land Use Analyst
#Alberta Biodiversity Monitoring Institute, Nov, 14, 2017
#------------------------------------------------------------------
#------------------------------------------------------------------
###################################################################

#Load libraries
library(raster)
library(rgdal)
library(ggplot2)
library(dplyr)
library(caret)
library(snow)
library(rgeos)
library(RPyGeo)
library(dismo)
library(gbm)
library(ggthemes)

tt1 <- Sys.time()

#location of input rasters
location <- "J:/LandCover/ProbabilityofWetArea/Boreal/processed"
#location of land and wetland shapefiles
location.wetlandLand <- "J:/LandCover/ProbabilityofWetArea/Boreal/training"
#location of outputs
outputs <- "J:/LandCover/ProbabilityofWetArea/Boreal/OUTPUTS"
#version of the outputs
OutNum <- "_v1"
#enter number of iteration of subsampling
iter <- 40
#set your python location
py.loc <- "C:/Python27/ArcGIS10.4"
#Set minimum distance of sample wetland points
MinDist <- 1000

#set location of a temporary raster dump
#this can take up 100-300BG per run but is deleted after
rasterOptions(maxmemory = 1e+09,tmpdir = "J:/RtmpRasterDump")



#DEFINE FUNCTIONS
###################################################################
#------------------------------------------------------------------
#------------------------------------------------------------------

#1) set up function to get random points give location of your land and wetland shapefile and the minimum distance
pts.gen <- function(directory, mindist){
	#set up ArcGIS environment
	w <- directory
	env <- rpygeo.build.env(python.path = py.loc, workspace = directory, overwrite=1)
	rpygeo.geoprocessor("CreateRandomPoints_management", c("", "pts.shp", constraining_feature_class = "Boreal.shp", 
		constraining_extent = "", number_of_points_or_field = 80000, minimum_allowed_distance = mindist), env = env, detect.required.extensions = T)
	pts <- readOGR(w, "pts")
	return(pts)
}

#2)multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#----------------------------------------------------------------
#----------------------------------------------------------------
#################################################################



#Name vars and get min and max for response curves
###################################################################
#------------------------------------------------------------------
#------------------------------------------------------------------
#set input varibles
wetland.tifs <- c("ARI.tif", "NDPOL.tif", "NDWI.tif", "REIP.tif", "TPI.tif", "TWI.tif")
wetland.colnames <- c("ARI", "NDPOL", "NDWI", "REIP", "TPI", "TWI", "wetland")
wetland.bricknames <- c("ARI", "NDPOL", "NDWI", "REIP", "TPI", "TWI")
fls <- wetland.tifs

setwd(location)
#extract min and max of all input rasters for response curves
r <- raster(fls[1])
min1 <- cellStats(r, 'min')
max1 <- cellStats(r, 'max')
r <- raster(fls[2])
min2 <- cellStats(r, 'min')
max2 <- cellStats(r, 'max')
r <- raster(fls[3])
min3 <- cellStats(r, 'min')
max3 <- cellStats(r, 'max')
r <- raster(fls[4])
min4 <- cellStats(r, 'min')
max4 <- cellStats(r, 'max')
r <- raster(fls[5])
min5 <- cellStats(r, 'min')
max5 <- cellStats(r, 'max')
r <- raster(fls[6])
min6 <- cellStats(r, 'min')
max6 <- cellStats(r, 'max')
minval <- rbind(min1, min2, min3, min4, min5, min6)
maxval <- rbind(max1, max2, max3, max4, max5, max6)
mm.df <- data.frame(wetland.bricknames, minval, maxval)

#binary training raster

wetland <- raster("J:/LandCover/ProbabilityofWetArea/Boreal/training/wetland.tif")
#----------------------------------------------------------------
#----------------------------------------------------------------
#################################################################



#BUILD MODEL 1
###################################################################
#------------------------------------------------------------------
#------------------------------------------------------------------
#define fit AUC and dev
AUC <- vector()
dev <- vector()
pts.model <- pts.gen(location.wetlandLand, MinDist)
dat <- data.frame(row=1:length(pts.model))
for (i in 1:length(fls)){
	r <- raster(fls[i])
	ext <- extract(r,pts.model)
	dat <- cbind(dat,ext)
	print(paste0("done ",fls[i]))
}
ext <- extract(wetland, pts.model)
dat <- cbind(dat,ext)
dat <- dat[,-1]
colnames(dat) <- wetland.colnames
dat <- na.omit(dat)
cor(dat)

#defint model list
fit <- list()
#build model
fit[[1]] <- gbm.step(dat, 1:length(fls), length(fls)+1, family = "bernoulli", tree.complexity = 5,
	learning.rate = 0.005, bag.fraction = 0.5, silent = TRUE, warnings = FALSE)
df.importance <- data.frame(summary(fit[[1]]))
df.importance <- arrange(df.importance, var)
df.importance <- df.importance[,2]
#model stats
AUC[1] <- fit[[1]]$cv.statistics$discrimination.mean
dev[1] <- (fit[[1]]$self.statistics$mean.null - fit$self.statistics$mean.resid) / fit$self.statistics$mean.null
	
response.df <- data.frame(dummy=c(1:1001))
for (n in wetland.bricknames){
	d <- plot.gbm(fit[[1]], i.var = n, return.grid=TRUE, type="response")
	get.min.max <- filter(mm.df, wetland.bricknames == n)
	mn <- get.min.max[,2]
	mx <- get.min.max[,3]
	xout <- seq(mn, mx, (mx-mn)/1000)
	d <- approx(d[,1], d[,2], xout = xout, rule=2)
	d <- as.data.frame(d)
	response.df <- cbind(response.df,d)
}
#----------------------------------------------------------------
#----------------------------------------------------------------
#################################################################






#BUILD MODEL ALL MODELS AND OUTPUT STATS
###################################################################
#------------------------------------------------------------------
#------------------------------------------------------------------
setwd(location)


for (i in 2:iter){
	pts.model <- pts.gen(location.wetlandLand, MinDist)
	dat <- data.frame(row=1:length(pts.model))
	for (j in 1:length(fls)){
		r <- raster(fls[j])
		ext <- extract(r,pts.model)
		dat <- cbind(dat,ext)
	}
	ext <- extract(wetland, pts.model)
	dat <- cbind(dat,ext)
	dat <- dat[,-1]
	colnames(dat) <- wetland.colnames
	dat <- na.omit(dat)
	
	#build model
	fit[[i]] <- gbm.step(dat, 1:length(fls), length(fls)+1, family = "bernoulli", tree.complexity = 5,
		learning.rate = 0.005, bag.fraction = 0.5, silent = TRUE, warnings = FALSE)
	v.importance <- as.data.frame(summary(fit[[i]]))
	v.importance <- arrange(v.importance, var)
	v.importance <- v.importance[,2]
	df.importance <- cbind(df.importance, v.importance)
	
	for (n in wetland.bricknames){
		d <- plot.gbm(fit[[i]], i.var = n, return.grid=TRUE, type="response")
		get.min.max <- filter(mm.df, wetland.bricknames == n)
		mn <- get.min.max[,2]
		mx <- get.min.max[,3]
		xout <- seq(mn, mx, (mx-mn)/1000)
		d <- approx(d[,1], d[,2], xout = xout, rule=2)
		d <- as.data.frame(d)
		response.df <- cbind(response.df,d)
	}
	
	#model stats
	AUC[i] <- fit[[i]]$cv.statistics$discrimination.mean
	dev[i] <- (fit[[i]]$self.statistics$mean.null - fit$self.statistics$mean.resid) / fit$self.statistics$mean.null
	
	print(paste0("done building model ", i))
}

importance <- rowMeans(df.importance)
imp.names <- c("ARI", "NDPOL", "NDWI", "REIP", "TPI", "TWI")
importance <- data.frame(imp.names, importance)
#----------------------------------------------------------------
#----------------------------------------------------------------
#################################################################





#GENERATE RESPONSE CURVES AND MODEL STATS
###################################################################
#------------------------------------------------------------------
#------------------------------------------------------------------
response.df <- response.df[,-1]
response.df.x <- response.df[,seq(1,length(response.df),2)] 
response.df.x <- response.df.x[,1:length(wetland.bricknames)]
response.df.y <- response.df[,seq(2,length(response.df),2)] 
for (i in 1:length(wetland.bricknames)){
	varname <- wetland.bricknames[i]
	collumns <- seq(i, length(response.df.y), length(wetland.bricknames)) 
	yvals <- response.df.y[,collumns]
	xvals <- response.df.x[,i]
	yvals.mean <- rowMeans(yvals)
	yvals.std <- apply(yvals,1,sd)
	yvals.neg.std <- yvals.mean - yvals.std
	yvals.add.std <- yvals.mean + yvals.std
	yvals.df <- cbind(xvals, yvals.mean, yvals.neg.std, yvals.add.std)
	yvals.df <- as.data.frame(yvals.df)
	if(i>0){
		xlim <- c(min(xvals), max(xvals))
	} else{
		xlim <- c(min(xvals), 1200)
	}
	g <- ggplot(yvals.df, aes(x=xvals, y=yvals.mean)) + 
		theme_minimal()+ 
		geom_ribbon(aes(ymin=yvals.neg.std, ymax=yvals.add.std), fill="#6baed6", alpha=0.35) +
		geom_line(colour="#08519c", size=1.8) +
		xlab(varname) + ylab("predicted probability") + xlim(xlim) +
		theme(axis.title.x = element_text(size=22), axis.title.y = element_text(size=20),axis.text = element_text(size=16))
	assign(paste0(varname,".plot"), g)
}

#output model stats
setwd(outputs)
AUC <- mean(AUC)
dev <- mean(dev)
stats <- cbind(AUC,dev)
write.csv(stats, paste0("wetlandModelStats", ".csv"), row.names=FALSE)

tiff(paste0("wetlandResponseCurves", ".tiff"), width = 1200, height = 1000)
multiplot(ARI.plot, NDPOL.plot, NDWI.plot, REIP.plot, TPI.plot, TWI.plot, cols=2)
dev.off()

tiff(paste0("wetlandVarImportance", ".tiff"), width = 1000, height = 700)
ggplot(importance,aes(x=reorder(imp.names,-importance),y=importance)) + geom_bar(stat="identity", show.legend=FALSE, fill= "grey60") + theme_minimal() +
	theme(axis.title.x = element_text(size=24), axis.title.y = element_text(size=24),axis.text = element_text(size=22), legend.text = element_text(size=20), legend.title = element_text(size=22)) + 
	labs(x = "Variables") + labs(y = "Importance" ) 
dev.off()
#----------------------------------------------------------------
#----------------------------------------------------------------
#################################################################



###################################################################
#SECTION 2 PREDICT WETLAND PROBABILITY BY TILE
###################################################################

#LOOP THROUGH PREDICTION OF RASTERS BASED ON MODEL FITS
###################################################################
#------------------------------------------------------------------
#------------------------------------------------------------------

PUs <- readOGR("J:/LandCover/CurrentSurfaceWater/Boreal", "Boreal_PUs")	
for (i in 1:length(PUs)){

	t1 <- Sys.time()
	
	#build raster brick
	setwd(location)
	fls <- wetland.tifs
	PU <- PUs[i,]
	#build raster brick
	r1 <- raster(fls[1])
	r1 <- crop(r1, PU)
	print("done cropping1")

	r2 <- raster(fls[2])
	r2 <- crop(r2, PU)

	r3 <- raster(fls[3])
	r3 <- crop(r3, PU)

	r4 <- raster(fls[4])
	r4 <- crop(r4, PU)
	
	r5 <- raster(fls[5])
	r5 <- crop(r5, PU)
	
	r6 <- raster(fls[6])
	r6 <- crop(r6, PU)

	
	r.b <- brick(r1,r2,r3,r4, r5, r6)
	names(r.b) <- wetland.bricknames
	
	#predict fit[[1]]
	beginCluster(20)
	r.b.p <- clusterR(r.b, raster::predict, args = list(model = fit[[1]], type = "response", n.trees = fit[[1]]$gbm.call$best.trees))
	endCluster()
	plot(r.b.p)
	
	#START PREDICTION OF RASTER STACK
	for (n in 2:length(fit)){
		beginCluster(20)
		prediction <- clusterR(r.b, raster::predict, args = list(model = fit[[n]], type = "response", n.trees = fit[[n]]$gbm.call$best.trees))
		endCluster()
		r.b.p <- stack(r.b.p, prediction)
		print(paste0("Done model prediction ", n))
	}
	
	wet <- calc(r.b.p, fun = mean)
	
	#Save wetland prediction rasters
	setwd(outputs)
	writeRaster(wet, paste0("WetlandProbability", i, ".tif"), datatype= "FLT4S")
	
	t2 <- Sys.time()
	t.diff <- difftime(t2,t1, units="hours")
	print(paste0("Done predicting PU ", i, " it took ", round(t.diff, 2), " hours"))
	
}
#------------------------------------------------------------------
#------------------------------------------------------------------
####################################################################
tt2 <- Sys.time()
t.diff <- difftime(tt2,tt1, units="hours")
print(paste0("total script took ", round(t.diff, 2), " hours"))
 







