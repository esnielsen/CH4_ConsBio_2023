#### R scripts for "Integrating environmental, evolutionary, and socioeconomic vulnerability to future-proof coastal conservation planning"
### written by Erica Nielsen (esnielsen@ucdavis.edu)
### uses data from Nielsen et al. 2021, Glob Change Biol & Nielsen et al.2021 J Biogeo 
### acronyms: CP=C. punctatus, PA=P.angulosus, SG=S.granularis

### not all of the biodiversity spatial layers are available.. will have to generate: climatic stability (termed 'similarity' in publication), climatic variability, human footprint, population change

######################################################
##### Part 1: Calculating climatic variables #####
######################################################

library(climateStability)
library(raster)
library(sdmpredictors)
library("xlsx")
library("readxl")
library(gstat)
library(raster)
library(tidyverse)
library(prioritizr)
library(rgdal)
library(rgeos)
library(mapview)
library(units)
library(scales)
library(assertthat)
library(gridExtra)
library(dplyr)
library(GGally)
library(fmsb)

#####  Generating climatic stability #####

##### DOWNLOAD 
##### CURRENT ENV

#Tmean <- load_layers("WC_bio1")
SSTmean <- load_layers("BO22_tempmean_ss")
SSSmean <- load_layers("BO22_salinitymean_ss")
#Pmean <- load_layers("WC_bio12")
#CHLmean <- load_layers("BO2_chlomean_ss")

##### DOWNLOAD 
##### FUTURE ENV

#Tmean.45.t <- load_layers("WC_bio1_mr45_2050")
SSTmean.45.t <- load_layers("BO2_RCP45_2050_tempmean_ss")
SSSmean.45.t <- load_layers("BO2_RCP45_2050_salinitymean_ss")
#Pmean.45.t <- load_layers("WC_bio12_mr45_2050")
#CHLmean.45.t <- load_layers("BO2_RCP45_2050_chlomean_ss")

#Tmean.85.t <- load_layers("WC_bio1_mr85_2050")
SSTmean.85.t <- load_layers("BO2_RCP85_2050_tempmean_ss")
SSSmean.85.t <- load_layers("BO2_RCP85_2050_salinitymean_ss")
#Pmean.85.t <- load_layers("WC_bio12_mr85_2050")
#CHLmean.85.t <- load_layers("BO2_RCP85_2050_chlomean_ss")

# Set species extents
SA.ext <- extent(14, 34, -35.5, -25)
#SA.ext <- extent(5, 45, -40, -10)
#SA.cp.ext <- extent(10, 40, -37, -22)
#SA.pa.ext <- extent(10, 40, -37, -19)
#SA.sg.ext <- extent(10, 40, -37, -12)

#Tmeanc <- crop(Tmean, SA.ext)
SSTmeanc <- crop(SSTmean, SA.ext)
SSSmeanc <- crop(SSSmean, SA.ext)
#Pmeanc <- crop(Pmean, SA.ext)
#CHLmeanc <- crop(CHLmean, SA.ext)

#Tmean.45.tc <- crop(Tmean.45.t, SA.ext)
SSTmean.45.tc <- crop(SSTmean.45.t, SA.ext)
SSSmean.45.tc <- crop(SSSmean.45.t, SA.ext)
#Pmean.45.tc <- crop(Pmean.45.t, SA.ext)
#CHLmean.45.tc <- crop(CHLmean.45.t, SA.ext)

#Tmean.85.tc <- crop(Tmean.85.t, SA.ext)
SSTmean.85.tc <- crop(SSTmean.85.t, SA.ext)
SSSmean.85.tc <- crop(SSSmean.85.t, SA.ext)
#Pmean.85.tc <- crop(Pmean.85.t, SA.ext)
#CHLmean.85.tc <- crop(CHLmean.85.t, SA.ext)

#Tmean.50 <- overlay(Tmean.45.tc, Tmean.85.tc, fun=mean)
SSTmean.50 <- overlay(SSTmean.45.tc, SSTmean.85.tc, fun=mean)
SSSmean.50 <- overlay(SSSmean.45.tc, SSSmean.85.tc, fun=mean)
#Pmean.50 <- overlay(Pmean.45.tc, Pmean.85.tc, fun=mean)
#CHLmean.50 <- overlay(CHLmean.45.tc, CHLmean.85.tc, fun=mean)

#writeRaster(Tmean.50, "Tmean.50.tif", format="GTiff", overwrite=TRUE)
writeRaster(SSTmean.50, "SSTmean.50.tif", format="GTiff", overwrite=TRUE)
writeRaster(SSSmean.50, "SSSmean.50.tif", format="GTiff", overwrite=TRUE)
#writeRaster(Pmean.50, "Pmean.50.tif", format="GTiff", overwrite=TRUE)
#writeRaster(CHLmean.50, "CHLmean.50.tif", format="GTiff", overwrite=TRUE)

#writeRaster(Tmeanc, "Tmean.c.tif", format="GTiff", overwrite=TRUE)
writeRaster(SSTmeanc, "SSTmean.c.tif", format="GTiff", overwrite=TRUE)
writeRaster(SSSmeanc, "SSSmean.c.tif", format="GTiff", overwrite=TRUE)
#writeRaster(Pmeanc, "Pmean.c.tif", format="GTiff", overwrite=TRUE)
#writeRaster(CHLmeanc, "CHLmean.c.tif", format="GTiff", overwrite=TRUE)

#Tdiff <- Tmean.50 - Tmeanc
SSTdiff <- SSTmean.50 - SSTmeanc
SSSdiff <- SSSmean.50 - SSSmeanc
#Pdiff <- Pmean.50 - Pmeanc
#CHLdiff <- CHLmean.50 - CHLmeanc

#T.d.norm <- rescale0to1(Tdiff) #rescale from the package climateStability
SST.d.norm <- rescale0to1(SSTdiff)
SSS.d.norm <- rescale0to1(SSSdiff)
#P.d.norm <- rescale0to1(Pdiff)
#CHL.d.norm <- rescale0to1(CHLdiff)

# this is sum of normalized difference between 2010 and 2050 for each clim var
All.mar.diff <- SST.d.norm + SSS.d.norm 

setwd("~/Desktop/CH4_ms/input_data")
writeRaster(All.mar.diff, "All.mar.diff.tif", format="GTiff", overwrite=TRUE)


####### PRESENT DAY CLIM VARIABLILTY #########


setwd("~/Desktop/CH4_ms/input_data/SDMs")

SSTvar <- load_layers("BO22_temprange_ss")
SSSvar <- load_layers("BO22_salinityrange_ss")

SSTvar.c <- crop(SSTvar, SA.ext)
SSSvar.c <- crop(SSSvar, SA.ext)

SST.v.norm <- rescale0to1(SSTvar.c)
SSS.v.norm <- rescale0to1(SSSvar.c)

All.mar.var <- SST.v.norm + SSS.v.norm

writeRaster(All.mar.var, "All.mar.var.tif", format="GTiff", overwrite=TRUE)

######################################################
##### Part 2: Prepping layers for prioritization #####
######################################################

######                      ######
###### Calculate Thresholds ######


### Human Footprint (FP) ###
setwd("~/Desktop/CH4_ms/shps")
HP <- raster::stack("HFP.reproj.tif")

#HPc <- crop(HP, SA.ext)

pu_ext <- extent(16.48, 32.98, -34.77, -26.87)
HPc <- crop(HP, pu_ext)
HPm <- mask(HPc, pu_data)

threshold <-
  raster::quantile(HPm,
                   probs = 0.25,
                   na.rm = TRUE)

### clamp data to threshold such that values ABOVE threshold become
### zero and other values remain the same
HPm[HPm > threshold] <- 0

### set it so areas with greater than zero are all 1's
HPm[HPm > 0] <- 1

#write raster
writeRaster(HPm, "all.new.HFP.tif", format="GTiff", overwrite=TRUE)


### Pop increase (PC) ###

# function to convert shpfile to raster
shp2raster <- function(shp, mask.raster, label, value, transform = FALSE, proj.from = NA,
                       proj.to = NA, map = TRUE) {
  require(raster, rgdal)
  # use transform==TRUE if the polygon is not in the same coordinate system as
  # the output raster, setting proj.from & proj.to to the appropriate
  # projections
  if (transform == TRUE) {
    proj4string(shp) <- proj.from
    shp <- spTransform(shp, proj.to)
  }
  # convert the shapefile to a raster based on a standardised background
  # raster
  r <- rasterize(shp, mask.raster)
  # set the cells associated with the shapfile to the specified value
  r[!is.na(r)] <- value
  # merge the new raster with the mask raster and export to the working
  # directory as a tif file
  r <- mask(merge(r, mask.raster), mask.raster, filename = label, format = "GTiff",
            overwrite = T)
  # plot map of new raster
  if (map == TRUE) {
    plot(r, main = label, axes = F, box = F)
  }
  names(r) <- label
  return(r)
}

#open shpfile
setwd("~/Desktop/CH4_ms/shps")

#here I kept only areas that will have less than 50% growth from 2011-2050
no.growth.2050 <- readOGR("pop.change.shp")

#here I kept only areas with increasing population sizes
growth.50<- readOGR("high.pop.change.shp")

#pop.growth.2050 <- readOGR("GreenBook_Population_Growth_2050.shp")

#create raster mask
#HPc <- crop(HP, SA.ext)
HPm[!is.na(HPm)] <- 0

pop.change.raster <- shp2raster(shp = growth.50, mask.raster = HPm, value = 1, label= "Pop Increase", transform = TRUE, proj.from = "+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs", proj.to = "+proj=longlat +datum=WGS84 +no_defs")


#no.pop.change.raster <- shp2raster(shp = pop.growth.2050[grepl(c("Decrease|No Change|Medium"), pop.growth.2050$Pop_Pressu),], mask.raster = HPm, value = 1, label= "No Pop Increase", transform = TRUE, proj.from = "+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs", proj.to = "+proj=longlat +datum=WGS84 +no_defs")
#pop.change.raster <- shp2raster(shp = pop.growth.2050[grepl(c("High|Extreme|Medium"), pop.growth.2050$Pop_Pressu),], mask.raster = HPm, value = 1, label= "Pop Increase", transform = TRUE, proj.from = "+proj=utm +zone=35 +south +datum=WGS84 +units=m +no_defs", proj.to = "+proj=longlat +datum=WGS84 +no_defs")

PC <- aggregate(pop.change.raster, fact=3)

threshold=1
PC[PC > threshold] <- 0

### set it so areas with greater than zero are all 1's
PC[PC > 0] <- 1

#write raster
writeRaster(PC, "all.new.PC.tif", format="GTiff", overwrite=TRUE)



### Genetic Diversity (GD) ###
setwd("~/Desktop/CH4_ms/input_data")
library("readxl")
library(gstat)

cp_snp_models <- read_excel("./cp.snp.models.xlsx", 
                            col_types = c("text", "numeric", "numeric", 
                                          "numeric", "text", "text", "numeric", 
                                          "numeric", "numeric", "text"))


#Interpolate He along coast
dsp <- SpatialPoints(cp_snp_models[,2:3], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
dsp <- SpatialPointsDataFrame(dsp, cp_snp_models)
TA <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +no_defs")
dta <- spTransform(dsp, TA)
gs <- gstat(formula=He~1, locations=dsp)
p <- raster(MyBinCP_curr)
idw <- interpolate(p, gs)
CP.he.idw.m <- mask(idw, MyBinCP_curr)
plot(CP.he.idw.m)

CP.he <- CP.he.idw.m

#Get threshold
GD.cp.thres <-
  raster::quantile(CP.he,
                   probs = 0.5,
                   na.rm = TRUE)

CP.he[CP.he < GD.cp.thres] <- 0

#Do for PA
dsp <- SpatialPoints(pa_snp_models[,2:3], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
dsp <- SpatialPointsDataFrame(dsp, pa_snp_models)
TA <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +no_defs")
dta <- spTransform(dsp, TA)
gs <- gstat(formula=He~1, locations=dsp)
p <- raster(MyBinCP_curr)
idw <- interpolate(p, gs)
PA.he.idw.m <- mask(idw, MyBinCP_curr)
plot(PA.he.idw.m)

PA.he <- PA.he.idw.m

GD.pa.thres <-
  raster::quantile(PA.he,
                   probs = 0.5,
                   na.rm = TRUE)

PA.he[PA.he < GD.pa.thres] <- 0

#Do for SG
dsp <- SpatialPoints(sg_snp_models[,2:3], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
dsp <- SpatialPointsDataFrame(dsp, sg_snp_models)
TA <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +no_defs")
dta <- spTransform(dsp, TA)
gs <- gstat(formula=He~1, locations=dsp)
p <- raster(MyBinCP_curr)
idw <- interpolate(p, gs)
SG.he.idw.m <- mask(idw, MyBinCP_curr)
plot(SG.he.idw.m)

SG.he <- SG.he.idw.m

GD.sg.thres <-
  raster::quantile(SG.he,
                   probs = 0.5,
                   na.rm = TRUE)

SG.he[SG.he < GD.sg.thres] <- 0

writeRaster(CP.he, "CP.he.tif", format="GTiff", overwrite=TRUE)
writeRaster(PA.he, "PA.he.tif", format="GTiff", overwrite=TRUE)
writeRaster(SG.he, "SG.he.tif", format="GTiff", overwrite=TRUE)

all.He <- stack(CP.he, PA.he, SG.he)
writeRaster(all.He, "all.He.tif", format="GTiff", overwrite=TRUE)

### Genomic Vulnerability (GV) ###
setwd("~/Desktop/CH4_ms/shps")
GV.cp <- raster::stack("CP.gen.off.o50.tif")
GV.pa <- raster::stack("PA.gen.off.o50.tif")
GV.sg <- raster::stack("SG.gen.off.o50.tif")

GV.cp <- crop(GV.cp, SA.ext)
GV.pa <- crop(GV.pa, SA.ext)
GV.sg <- crop(GV.sg, SA.ext)

GV.cp.thres <-
  raster::quantile(GV.cp,
                   probs = 0.25,
                   na.rm = TRUE) # 0.0193385 


GV.pa.thres <-
  raster::quantile(GV.pa,
                   probs = 0.25,
                   na.rm = TRUE) # 0.006811526

GV.sg.thres <-
  raster::quantile(GV.sg,
                   probs = 0.25,
                   na.rm = TRUE) # 0.01394523

### clamp data to threshold such that values ABOVE threshold become
### zero and other values remain the same
GV.cp[GV.cp > GV.cp.thres] <- 0
GV.pa[GV.pa > GV.cp.thres] <- 0
GV.sg[GV.sg > GV.cp.thres] <- 0

### set it so areas with greater than zero are all 1's
GV.cp[GV.cp > 0] <- 1
GV.pa[GV.pa > 0] <- 1
GV.sg[GV.sg > 0] <- 1

GV.cp.ext <- extend(GV.cp, SA.ext, value=NA, snap="near")
GV.pa.ext <- extend(GV.pa, SA.ext, value=NA, snap="near")
GV.sg.ext <- extend(GV.sg, SA.ext, value=NA, snap="near")

GV.all <- stack(GV.cp.ext, GV.pa.ext, GV.sg.ext)
writeRaster(GV.cp.ext, "GV.cp.ext.tif", format="GTiff", overwrite=TRUE)
writeRaster(GV.pa.ext, "GV.pa.ext.tif", format="GTiff", overwrite=TRUE)
writeRaster(GV.sg.ext, "GV.sg.ext.tif", format="GTiff", overwrite=TRUE)

writeRaster(GV.all, "all.GV.tif", format="GTiff", overwrite=TRUE)

### Eco Regions ###
setwd("~/Desktop/CH4_ms/shps")
Agul <- raster::stack("Agulhas.tif")
Del <- raster::stack("Delagoa.tif")
Nam <- raster::stack("Namaqua.tif")
Nat <- raster::stack("Natal.tif")
SW <- raster::stack("Southwest.tif")
SB <- raster::stack("SouthernBeng.tif")
ND <- raster::stack("NatalDel.tif")

Bioreg <- stack(Agul, Del, Nam, Nat, SW, SB, ND)

writeRaster(Bioreg, "all.bioreg.tif", format="GTiff", overwrite=TRUE)

### Clim Stab (CS) ###
CS <- raster::stack("All.mar.diff.tif")

threshold <-
  raster::quantile(CS,
                   probs = 0.75,
                   na.rm = TRUE)
### clamp data to threshold such that values ABOVE threshold become
### zero and other values remain the same
CS2[CS > threshold] <- 0

writeRaster(CS, "all.new.clim.stab.tif", format="GTiff", overwrite=TRUE)

### Clim Var (CV) ###
CV <- raster::stack("All.mar.var.tif")

threshold <-
  raster::quantile(CV,
                   probs = 0.75,
                   na.rm = TRUE)

### clamp data to threshold such that values ABOVE threshold become
### zero and other values remain the same
CV[CV < threshold] <- 0

writeRaster(CV, "all.clim.var.tif", format="GTiff", overwrite=TRUE)


##############################################
#####    Part 3: Running prioritzr     #####
##############################################

library(tidyverse)
library(prioritizr)
library(rgdal)
library(raster)
library(rgeos)
library(mapview)
library(units)
library(scales)
library(assertthat)
library(gridExtra)

## setting up PU layer (this was made in QGIS with Qmarxan plug-in)
setwd("~/Desktop/CH4_ms/shps")
pu_data <- readOGR("pu.0.1.edt.shp")

setwd("~/Desktop/CH4_ms/Cumulative impacts/Cumulative impacts")
ras_CI <-raster("NBA2018_marine_cumulative_pressure.tif")
v1 <- raster::extract( ras_CI, pu_data, fun=mean, na.rm=TRUE)

colnames(pu_data@data)[7] = "locked_in"

pu_ID <- as.data.frame(pu_data@data[1])

CI.df <- cbind(pu_ID, v1)
pu_data_2 <-merge(pu_data,CI.df,by="puid")

names(pu_data_2)[8] <- "cost"

#give mpa areas (classified as 1's) label of locked in
mpa.df <- as.data.frame(pu_data_2@data[7])
mpa.df$mpa [mpa.df$locked_in == 1] <- "TRUE"
mpa.df$mpa [mpa.df$locked_in == 0] <- "FALSE"
names(mpa.df)[2] <- "locked_in"
MPA.df <- cbind(CI.df$puid, mpa.df)
names(MPA.df)[1] <- "puid"
df <- MPA.df[ -c(2) ]

pu_data_3 <-merge(pu_data_2,df,by="puid")

#give non-mpa areas (classified as 0's) label of locked out
mpa.df.2 <- as.data.frame(pu_data_3@data[7])
mpa.df.2$mpa [mpa.df.2$locked_in.x == 1] <- "FALSE"
mpa.df.2$mpa [mpa.df.2$locked_in.x == 0] <- "TRUE"
names(mpa.df.2)[2] <- "locked_out"
MPA.df <- cbind(CI.df$puid, mpa.df.2)
names(MPA.df)[1] <- "puid"
df <- MPA.df[ -c(2) ]

pu_data_4 <-merge(pu_data_3,df,by="puid")
pu_data_4 <- pu_data_4[, -(7)] 
names(pu_data_4)[9] <- "locked_in"

writeOGR(pu_data_4, ".", "pu_layer_edt", 
         driver = "ESRI Shapefile", overwrite_layer = TRUE) 




########## Setting up problems
# all scenarios have:
#  -minimum set objective
#  -target = 20% each feature (relative targets)
#  -binary decision type
# add add_boundary_penalties()

### create feature layer
setwd("~/Desktop/CH4_ms/input_data")
BR_all <-stack("all.bioreg.tif")
CS_old <-stack("all.clim.stab.tif")
GV_all <-stack("all.GV.tif")
GD_all <-stack("all.He.tif")
HFP_all <-stack("all.new.HFP.tif")
PC_all <-stack("all.new.PC.tif")
CV_all <-stack("all.clim.var.tif")
CS_all <-stack("all.new.clim.stab.tif")

GV_cp <-stack("GV.cp.ext.tif")
GV_pa <-stack("GV.pa.ext.tif")
GV_sg <-stack("GV.sg.ext.tif")

D <- stack("dummy.layer.tif")

#Re-sample the layers and stack into one raster stack [extent=SA.ext <- extent(14, 34, -35.5, -25)], [res=0.08333333]

t<-extent(16.39254, 32.90954, -34.90652, -26.85552 ) #layer extent from Bioreg
m<-extent(14, 34, -35.5, -25) #layer extent from CS_all
#no need to edit the following 6 lines
extent_list<-list(t, m)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)

reso<-res(CS_all) #choose layer you want to keep resolution
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=CS_all@crs) #choose layer crs you want to keep

BR_resamp <-resample(BR_all, s2, method="ngb") #resample by s2


#for GD_all
GD_allc <- crop(GD_all, GV_all)

#For HFP
t<-extent(13.99652, 33.99963, -35.50268, -24.9966) #layer extent from HFP
m<-extent(14, 34, -35.5, -25) #layer extent from CS_all
#no need to edit the following 6 lines
extent_list<-list(t, m)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)

reso<-res(CS_all) #choose layer you want to keep resolution
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=CS_all@crs) #choose layer crs you want to keep

HFP_resamp <-resample(HFP_all, s2, method="ngb") #resample by s2

#For PC
t<-extent(13.99652, 33.99963, -35.50268, -24.9966) #layer extent from HFP
m<-extent(14, 34, -35.5, -25) #layer extent from CS_all
#no need to edit the following 6 lines
extent_list<-list(t, m)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)

reso<-res(CS_all) #choose layer you want to keep resolution
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=CS_all@crs) #choose layer crs you want to keep

PC_resamp <-resample(PC_all, s2, method="ngb") #resample by s2


#For D
t<-extent(16.48, 32.997, -34.838, -26.87) #layer extent from HFP
m<-extent(14, 34, -35.5, -25) #layer extent from CS_all
#no need to edit the following 6 lines
extent_list<-list(t, m)
extent_list<-lapply(extent_list, as.matrix)
matrix_extent<-matrix(unlist(extent_list), ncol=length(extent_list))
rownames(matrix_extent)<-c("xmin", "ymin", "xmax", "ymax")
best_extent<-extent(min(matrix_extent[1,]), max(matrix_extent[3,]), min(matrix_extent[2,]), max(matrix_extent[4,]))
ranges<-apply(as.matrix(best_extent), 1, diff)

reso<-res(CS_all) #choose layer you want to keep resolution
nrow_ncol<-ranges/reso
s2<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=CS_all@crs) #choose layer crs you want to keep

D_resamp <-resample(D, s2, method="ngb") #resample by s2


ALL_p_f <- stack(BR_resamp, CV_all, GD_allc, HFP_resamp, CS_all, GV_cp, GV_pa, GV_sg, PC_resamp)

ALL_all <- stack(BR_resamp, CV_all, CS_all, GD_allc, GV_cp, GV_pa, GV_sg, HFP_resamp, PC_resamp)

env.all <- stack(CV_all, CS_all)
gen.all <- stack(GD_allc, GV_cp, GV_pa, GV_sg)
soc.all <- stack(HFP_resamp, PC_resamp)

pres.all <- stack(CV_all, GD_allc, HFP_resamp)
fut.all <- stack(CS_all, GV_cp, GV_pa, GV_sg, PC_resamp)

writeRaster(ALL_all, "ALL_all.tif", format="GTiff", overwrite=TRUE)
writeRaster(env.all, "env.all.tif", format="GTiff", overwrite=TRUE)
writeRaster(gen.all, "gen.all.tif", format="GTiff", overwrite=TRUE)
writeRaster(soc.all, "soc.all.tif", format="GTiff", overwrite=TRUE)

writeRaster(pres.all, "pres.all.tif", format="GTiff", overwrite=TRUE)
writeRaster(fut.all, "fut.all.tif", format="GTiff", overwrite=TRUE)


################# WORKING WITH PRIORITZR
setwd("~/Desktop/CH4_ms/input_data/new_inputs")

### get input feature rasters
env.all <-stack("env.all.tif")
gen.all <-stack("gen.all.tif")
soc.all <-stack("soc.all.tif")

pres.all <-stack("pres.all.tif")
fut.all <-stack("fut.all.tif")

# get pu_data layer
setwd("~/Desktop/CH4_ms/shps")
pu_data <- readOGR("pu_layer_edt.shp")
names(pu_data)[2] <- "pu_status"
names(pu_data)[8] <- "locked_in"
names(pu_data)[9] <- "locked_out"
pu_data$locked_in <- as.logical(pu_data$locked_in)
pu_data$locked_out <- as.logical(pu_data$locked_out)
print(pu_data)
plot(pu_data)
spplot(pu_data, "cost")


######### p0 = quasi gap analysis (doesn't lock out areas)

p0 <- problem(pu_data, ALL_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_default_solver(gap = 0)%>% 
  add_locked_in_constraints("locked_in")

s0 <- solve(p0)

features <- eval_feature_representation_summary(p0, s0[, "solution_1"])
targets <- eval_target_coverage_summary(p0, s0[, "solution_1"])
writexl::write_xlsx(targets,"20perc_targets.gap.xlsx")

pdf("existing.mpas.sol.pdf")
spplot(s0, "solution_1", col.regions = c("white", "darkgreen"), main = "s0")
dev.off()


### Actual Gap Analysis
p0 <- problem(pu_data, ALL_all, cost_column = "cost")
print(p0)

abundance_data <- feature_abundances(p0)
print(abundance_data)

# add new column with feature abundances in km^2
abundance_data$absolute_abundance_km2 <-
  (abundance_data$absolute_abundance * prod(res(ALL_all))) %>%
  set_units(m^2) %>%
  set_units(km^2)

print(abundance_data)

# plot histogram of the features' abundances
hist(abundance_data$absolute_abundance_km2, main = "Feature abundances")

# create column in planning unit data with binary values (zeros and ones)
# indicating if a planning unit is covered by protected areas or not
pu_data$pa_status <- as.numeric(pu_data$locked_in)

# calculate feature representation by protected areas
repr_data <- eval_feature_representation_summary(p0, pu_data[, "pa_status"])
print(repr_data)

write.table(repr_data, "new.gap.anal.txt")

######### p1 = Present_all

p1 <- problem(pu_data, pres.all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s1 <- solve(p1)

spplot(s1, "solution_1", col.regions = c("white", "darkgreen"), main = "P")

solution_columns <- which(grepl("solution", names(s1)))
s1$selection_frequencies <- rowSums(as.matrix(s1@data[, solution_columns]))
P_all.sf <- as.data.frame(s1$selection_frequencies)


######### p2 = Future_all

p2 <- problem(pu_data, fut.all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s2 <- solve(p2)

spplot(s2, "solution_1", col.regions = c("white", "darkgreen"), main = "F")

solution_columns <- which(grepl("solution", names(s2)))
s2$selection_frequencies <- rowSums(as.matrix(s2@data[, solution_columns]))
F_all.sf <- as.data.frame(s2$selection_frequencies)

######### p3 = CV

p3 <- problem(pu_data, CV_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s3 <- solve(p3)
spplot(s3, "solution_1", col.regions = c("white", "darkgreen"), main = "BR")

solution_columns <- which(grepl("solution", names(s3)))
s3$selection_frequencies <- rowSums(as.matrix(s3@data[, solution_columns]))
CV_all.sf <- as.data.frame(s3$selection_frequencies)

######### p4 = CS

p4 <- problem(pu_data, CS_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s4 <- solve(p4)
spplot(s4, "solution_1", col.regions = c("white", "darkgreen"), main = "CS")

solution_columns <- which(grepl("solution", names(s4)))
s4$selection_frequencies <- rowSums(as.matrix(s4@data[, solution_columns]))
CS_all.sf <- as.data.frame(s4$selection_frequencies)

######### p5 = GD

p5 <- problem(pu_data, GD_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s5 <- solve(p5)
spplot(s5, "solution_1", col.regions = c("white", "darkgreen"), main = "GD")

solution_columns <- which(grepl("solution", names(s5)))
s5$selection_frequencies <- rowSums(as.matrix(s5@data[, solution_columns]))
GD_all.sf <- as.data.frame(s5$selection_frequencies)

######### p6 = GV

p6 <- problem(pu_data, GV_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s6 <- solve(p6)
spplot(s6, "solution_1", col.regions = c("white", "darkgreen"), main = "GV")

solution_columns <- which(grepl("solution", names(s6)))
s6$selection_frequencies <- rowSums(as.matrix(s6@data[, solution_columns]))
GV_all.sf <- as.data.frame(s6$selection_frequencies)

######### p7 = HFP

p7 <- problem(pu_data, HFP_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s7 <- solve(p7)
spplot(s7, "solution_1", col.regions = c("white", "darkgreen"), main = "HFP")

solution_columns <- which(grepl("solution", names(s7)))
s7$selection_frequencies <- rowSums(as.matrix(s7@data[, solution_columns]))
HFP_all.sf <- as.data.frame(s7$selection_frequencies)

######### p8 = PC

p8 <- problem(pu_data, PC_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s8 <- solve(p8)
spplot(s8, "solution_1", col.regions = c("white", "darkgreen"), main = "PC")

solution_columns <- which(grepl("solution", names(s8)))
s8$selection_frequencies <- rowSums(as.matrix(s8@data[, solution_columns]))
PC_all.sf <- as.data.frame(s8$selection_frequencies)


#export selection freqs for all scenarios
SF.all.feat_0.2 <- cbind(P_all.sf, F_all.sf, BR_all.sf, CS_all.sf, GD_all.sf, GV_all.sf, HFP_all.sf, PC_all.sf)

pu_id <- as.data.frame(pu_data@data$puid)
SF.all.feat_0.2 <- cbind(pu_id, SF.all.feat_0.2)

writexl::write_xlsx(SF.all.feat_0.2,"SF.all_0.2.xlsx")


### Run pres and fut all without portfolio 100 solutions

######### p9 = Present_all

p9 <- problem(pu_data, pres.all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

s9 <- solve(p9)

P_pu_sel <- as.data.frame(s9@data$solution_1)

######### p10 = Future_all

p10 <- problem(pu_data, fut.all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

s10 <- solve(p10)

F_pu_sel <- as.data.frame(s10@data$solution_1)

#export dataframe
P_F_sols <- cbind(pu_id, P_pu_sel, F_pu_sel)

writexl::write_xlsx(P_F_sols,"new_P_F_sols.xlsx")

######### p11 = ALL_all

p11 <- problem(pu_data, ALL_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s11 <- solve(p11)

############ Get 100 solution feature representation
#this will create a df with % held per feature, repeated over the 100 solutions
#this is done with p11 which is ALL_all
y  <- NULL; 
for(i in 1:100){tmp <- eval_feature_representation_summary(p11, s11[, paste0("solution_",i)])
y <- rbind(y, tmp)}

#create column z, which gives solution number (change the "each" number to # of feat)
z<-rep(1:100,each=16)
df<-data.frame(z)

sol.feats <- cbind(y,df)


#filter the df to only include bioregions
bioreg.feat <- sol.feats[grep("bioreg", sol.feats$feature), ]

bioreg.feat <- bioreg.feat %>%
  group_by(z) %>%
  mutate(Bioreg_held = mean(relative_held))

#do for He
He.feat <- sol.feats[grep("He", sol.feats$feature), ]

He.feat <- He.feat %>%
  group_by(z) %>%
  mutate(GD_held = mean(relative_held))

#do for GV
GV.feat <- sol.feats[grep("GV", sol.feats$feature), ]

GV.feat <- GV.feat %>%
  group_by(z) %>%
  mutate(GV_held = mean(relative_held))


GV.ft<-GV.feat[ -c(1:5) ]
He.ft<-He.feat[ -c(1:5) ]
BR.ft<-bioreg.feat[ -c(1:5) ]

df_merge <- merge(GV.ft, He.ft, by="z")
he_gv_br_merge <- merge(df_merge, BR.ft, by="z")

feat.held <-unique(he_gv_br_merge)

#filter for HFP
HFP.feat <- sol.feats[grep("HFP", sol.feats$feature), ]

names(HFP.feat)[names(HFP.feat) == 'relative_held'] <- 'HFP_held'

#filter for PC
PC.feat <- sol.feats[grep("PC", sol.feats$feature), ]

names(PC.feat)[names(PC.feat) == 'relative_held'] <- 'PC_held'

#filter for CS
CS.feat <- sol.feats[grep("clim.stab", sol.feats$feature), ]

names(CS.feat)[names(CS.feat) == 'relative_held'] <- 'CS_held'

HFP.ft<-HFP.feat[ -c(1:4) ]
PC.ft<-PC.feat[ -c(1:4) ]
CS.ft<-CS.feat[ -c(1:4) ]

pc_hfp_merge <- merge(PC.ft, HFP.ft, by="z")

pc_hfp_cs_merge <- merge(CS.ft, pc_hfp_merge, by="z")

all_feat_merge <- merge(pc_hfp_cs_merge, feat.held, by="z")

writexl::write_xlsx(all_feat_merge,"P11_feats_held.xlsx")

writexl::write_xlsx(all_feat_merge,"all_feats_rel_held.xlsx")

names(all_feat_merge)[2] <- 'CS'
names(all_feat_merge)[3] <- 'PC'
names(all_feat_merge)[4] <- 'HFP'
names(all_feat_merge)[5] <- 'GV'
names(all_feat_merge)[6] <- 'GD'
names(all_feat_merge)[7] <- 'BR'

######### Plotting parallel coordinates

#packages = c('GGally', 'plotly', 'parcoords', 'tidyverse')

#for(p in packages){
#  if(!require(p, character.only = T)){
#    install.packages(p)
#  }
#}

### Using P11_feat input here!! (latest one)

#create parallele coords plot with scale= "std": univariately, subtract mean and divide by standard deviation
pdf("parallel_coords_ALL.pdf")

p <- ggparcoord(data = P11_feats_held,
           columns = 2:7,
           alphaLines = 0.2,
           showPoints = TRUE,
           boxplot = TRUE,
           shadeBox = 3)

p + theme_minimal()

dev.off()

#create parallele coords plot with scale= "uniminmax": univariately, scale so the minimum of the variable is zero, and the maximum is one
p <- ggparcoord(data = all_feat_merge,
                columns = 2:7,
                alphaLines = 0.2,
                showPoints = TRUE,
                boxplot = TRUE,
                scale = "uniminmax")

p + theme_minimal()


######## Getting irreplacability

## No portfolio - CV
cv <- problem(pu_data, CV_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

s_cv <- solve(cv)

rc.cv <- cv %>% eval_replacement_importance(s_cv[, "solution_1"])

CV_rc <- as.data.frame(rc.cv$rc)

## No portfolio - HP
hp <- problem(pu_data, HFP_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

s_hp <- solve(hp)

rc.hp <- hp %>% eval_replacement_importance(s_hp[, "solution_1"])

HP_rc <- as.data.frame(rc.hp$rc)

## No portfolio - GD
gd <- problem(pu_data, GD_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

s_gd <- solve(gd)

rc.gd <- gd %>% eval_replacement_importance(s_gd[, "solution_1"])

GD_rc <- as.data.frame(rc.gd$rc)

#merge for past, save as excel
BR_HP <-cbind(BR_rc, HP_rc)

br_hp_gd <- cbind(BR_HP, GD_rc)

writexl::write_xlsx(br_hp_gd,"br_hp_gd.xlsx")


#### For future

## No portfolio - CS
cs <- problem(pu_data, CS_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

s_cs <- solve(cs)

rc.cs <- cs %>% eval_replacement_importance(s_cs[, "solution_1"])

CS_rc <- as.data.frame(rc.cs$rc)

## No portfolio - GV
gv <- problem(pu_data, GV_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

s_gv <- solve(gv)

rc.gv <- gv %>% eval_replacement_importance(s_gv[, "solution_1"])

GV_rc <- as.data.frame(rc.gv$rc)


## No portfolio - PC
pc <- problem(pu_data, PC_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

s_pc <- solve(pc)

rc.pc <- pc %>% eval_replacement_importance(s_pc[, "solution_1"])

PC_rc <- as.data.frame(rc.pc$rc)



#merge for future- export of excel
CS_GV <-cbind(CS_rc, GV_rc)

cs_gv_pc <- cbind(CS_GV, PC_rc)

writexl::write_xlsx(cs_gv_pc,"cs_gv_pc.xlsx")

#### HERE is where I took the irreplacability values (with selection freqs) of the 
#### 6 independ feat scens, and comapred to the binary solutions of P_all or F_all to 
#### get df with pres/abs of solution, and what feat each pres PU is most important for
setwd("~/Desktop/CH4_ms/prioritzr_sols/new_sols")

P_sols <- read_excel("new_P_F_sols.xlsx", sheet = "Sheet4")

#give 4 classes numeric values
P_sols$most_imp_P[P_sols$most_imp_P == "NA"] <- "0"
P_sols$most_imp_P[P_sols$most_imp_P == "CV"] <- "1"
P_sols$most_imp_P[P_sols$most_imp_P == "GD"] <- "2"
P_sols$most_imp_P[P_sols$most_imp_P == "HF"] <- "3"

#change class to num
P_sols$most_imp_P <- as.numeric(P_sols$most_imp_P)

#merge these values with pu_data df
pu_p_imp <- merge(pu_data, P_sols, by = "puid")


#Do same for future feats
F_sols <- read_excel("new_P_F_sols.xlsx", sheet = "Sheet5")

F_sols$most_imp_F[F_sols$most_imp_F == "NA"] <- "0"
F_sols$most_imp_F[F_sols$most_imp_F == "CS"] <- "1"
F_sols$most_imp_F[F_sols$most_imp_F == "GV"] <- "2"
F_sols$most_imp_F[F_sols$most_imp_F == "PC"] <- "3"

F_sols$most_imp_F <- as.numeric(F_sols$most_imp_F)

pu_p_imps <- merge(pu_p_imp, F_sols, by = "puid")

writeOGR(pu_p_imps, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols", layer = "pu_p.f_imp_vals.shp",
         driver = "ESRI Shapefile" )


#read shp
setwd("~/Desktop/CH4_ms/prioritzr_sols/new_sols")
imp_shp <-readOGR(".","pu_p.f_imp_vals.shp")

#rasterise for most important Present features
ext  = extent(SA.ext <- extent(14, 34, -35.5, -25))
gridsize  = 0.008333333
r         = raster(ext, res = gridsize)
Pvals.r <- rasterize(imp_shp, field="most_imp_P", r)

writeRaster(Pvals.r, "Pvals.r.tif", format="GTiff", overwrite=TRUE)


#rasterise for most important Future features

r         = raster(ext, res = gridsize)
Fvals.r <- rasterize(imp_shp, field="most_imp_F", r)

writeRaster(Fvals.r, "Fvals.r.tif", format="GTiff", overwrite=TRUE)

#save ALL_all selection freq as shp (to then rasterize)

solution_columns <- which(grepl("solution", names(s11)))
s11$selection_frequencies <- rowSums(as.matrix(s11@data[, solution_columns]))
ALL_all.sf <- as.data.frame(s11$selection_frequencies)

pu_all_sf <- cbind(pu_data, ALL_all.sf)
#pu_all_sf <- merge(pu_data, ALL_all.sf, by = "puid")

writeOGR(pu_all_sf, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols", layer = "pu_all_sf.shp",
         driver = "ESRI Shapefile" )

all_sf <-readOGR(".","pu_all_sf.shp")

r         = raster(ext, res = gridsize)
all_sf.r <- rasterize(all_sf, field="s11_sl_", r)

writeRaster(all_sf.r, "ALL_all_sf.tif", format="GTiff", overwrite=TRUE)


############## Pareto Fronts
#climate
show_front <- function(pref) {
  plot(all_feats_rel_held$Bioreg_held, all_feats_rel_held$CS_held)
  sky <- psel(all_feats_rel_held, pref)
  plot_front(all_feats_rel_held, pref, col = rgb(0, 0, 1))
  points(sky$Bioreg_held, sky$CS_held, lwd = 3)
}

show_front(high(Bioreg_held) * high(CS_held))

pdf("climate.pareto.pdf")
show_front(high(Bioreg_held) | high(CS_held))
dev.off()


#genetic
show_front <- function(pref) {
  plot(all_feats_rel_held$GD_held, all_feats_rel_held$GV_held)
  sky <- psel(all_feats_rel_held, pref)
  plot_front(all_feats_rel_held, pref, col = rgb(0, 0, 1))
  points(sky$GD_held, sky$GV_held, lwd = 3)
}

show_front(high(GD_held) * high(GV_held))

pdf("genetic.pareto.pdf")
show_front(high(GD_held) | high(GV_held))
dev.off()

#social
show_front <- function(pref) {
  plot(all_feats_rel_held$HFP_held, all_feats_rel_held$PC_held)
  sky <- psel(all_feats_rel_held, pref)
  plot_front(all_feats_rel_held, pref, col = rgb(0, 0, 1))
  points(sky$HFP_held, sky$PC_held, lwd = 3)
}

show_front(high(GD_held) * high(GV_held))

pdf("social.pareto.pdf")
show_front(high(HFP_held) | high(PC_held))
dev.off()



# to create pareto front plot with ggplot:

p <- high(Bioreg_held) * high(CS_held) 

res <- psel(all_feats_rel_held, p, top = nrow(all_feats_rel_held))

gp <- ggplot(res, aes(x = Bioreg_held, y = CS_held, color = factor(.level))) + geom_point(size = 3)



# Running future scenario with no targets for future

gen <- problem(pu_data, ALL_p_f, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_manual_targets(data.frame(
    feature = names(ALL_p_f)[1:11], type = "relative",
    sense = ">=", target = 0.2))%>%
  add_binary_decisions()%>% 
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

s12 <- solve(p12)

features <- eval_feature_representation_summary(p12, s12[, "solution_1"])
targets <- eval_target_coverage_summary(p12, s12[, "solution_1"])
writexl::write_xlsx(features,"pres4fut.held.gap.xlsx")


#### Get pareto front for gen/clim/soc objectives

ALL_clim <- stack(BR_resamp, CS_all, GD_allc, HFP_resamp,  GV_cp, GV_pa, GV_sg, PC_resamp)
ALL_gen <- stack(GD_allc, GV_cp, GV_pa, GV_sg, BR_resamp, CS_all, HFP_resamp,  PC_resamp)
ALL_soc <- stack(HFP_resamp, PC_resamp, BR_resamp, CS_all, GD_allc, GV_cp, GV_pa, GV_sg)


######### gen

pf <- problem(pu_data, ALL_p_f, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_manual_targets(data.frame(
    feature = names(ALL_soc)[12:16], type = "relative",
    sense = ">=", target = 0.2))%>%
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s.pf <- solve(pf)

############ Get 100 solution feature representation
#this will create a df with % held per feature, repeated over the 100 solutions
#this is done with p11 which is ALL_all
y  <- NULL; 
for(i in 1:100){tmp <- eval_feature_representation_summary(pf, s.pf[, paste0("solution_",i)])
y <- rbind(y, tmp)}

#create column z, which gives solution number (change the "each" number to # of feat)
z<-rep(1:100,each=16)
df<-data.frame(z)

sol.feats <- cbind(y,df)


#filter the df to only include bioregions
bioreg.feat <- sol.feats[grep("bioreg", sol.feats$feature), ]

bioreg.feat <- bioreg.feat %>%
  group_by(z) %>%
  mutate(Bioreg_held = mean(relative_held))

#do for He
He.feat <- sol.feats[grep("He", sol.feats$feature), ]

He.feat <- He.feat %>%
  group_by(z) %>%
  mutate(GD_held = mean(relative_held))

#do for GV
GV.feat <- sol.feats[grep("GV", sol.feats$feature), ]

GV.feat <- GV.feat %>%
  group_by(z) %>%
  mutate(GV_held = mean(relative_held))


GV.ft<-GV.feat[ -c(1:5) ]
He.ft<-He.feat[ -c(1:5) ]
BR.ft<-bioreg.feat[ -c(1:5) ]

df_merge <- merge(GV.ft, He.ft, by="z")
he_gv_br_merge <- merge(df_merge, BR.ft, by="z")

feat.held <-unique(he_gv_br_merge)

#filter for HFP
HFP.feat <- sol.feats[grep("HFP", sol.feats$feature), ]

names(HFP.feat)[names(HFP.feat) == 'relative_held'] <- 'HFP_held'

#filter for PC
PC.feat <- sol.feats[grep("PC", sol.feats$feature), ]

names(PC.feat)[names(PC.feat) == 'relative_held'] <- 'PC_held'

#filter for CS
CS.feat <- sol.feats[grep("clim.stab", sol.feats$feature), ]

names(CS.feat)[names(CS.feat) == 'relative_held'] <- 'CS_held'

HFP.ft<-HFP.feat[ -c(1:4) ]
PC.ft<-PC.feat[ -c(1:4) ]
CS.ft<-CS.feat[ -c(1:4) ]

pc_hfp_merge <- merge(PC.ft, HFP.ft, by="z")

pc_hfp_cs_merge <- merge(CS.ft, pc_hfp_merge, by="z")

all_feat_merge <- merge(pc_hfp_cs_merge, feat.held, by="z")

writexl::write_xlsx(all_feat_merge,"just_fut_feats_held.xlsx")

######
###### plotting radar charts
######

scenario_radar_df <- read_excel("scenario.radar.df.xlsx")

sc_df <-mutate_if(scenario_radar_df, is.numeric, ~ . * 100)

#then edited in excel, giving max and min rows (100, and 0s respectively)

data2 <- read.xlsx2("scen.4.radar.xlsx",  # Column 1 as row names
                    +                     sheetIndex = 1)

#didnt work because values were characters
sapply(data2, class)

#change df to numeric
data2[ , i] <- apply(data2[ , i], 2,            # Specify own function within apply
                       function(x) as.numeric(as.character(x)))

clim_data <- data2[c(1, 2, 3), ]
radarchart(clim_data)

#function to plot nice radar plots
library(fmsb)
create_beautiful_radarchart <- function(data, color = "green2", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
  
}

#read data for radar plot w/ MPA overlap per scenario
overlap_radar <- read_excel("new_sols/overlap.radar.xlsx")
overlap_radar <- overlap_radar[,2:7]

pdf("MPA.radar.plot.new.axes.pdf")
op <- par(mar = c(1, 2, 2, 1))
create_beautiful_radarchart(overlap_radar)
par(op)
dev.off()


# plot objectives
# Define colors and titles
#colors <- c("#E7B800", "#FC4E07", "#66FF00", "#00688B","#EE1289")
#titles <- c("GD", "GV", "HF", "PC", "CV", "CS")

# Reduce plot margin using par()
# Split the screen in 3 parts
#op <- par(mar = c(1, 1, 1, 1))
#par(mfrow = c(2,3))

# Create the radar chart
#for(i in 1:6){
#  create_beautiful_radarchart(
#    data = data2[c(1, 2, i+2), ], caxislabels = c(0, 5, 10, 15, 20),
#    color = colors[i], title = titles[i]
 # )
#}

#par(op)

######
### Getting irreplacability for radar/pareto plots
######

gen <- problem(pu_data, gen.all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% add_default_solver(gap = 0, verbose = FALSE)

gen <- solve(gen)

gen_pu_sel <- as.data.frame(gen@data$solution_1)

#export dataframe
obj_sols <- cbind(pu_id, gen_pu_sel, clim_pu_sel, soc_pu_sel, all_pu_sel)

writexl::write_xlsx(obj_sols,"obj_sols.xlsx")



### editting SF dataframe to get % captured by present day solution
SF_all_0_2 <- read_excel("Desktop/CH4_ms/prioritzr_sols/SF.all_0.2.xlsx")

SF_all_0_2 <- SF_all_0_2 %>% mutate(CSsol =
                                      case_when(CS_sf >= 55 ~ "1", CS_sf < 55 ~ "0",)
)

nrow(sols_for_pi_charts[sols_for_pi_charts$PC_caught > 1, ])

View(sols_for_pi_charts2) 

#creating upset plot
library(UpSetR)

individ_sols <- read_excel("Desktop/CH4_ms/prioritzr_sols/individ_sols.xlsx")

sols <- as.data.frame(individ_sols)

upset(sols, sets = c("BR", "HF", "GD", "CS", "PC", "GV"), order.by = "degree")

####
######## Mapping ALL scen with MPAs ontop ########
####
library(raster)
library(tmap)
library(rgdal)
library(rgeos)
library(sf)

setwd("~/Desktop/CH4_ms/prioritzr_sols/new_sols")
ALL <- raster::stack("ALL_all_sf.tif")

setwd("~/Desktop/CH4_ms/shps")
mpas <- readOGR("coast.MPAs.shp")

AF <- readOGR("Africa.shp")
SA.ext <- extent(14, 34, -35.5, -25)
map4 <- crop(AF, SA.ext)


#go from SF to binary solutions
ALL[ALL > 50] <- 100
ALL[ALL < 50] <- 0


t.ras <- tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(ALL) + 
 tm_raster(title = "Selection", palette = c("thistle1", "darkblue"))+
  tm_legend(outside = TRUE) 

t.ras + tm_shape(mpas) + tm_polygons(col="white", border.col="green2") + tm_scale_bar(breaks = c(0, 100, 200), text.size = 1) +
  tm_compass(type = "4star", size = 1, position = c("left", "top"))


#Mapping individ scenarios with MPAs ontop
setwd("~/Desktop/CH4_ms/prioritzr_sols/new_sols")
CV <- readOGR("pu_CV_sol.shp.shp")

tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(CV) + tm_fill("CV", palette = "darkblue")+ tm_shape(mpas) + tm_polygons(col="white", border.col="green2") + tm_scale_bar(breaks = c(0, 100, 200), text.size = 1) +
  tm_compass(type = "4star", size = 1, position = c("left", "top"))


#create P/F by objective map
setwd("~/Desktop/CH4_ms/prioritzr_sols/new_sols")
Pval <- raster::stack("Pvals.r.tif")

p.ras <- tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(Pval) + 
  tm_raster(col = "Pvals.r", style = "fixed", breaks = c(0, 1, 2, 3, 4), palette = c("lightgrey", "#00688B", "#EE1289", "#FFC125"))+
  tm_legend(outside = TRUE) 

pdf("P.val.map.pdf")
p.ras
dev.off()

Fval <- raster::stack("Fvals.r.tif")

f.ras <- tm_shape(Fval) + 
  tm_raster(col = "Fvals.r", style = "fixed", breaks = c(0, 1, 2, 3, 4), palette = c("lightgrey", "#00688B", "#EE1289", "#FFC125"))+
  tm_legend(outside = TRUE) 

pdf("F.val.map.pdf")
f.ras
dev.off()

#####
######### create cum sol map ######### 
#####

plas.col <- plasma(6, direction = -1)

cum_sols <- read_excel("cum_sols.xlsx", sheet = "Sheet2", 
                        col_types = c("numeric", "numeric"))

pu_p_imp <- merge(pu_data, cum_sols, by = "puid")
writeOGR(pu_p_imp, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols", layer = "pu_cum_obj",
         driver = "ESRI Shapefile" )

setwd("~/Desktop/CH4_ms/prioritzr_sols/new_sols")
cum_sf <- read_sf("pu_cum_obj.shp")

cum.ras <- tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(cum_sf) +
  tm_polygons(col = "cum_sol", palette = plas.col, style = "fixed", breaks = c(0, 1, 2, 3, 4, 5, 6)) + tm_scale_bar(breaks = c(0, 100, 200), text.size = 0.5, position = c("left", "bottom")) +
  tm_compass(type = "4star", size = 1, position = c("left", "top"))

pdf("cum_sols_map.pdf")
cum.ras
dev.off()

######### 
### Running scenarios with % endangered
######### 
setwd("~/Desktop/CH4_ms/shps")
endg <- readOGR("SA.vuln.18.shp")

endg2 <- spTransform(endg, crs(pu_data))
r <- raster(endg2,nrows=126,ncols=240)
r <- rasterize(endg2,r,fun="first")

p11 <- problem(pu_data, ALL_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

p_all <- p11 %>% add_locked_out_constraints(r)
s_all <- solve(p_all)

solution_columns <- which(grepl("solution", names(s_all)))
s_all$selection_frequencies <- rowSums(as.matrix(s_all@data[, solution_columns]))
eng_all.sf <- as.data.frame(s_all$selection_frequencies)

#spplot(s_all, "solution_1", col.regions = c("white", "darkgreen"), main = "s0")


pu_all_e_sf <- cbind(pu_data, eng_all.sf)
#pu_all_sf <- merge(pu_data, ALL_all.sf, by = "puid")

writeOGR(pu_all_e_sf, dsn = "~/Desktop/CH4_ms/prioritzr_sols", layer = "pu_all_e_sf.shp",
         driver = "ESRI Shapefile" )

all_sf_e <-readOGR(".","pu_all_e_sf.shp")
ext  = extent(SA.ext <- extent(14, 34, -35.5, -25))
r         = raster(ext, res = gridsize)
all_sf_e.r <- rasterize(all_sf_e, field="s_ll_s_", r)

writeRaster(all_sf_e.r, "ALL_endg_sf.tif", format="GTiff", overwrite=TRUE)

 library(raster)
 library(tmap)
 library(rgdal)
 library(rgeos)
 library(sf)
 setwd("~/Desktop/CH4_ms/shps")

 AF <- readOGR("Africa.shp")
 SA.ext <- extent(14, 34, -35.5, -25)
 map4 <- crop(AF, SA.ext)
 t.ras <- tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(all_sf_e.r) + 
  tm_raster(style = "cont", title = "Elevation (m)", palette = c("thistle1", "darkblue"))+
  tm_legend(outside = TRUE) 

 pdf("ALL.endg.solmap.pdf")
 t.ras
 dev.off()
 
#plot with soc vuln
setwd("~/Desktop/CH4_ms/shps/GreenBook")
SV <- readOGR("GreenBook_SocioEconomic_Vulnerability_2011.shp")
 
t.shp <- tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(SV) + tm_polygons(col = "SEV11")
  
######### 
#### plot SV with bioregions
######### 
setwd("~/Desktop/CH4_ms/shps")
BZ <- readOGR("biozones.coast.shp")

t.shp<-tm_shape(map4, bbox = c(14, 34, -40, -25)) + tm_polygons(col="white", border.col="grey") + tm_shape(SV, bbox = c(14, 34, -40, -25)) + tm_polygons(col = "SEV11", palette=cols)

t.bz.sv <-  t.shp + tm_shape(BZ, bbox = c(14, 34, -40, -25)) + tm_polygons(col = "BIOREGION")

pdf("SV.BZ.map.pdf")
t.bz.sv 
dev.off()

#make inset map for fig 1 (zoom in on delagoa bioregion)
cols= viridis(5)
t.shp<-tm_shape(map4, bbox = c(30, 34, -30, -26)) + tm_polygons(col="white", border.col="grey") + tm_shape(SV, bbox = c(30, 34, -30, -26)) + tm_polygons(col = "SEV11", palette=cols, legend.show = FALSE)

bz.col<-brewer.pal(n = 5, name = "Pastel2")
t.bz.sv2 <-  t.shp + tm_shape(BZ, bbox = c(30, 34, -30, -26)) + tm_polygons(col = "BIOREGION", palette=bz.col, legend.show = FALSE)

pdf("SV.BZ.inset.map.pdf")
t.bz.sv2 
dev.off()

#t.ras <-  t.shp + tm_shape(all_sf_e.r) + tm_raster(style = "cont", title = "Elevation (m)", palette = c("white", "darkblue"))+tm_legend(outside = TRUE)

pdf("ALL.endg.SV.solmap.pdf")
t.ras
dev.off()
 
tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(endg) + tm_polygons(col = "#EE1289")

cols= viridis(5)
t.shp<-tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(SV) + tm_polygons(col = "SEV11", palette=cols)

######### 
##### Plot Gen importance maps #####
######### 
new_gen_pu_sol <- read_excel("new_gen_pu_sol.xlsx",  sheet = "Sheet2")
gen.sols <- merge(pu_p_imp, new_gen_pu_sol, by = "puid")
writeOGR(gen.sols, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols", layer = "gen.sols.2.shp",
         driver = "ESRI Shapefile" )

setwd("~/Desktop/CH4_ms/prioritzr_sols/new_sols")
gen_shp <-readOGR(".","gen.sols.2.shp")

GDcols <-c( "#F7F7F7", "#8DEEEE", "#008B8B", "#440154FF")
New.cols<-c( "white", "darkorchid2", "darkgoldenrod1", "cyan4")
GVcols <-c("")
tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(gen_shp) + tm_polygons("GV_spp", border.col="lightgrey", palette = GDcols, style = "fixed", breaks = c(0, 1, 2, 3, 4))

######### 
#### Venn Diagrams for spp GV GD comparison ####
######### 
library(ggVennDiagram)

new_gen_pu_sol <- read_excel("Desktop/CH4_ms/prioritzr_sols/new_sols/new_gen_pu_sol.xlsx", 
                       sheet = "Sheet3")

gen.df <- as.data.frame(new_gen_pu_sol)

g.df<-gen.df %>% mutate_if(is.numeric, ~1 * (. > 50))

ggVennDiagram(lapply(g.df[,1:3], function(x) which(x == 1)),label = "count", cex=2)+ scale_color_brewer(palette = "Blues")


## heirarchical clustering of selection freqs 
SF_new_0_2 <- read_excel("SF.new_0.2.xlsx", 
                         +     sheet = "Sheet3")
distance_mat <- dist(SF_new_0_2, method = 'euclidean')
Hierar_cl <- hclust(distance_mat, method = "average")
plot(Hierar_cl)

######### 
###### Gap analysis with mult PU sizes
######### 
setwd("~/Desktop/CH4_ms/shps")
prt <- readOGR("SA.prot_18.shp")

prt2 <- spTransform(prt, crs(pu_data))
pr <- raster(prt2,nrows=126,ncols=240)
pr <- rasterize(prt2,pr,fun="first")

p0.1 <- problem(pu_data, ALL_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_binary_decisions()%>% 
  add_default_solver(gap = 0)

p_0.1 <- p0.1 %>% add_locked_in_constraints(pr)
s_0.1 <- solve(p_0.1)

#features <- eval_feature_representation_summary(p0.1, s0.1[, "solution_1"])
targets <- eval_target_coverage_summary(p_0.1, s_0.1[, "solution_1"])
writexl::write_xlsx(targets,"pu.0.1.targets.gap.xlsx")


## do for other layer
setwd("~/Desktop/CH4_ms/shps")
pu_data <- readOGR("pu.test.0.05.shp")

setwd("~/Desktop/CH4_ms/Cumulative impacts/Cumulative impacts")
ras_CI <-raster("NBA2018_marine_cumulative_pressure.tif")
v1 <- raster::extract( ras_CI, pu_data, fun=mean, na.rm=TRUE)

pu_ID <- as.data.frame(pu_data@data[1])

CI.df <- cbind(pu_ID, v1)
pu_data_5 <-merge(pu_data,CI.df,by="puid")

names(pu_data_5)[7] <- "cost"
names(pu_data_5)[2] <- "pu_status"
spplot(pu_data_5, "cost")

p0.05 <- problem(pu_data_5, ALL_all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_binary_decisions()%>% 
  add_default_solver(gap = 0)

p_0.05 <- p0.05 %>% add_locked_in_constraints(pr)
s_0.05 <- solve(p_0.05)

#features <- eval_feature_representation_summary(p_0.05, s_0.05[, "solution_1"])
targets <- eval_target_coverage_summary(p_0.025, s_0.025[, "solution_1"])
writexl::write_xlsx(features,"pu.0.025.targets.gap.xlsx")



#colnames(pu_data@data)[7] = "locked_in"
#writeOGR(pu_data_4, ".", "pu_layer_edt", 
#         driver = "ESRI Shapefile", overwrite_layer = TRUE)


######### 
#### Kappa Stats ####
######### 
library(psych)

cum_sols <- read_excel("cum_sols.xlsx")
D_df<- as.data.frame(cum_sols[,c(2,3)])
cohen.kappa(D_df)

kappa_stats <- read_excel("kappa.stats.xlsx", 
                             sheet = "Sheet3", col_types = c("text", 
                                                                     "text", "numeric"))

p <- ggplot(data = kappa_stats, aes(x = scenarios, y = kappa)) + geom_point()
p + facet_wrap(~comparison, scale="free")+ scale_y_continuous(limits=c(0,1))+theme_minimal()


######
### Creating GD ~ GV biplots ### 
######
library(ggplot2)
library(patchwork)

GD_GV_biplots <- read_excel("GD.GV.biplots.xlsx")

pdf("GD.GV.biplot.pdf")
par(mfrow = c(2, 2))

PA<-ggplot(GD_GV_biplots, aes(x=pa_gd, y=pa_gv)) +
  geom_point() + ggtitle("P. angulosus") + xlab("Genomic Diversity") + ylab("Genomic Vulnerability") + theme_minimal()+theme(plot.title =element_text(face="italic"))

CP<-ggplot(GD_GV_biplots, aes(x=cp_gd, y=cp_gv)) +
  geom_point() + ggtitle("C. punctatus") + xlab("Genomic Diversity") + ylab("Genomic Vulnerability") + theme_minimal()+theme(plot.title =element_text(face="italic"))

SG<-ggplot(GD_GV_biplots, aes(x=sg_gd, y=sg_gv)) +
  geom_point() + ggtitle("S. granularis") + xlab("Genomic Diversity") + ylab("Genomic Vulnerability") + theme_minimal()+theme(plot.title =element_text(face="italic"))

ALL<-ggplot(GD_GV_biplots, aes(x=all_gd, y=all_gv)) +
  geom_point() + ggtitle("All species")+ xlab("Genomic Diversity") + ylab("Genomic Vulnerability")+ theme_minimal()

CP + PA + SG + ALL
dev.off()


### Getting % overlap between F, P, and ALL with endangered areas

setwd("~/Desktop/CH4_ms/prioritzr_sols/new_sols")

# make solution layers for F & P
######### p1 = Present_all

p1 <- problem(pu_data, pres.all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s1 <- solve(p1)

spplot(s1, "solution_1", col.regions = c("white", "darkgreen"), main = "P")

#save ALL_all selection freq as shp (to then rasterize)

solution_columns <- which(grepl("solution", names(s1)))
s1$selection_frequencies <- rowSums(as.matrix(s1@data[, solution_columns]))
P_all.sf <- as.data.frame(s1$selection_frequencies)

pu_P_sf <- cbind(pu_data, P_all.sf)
#pu_all_sf <- merge(pu_data, ALL_all.sf, by = "puid")

writeOGR(pu_P_sf, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols", layer = "pu_P_sf.shp",
         driver = "ESRI Shapefile" )

P_sf <-readOGR(".","pu_P_sf.shp")
ext  = extent(SA.ext <- extent(14, 34, -35.5, -25))
gridsize  = 0.008333333
r = raster(ext, res = gridsize)
P_sf.r <- rasterize(P_sf, field="s1_slc_", r)

writeRaster(P_sf.r, "P_all_sf.tif", format="GTiff", overwrite=TRUE)

P <- raster::stack("P_all_sf.tif")

#remove those with values less than 50
tmpfilter <- P < 50

#mask with <50 cells to get raster of solution
filtered_P <- mask(P, tmpfilter, maskvalue=1)

#turn into shapefile
Psol_shp = rasterToPolygons(filtered_P)

#write shapefile
writeOGR(Psol_shp, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols", layer = "Psol.shp",
         driver = "ESRI Shapefile" )


######### p2 = Future_all

p2 <- problem(pu_data, fut.all, cost_column = "cost") %>%
  add_min_set_objective()%>%
  add_relative_targets(0.20) %>% # 20% representation targets
  add_binary_decisions()%>% 
  add_gap_portfolio(number_solutions = 100, pool_gap = 0.2)

s2 <- solve(p2)

solution_columns <- which(grepl("solution", names(s2)))
s2$selection_frequencies <- rowSums(as.matrix(s2@data[, solution_columns]))
F_all.sf <- as.data.frame(s2$selection_frequencies)

pu_F_sf <- cbind(pu_data, F_all.sf)
#pu_all_sf <- merge(pu_data, ALL_all.sf, by = "puid")

writeOGR(pu_F_sf, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols", layer = "pu_F_sf.shp",
         driver = "ESRI Shapefile" )

F_sf <-readOGR(".","pu_F_sf.shp")
ext  = extent(SA.ext <- extent(14, 34, -35.5, -25))
gridsize  = 0.008333333
r = raster(ext, res = gridsize)
F_sf.r <- rasterize(F_sf, field="s2_slc_", r)

writeRaster(F_sf.r, "F_all_sf.tif", format="GTiff", overwrite=TRUE)

F <- raster::stack("F_all_sf.tif")

#remove those with values less than 50
tmpfilter <- F < 50

#mask with <50 cells to get raster of solution
filtered_F <- mask(F, tmpfilter, maskvalue=1)

#turn into shapefile
Fsol_shp = rasterToPolygons(filtered_F)

#write shapefile
writeOGR(Fsol_shp, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols", layer = "Fsol.shp",
         driver = "ESRI Shapefile" )


##### Getting % of solution overlapping with endangered habs 
end_ov_F <- read.csv("end_overlap_F.csv")
F_df <-end_ov_F[,2:3]
100*colSums(F_df[-1] != 0)/nrow(F_df)


##### Getting % of solution overlapping with MPAs
mpa_ov_F <- read.csv("mpa_overlay_F.csv")
mpa_F_df <-mpa_ov_F[,12:13]
colSums(mpa_F_df[-1] != 0)


##### Get % overlap per spp and GV/GD
> setwd("~/Desktop/CH4_ms/prioritzr_sols/new_sols")
> ov_gd_cp <- read.csv("GD_CP_ovlp.csv")
> View(ov_gd_cp)
> ov_gd_cp <-ov_gd_cp[,12:13]
> 100*colSums(ov_gd_cp[-1] >= 50)/nrow(ov_gd_cp)


############ Plotting raw input files
#first make pu layer an sf
pu_sf<-st_as_sf(pu_data)

#PC is different than above
PC <- raster::stack("all.PC.tif")
PC <- aggregate(PC, fact=3)

#crop each individ layer
CS_m <- mask(CS_all, pu_sf)
CV_m <- mask(CV_all, pu_sf)
HP_m <- mask(soc.all$soc.all_1, pu_sf)
PC_m <- mask(PC, pu_sf)
GD_m <- mask(GD_all, pu_sf)
GV_m <- mask(GV_all, pu_sf)

#for those that have multiple values, make >0 all one color
## Assign values, based on your condition
values(CS_m) <- as.numeric(values(CS_m) > 0)
values(CV_m) <- as.numeric(values(CV_m) > 0)
values(HP_m) <- as.numeric(values(HP_m) > 0)
values(PC_m) <- as.numeric(values(PC_m) > 0)
values(GD_m) <- as.numeric(values(GD_m) > 0)
values(GV_m) <- as.numeric(values(GV_m) > 0)

## Create a Color Function
cpal <- colorRampPalette(c("lightgrey", "blue"))

## Plot with raster-package
plot(CS_m, col=cpal(2))
plot(CV_m, col=cpal(2))
plot(HP_m, col=cpal(2))
plot(PC_m, col=cpal(2))
plot(GD_m, col=cpal(2))
plot(GV_m, col=cpal(2))

pdf("clim.input.maps.pdf")
plot(CS_m, col=cpal(2))
plot(CV_m, col=cpal(2))
dev.off()

pdf("soc.input.maps.pdf")
plot(HP_m, col=cpal(2))
plot(PC_m, col=cpal(2))
dev.off()

pdf("GD.input.maps.pdf")
plot(GD_m, col=cpal(2))
dev.off()


######## Create SF plots per scenario and plot with MPA overlay

#read SF
SF_new_0_2 <- read_excel("~/Desktop/CH4_ms/prioritzr_sols/new_sols/SF.new_0.2.xlsx")

pu_sfs <- cbind(pu_data, SF_new_0_2)

cols= c("floralwhite", "lightblue1", "deepskyblue", "darkblue")

Psf_mpa<-tm_shape(map4) + tm_polygons(col="white", border.col="grey")  + tm_shape(pu_sfs) + tm_polygons(col = "P_new", palette=cols) + tm_shape(mpas) + tm_polygons(col="white", border.col="green2", alpha = .5)
Fsf_mpa<-tm_shape(map4) + tm_polygons(col="white", border.col="grey")  + tm_shape(pu_sfs) + tm_polygons(col = "F_new", palette=cols) + tm_shape(mpas) + tm_polygons(col="white", border.col="green2", alpha = .5)

CSsf_mpa<-tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(pu_sfs) + tm_polygons(col = "CS_all_sf", palette=cols)  + tm_shape(mpas) + tm_polygons(col="white", border.col="green2", alpha = .5)
CVsf_mpa<-tm_shape(map4) + tm_polygons(col="white", border.col="grey")  + tm_shape(pu_sfs) + tm_polygons(col = "CV_all_sf", palette=cols)+ tm_shape(mpas) + tm_polygons(col="white", border.col="green2", alpha = .5)

HPsf_mpa<-tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(pu_sfs) + tm_polygons(col = "HF_new", palette=cols)  + tm_shape(mpas) + tm_polygons(col="white", border.col="green2", alpha = .5)
PCsf_mpa<-tm_shape(map4) + tm_polygons(col="white", border.col="grey")  + tm_shape(pu_sfs) + tm_polygons(col = "PC_new", palette=cols)+ tm_shape(mpas) + tm_polygons(col="white", border.col="green2", alpha = .5)

GDsf_mpa<-tm_shape(map4) + tm_polygons(col="white", border.col="grey") + tm_shape(pu_sfs) + tm_polygons(col = "GD_sf", palette=cols)  + tm_shape(mpas) + tm_polygons(col="white", border.col="green2", alpha = .5)
GVsf_mpa<-tm_shape(map4) + tm_polygons(col="white", border.col="grey")  + tm_shape(pu_sfs) + tm_polygons(col = "GV_sf", palette=cols)+ tm_shape(mpas) + tm_polygons(col="white", border.col="green2", alpha = .5)

pdf("scen.sfs.mpa.maps.pdf")
Psf_mpa
Fsf_mpa
CSsf_mpa
CVsf_mpa
HPsf_mpa
PCsf_mpa
GDsf_mpa
GVsf_mpa
dev.off()

######### Get percent overlap between ea. spp's GD/GV and all GD/GV

new_GV_all_sf <- read_excel("~/Desktop/CH4_ms/prioritzr_sols/new_sols/new_GV.all_sf.xlsx")

CP_GV_sf <- cbind(pu_data, new_GV_all_sf$CP_GV)
#pu_all_sf <- merge(pu_data, ALL_all.sf, by = "puid")

writeOGR(CP_GV_sf, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols/rebut_outputs", layer = "GV_CP_sf.shp",
         driver = "ESRI Shapefile" )

CP_GV_sf <-readOGR(".","GV_CP_sf.shp")
ext  = extent(SA.ext <- extent(14, 34, -35.5, -25))
gridsize  = 0.008333333
r = raster(ext, res = gridsize)
CP_GV_sf.r <- rasterize(CP_GV_sf, field="c_3__0_", r)

writeRaster(CP_GV_sf.r, "CP_GV_sf.r.tif", format="GTiff", overwrite=TRUE)

F <- raster::stack("CP_GV_sf.r.tif")

#remove those with values less than 50
tmpfilter <- F < 50

#mask with <50 cells to get raster of solution
filtered_F <- mask(F, tmpfilter, maskvalue=1)

#turn into shapefile
CP_GV_sol_shp = rasterToPolygons(filtered_F)

#write shapefile
writeOGR(CP_GV_sol_shp, dsn = "~/Desktop/CH4_ms/prioritzr_sols/new_sols/rebut_outputs", layer = "CP_GV.sol.shp",
         driver = "ESRI Shapefile" )


## making table of # of endangered/MPAs captured per scenario
# open data table (sheet1=endangered areas / sheet2=MPAs)
all_scen_overlap_table <- read_excel("all.scen.overlap.table.xlsx")

apply(all_scen_overlap_table,2,function(x) sum(x > 0))



