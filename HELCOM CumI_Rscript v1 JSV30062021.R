#Cumi script, version 1 by Josefine Egekvist, DTU Aqua, Denmark. Email: jsv@aqua.dtu.dk

library(rgdal)
library(dplyr)
library(sf)
library(sp)
library(ggplot2)
library(rgeos)
library(raster)
library(stars)
library(lwgeom)

st_or = function(x, y, dim = 2) {
  
  # st_erase to get the remainder of the intersection (taken from ?st_difference)
  st_erase = function(x, y) st_difference(x, st_union(st_combine(y)))
  
  # we need st_dump to extract polygons from a potential GEOMETRYCOLLECTION
  st_dump = function(x, dim) {
    dims = sapply(x, st_dimension)
    x[dims == dim, ]
  }
  
  # get overlap via intersection
  overlap = sf::st_intersection(x, y)
  
  # extract polygons (if dimm = 2)
  overlap = overlap[st_dimension(overlap) == dim, ]
  
  gc = which(sapply(seq(nrow(overlap)), function(i) {
    inherits(overlap[i, ], "GEOMETRYCOLLECTION")
  }))
  
  if (length(gc) > 0) {
    dmp = st_dump(overlap, dim = dim)
    overlap = rbind(overlap[-gc, ], dmp)
  }
  
  # get the non-intersecting parts and set missing attribute to NA
  #### I have modified this to to the erasure and insertion of NA
  #### values before the rbind to ensure they have the same columns
  diff1 = st_erase(x, y)
  diff2 = st_erase(y, x)
  
  diff1[, setdiff(names(overlap), names(diff1))] = NA
  diff2[, setdiff(names(overlap), names(diff2))] = NA
  
  cb = rbind(diff1,diff2)
  
  # return combined geometry set with attributes
  return(rbind(overlap, cb))
}

pathdir <- "Q:/dfad/users/joeg/home/VMS/200921_HELCOM_Cumi/Data_layers/Fishing_intensity"
setwd(pathdir)
#Cumi
#Read VMS by year, 1 km HELCOM grid
VMS2011 <- readOGR(dsn = paste(pathdir,"2011_identity_dissolve.shp",sep="/"))
VMS2012 <- readOGR(dsn = paste(pathdir,"2012_identity_dissolve.shp",sep="/"))
VMS2013 <- readOGR(dsn = paste(pathdir,"2013_identity_dissolve.shp",sep="/"))
VMS2014 <- readOGR(dsn = paste(pathdir,"2014_identity_dissolve.shp",sep="/"))
VMS2015 <- readOGR(dsn = paste(pathdir,"2015_identity_dissolve.shp",sep="/"))
VMS2016 <- readOGR(dsn = paste(pathdir,"2016_identity_dissolve.shp",sep="/"))

VMS2011_1 <- subset(VMS2011, select= -c(intensity, intensityn, i_2011, i_n_2011, t, t_2011, t_2012, i_2012, i_n_2012,t_2013, i_2013, i_n_2013, t_2014, i_2014, i_n_2014,
                                          t_2015, i_2015, i_n_2015,t_2016, i_2016, i_n_2016))

VMS2012_1 <- subset(VMS2012, select= -c(intensity, intensityn, i_2012, i_n_2012, t, t_2012))
VMS2013_1 <- subset(VMS2013, select= -c(intensity, intensityn, i_2013, i_n_2013, t, t_2013))
VMS2014_1 <- subset(VMS2014, select= -c(intensity, intensityn, i_2014, i_n_2014, t, t_2014))
VMS2015_1 <- subset(VMS2015, select= -c(intensity, intensityn, i_2015, i_n_2015, t, t_2015))
VMS2016_1 <- subset(VMS2016, select= -c(intensity, intensityn, i_2016, i_n_2016, t, t_2016))


#Clip with Danish EEZ
pathdir <- "Q:/gis/Dynamisk/GEOdata/BasicLayers/Boundaries/EEZ/MaritimeGr?nser"
EEZ <- readOGR(dsn = paste(pathdir,"Con_exe_tes_merge_dissolve.shp",sep="/"))
EEZ@proj4string
VMS_crs <- st_crs(VMS2011_DK)
EEZ_sf <- st_as_sf(EEZ)
EEZ_crs <- st_transform(EEZ_sf,VMS_crs)

VMS2011_sf <- st_as_sf(VMS2011_1)
VMS2011_DK <-st_intersection(VMS2011_sf,EEZ_crs)

VMS2012_sf <- st_as_sf(VMS2012_1)
VMS2012_DK <-st_intersection(VMS2012_sf,EEZ_crs)

VMS2013_sf <- st_as_sf(VMS2013_1)
VMS2013_DK <-st_intersection(VMS2013_sf,EEZ_crs)

VMS2014_sf <- st_as_sf(VMS2014_1)
VMS2014_DK <-st_intersection(VMS2014_sf,EEZ_crs)

VMS2015_sf <- st_as_sf(VMS2015_1)
VMS2015_DK <-st_intersection(VMS2015_sf,EEZ_crs)

VMS2016_sf <- st_as_sf(VMS2016_1)
VMS2016_DK <-st_intersection(VMS2016_sf,EEZ_crs)

pathdir <- "Q:/dfad/users/joeg/home/VMS/200921_HELCOM_Cumi/Data_layers/Fishing_intensity"
setwd(pathdir)

save(VMS2011_DK,file="VMS2011_DK.RData")  
save(VMS2012_DK,file="VMS2012_DK.RData")  
save(VMS2013_DK,file="VMS2013_DK.RData")  
save(VMS2014_DK,file="VMS2014_DK.RData")  
save(VMS2015_DK,file="VMS2015_DK.RData")  
save(VMS2016_DK,file="VMS2016_DK.RData")  

##############################################
pathdir <- "Q:/dfad/users/joeg/home/VMS/200921_HELCOM_Cumi/Data_layers/Fishing_intensity"
setwd(pathdir)
load("VMS2011_DK.RData")  
load("VMS2012_DK.RData")  
load("VMS2013_DK.RData")  
load("VMS2014_DK.RData")  
load("VMS2015_DK.RData")  
load("VMS2016_DK.RData")  

#Frequency

VMS2011_DK$freq_2011 <- cbind(VMS2011_DK$freq_2011, ifelse(VMS2011_DK$SUM_SurfSa>=0.05,1,0))
VMS2012_DK$freq_2012 <- cbind(VMS2012_DK$freq_2012, ifelse(VMS2012_DK$SUM_SurfSa>=0.05,1,0))
VMS2013_DK$freq_2013 <- cbind(VMS2013_DK$freq_2013, ifelse(VMS2013_DK$SUM_SurfSa>=0.05,1,0))
VMS2014_DK$freq_2014 <- cbind(VMS2014_DK$freq_2014, ifelse(VMS2014_DK$SUM_SurfSa>=0.05,1,0))
VMS2015_DK$freq_2015 <- cbind(VMS2015_DK$freq_2015, ifelse(VMS2015_DK$SUM_SurfSa>=0.05,1,0))
VMS2016_DK$freq_2016 <- cbind(VMS2016_DK$freq_2016, ifelse(VMS2016_DK$SUM_SurfSa>=0.05,1,0))

VMS2011_freq <- subset(VMS2011_DK, select= -c(SUM_SurfSa, SUM_SubSar, dis ))
VMS2012_freq <- subset(VMS2012_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2012, int_subs_2012, dis ))
VMS2013_freq <- subset(VMS2013_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2013, int_subs_2013, dis ))
VMS2014_freq <- subset(VMS2014_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2014, int_subs_2014, dis ))
VMS2015_freq <- subset(VMS2015_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2015, int_subs_2015, dis ))
VMS2016_freq <- subset(VMS2016_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2016, int_subs_2016, dis ))

st_geometry(VMS2012_freq) <- NULL # remove geometry, coerce to data.frame
st_geometry(VMS2013_freq) <- NULL # remove geometry, coerce to data.frame
st_geometry(VMS2014_freq) <- NULL # remove geometry, coerce to data.frame
st_geometry(VMS2015_freq) <- NULL # remove geometry, coerce to data.frame
st_geometry(VMS2016_freq) <- NULL # remove geometry, coerce to data.frame

freq1 <- merge(VMS2011_freq, VMS2012_freq, by = "FID_HELCOM")
freq2 <- merge(freq1, VMS2013_freq, by="FID_HELCOM")
freq3 <- merge(freq2, VMS2014_freq, by="FID_HELCOM")
freq4 <- merge(freq3, VMS2015_freq, by="FID_HELCOM")
freq5 <- merge(freq4, VMS2016_freq, by="FID_HELCOM")

freq5$freq_sum <- freq5$freq_2011 + freq5$freq_2012 + freq5$freq_2013 + freq5$freq_2014 + freq5$freq_2015 + freq5$freq_2016 

VMS2011t2016_freq <- subset(freq5, select= -c(freq_2011, EEZ.x, EEZ.y, freq_2012, EEZ.x.1, freq_2013, EEZ.y.1, freq_2014, EEZ.x.2, freq_2015, EEZ.y.2, freq_2016))
VMS2011t2016_freq$freq_cat[VMS2011t2016_freq$freq_sum >0 & VMS2011t2016_freq$freq_sum <= 3 ] <- "occasional"
VMS2011t2016_freq$freq_cat[VMS2011t2016_freq$freq_sum == 4 ] <- "regular"
VMS2011t2016_freq$freq_cat[VMS2011t2016_freq$freq_sum == 5 ] <- "frequent"
VMS2011t2016_freq$freq_cat[VMS2011t2016_freq$freq_sum == 6 ] <- "persistent"

#Intensity - as it is not specified how to combine the years, an average is made for surface og subsurface abrasion
VMS2011_DK$SurfSAR2011 <- VMS2011_DK$SUM_SurfSa
VMS2011_DK$SubSAR2011 <- VMS2011_DK$SUM_SubSar
VMS2011_int <- subset(VMS2011_DK, select= -c(SUM_SurfSa, SUM_SubSar,  freq_2011, EEZ, dis ))

VMS2012_DK$SurfSAR2012 <- VMS2012_DK$SUM_SurfSa
VMS2012_DK$SubSAR2012 <- VMS2012_DK$SUM_SubSar
VMS2012_int <- subset(VMS2012_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2012, int_subs_2012, freq_2012, EEZ, dis ))

VMS2013_DK$SurfSAR2013 <- VMS2013_DK$SUM_SurfSa
VMS2013_DK$SubSAR2013 <- VMS2013_DK$SUM_SubSar
VMS2013_int <- subset(VMS2013_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2013, int_subs_2013, freq_2013, EEZ, dis ))

VMS2014_DK$SurfSAR2014 <- VMS2014_DK$SUM_SurfSa
VMS2014_DK$SubSAR2014 <- VMS2014_DK$SUM_SubSar
VMS2014_int <- subset(VMS2014_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2014, int_subs_2014, freq_2014, EEZ, dis ))

VMS2015_DK$SurfSAR2015 <- VMS2015_DK$SUM_SurfSa
VMS2015_DK$SubSAR2015 <- VMS2015_DK$SUM_SubSar
VMS2015_int <- subset(VMS2015_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2015, int_subs_2015, freq_2015, EEZ, dis ))

VMS2016_DK$SurfSAR2016 <- VMS2016_DK$SUM_SurfSa
VMS2016_DK$SubSAR2016 <- VMS2016_DK$SUM_SubSar
VMS2016_int <- subset(VMS2016_DK, select= -c(SUM_SurfSa, SUM_SubSar, int_surf_2016, int_subs_2016, freq_2016, EEZ, dis ))

st_geometry(VMS2012_int) <- NULL # remove geometry, coerce to data.frame
st_geometry(VMS2013_int) <- NULL # remove geometry, coerce to data.frame
st_geometry(VMS2014_int) <- NULL # remove geometry, coerce to data.frame
st_geometry(VMS2015_int) <- NULL # remove geometry, coerce to data.frame
st_geometry(VMS2016_int) <- NULL # remove geometry, coerce to data.frame

int1 <- merge(VMS2011_int, VMS2012_int, by = "FID_HELCOM")
int2 <- merge(int1, VMS2013_int, by="FID_HELCOM")
int3 <- merge(int2, VMS2014_int, by="FID_HELCOM")
int4 <- merge(int3, VMS2015_int, by="FID_HELCOM")
int5 <- merge(int4, VMS2016_int, by="FID_HELCOM")

int5$SurfSAR_mean <- (int5$SurfSAR2011 +int5$SurfSAR2012 + int5$SurfSAR2013 + int5$SurfSAR2014 + int5$SurfSAR2015 + int5$SurfSAR2016)/6
int5$SubSAR_mean <- (int5$SubSAR2011 + int5$SubSAR2012 + int5$SubSAR2013 + int5$SubSAR2014 + int5$SubSAR2015 + int5$SubSAR2016)/6

VMS2011t2016_int <- subset(int5, select= -c(SurfSAR2011,SubSAR2011,SurfSAR2012,SubSAR2012,SurfSAR2013,SubSAR2013,SurfSAR2014,SubSAR2014,SurfSAR2015,SubSAR2015,SurfSAR2016,SubSAR2016))

st_write(int5,"Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Fishing_intensity\\VMS2011t2016_int2.shp", delete_dsn = T)

VMS2011t2016_int$SurfInt[VMS2011t2016_int$SurfSAR_mean>=0.05 & VMS2011t2016_int$SurfSAR_mean<0.33] <- "very low"
VMS2011t2016_int$SurfInt[VMS2011t2016_int$SurfSAR_mean>=0.33 & VMS2011t2016_int$SurfSAR_mean<0.66] <- "low"
VMS2011t2016_int$SurfInt[VMS2011t2016_int$SurfSAR_mean>=0.66 & VMS2011t2016_int$SurfSAR_mean<2] <- "moderate"
VMS2011t2016_int$SurfInt[VMS2011t2016_int$SurfSAR_mean>=2 ] <- "high"

VMS2011t2016_int$SubInt[VMS2011t2016_int$SubSAR_mean>=0.05 & VMS2011t2016_int$SubSAR_mean<0.33] <- "very low"
VMS2011t2016_int$SubInt[VMS2011t2016_int$SubSAR_mean>=0.33 & VMS2011t2016_int$SubSAR_mean<0.66] <- "low"
VMS2011t2016_int$SubInt[VMS2011t2016_int$SubSAR_mean>=0.66 & VMS2011t2016_int$SubSAR_mean<2] <- "moderate"
VMS2011t2016_int$SubInt[VMS2011t2016_int$SubSAR_mean>=2 ] <- "high"

st_geometry(VMS2011t2016_int) <- NULL # remove geometry, coerce to data.frame
FishingPressure <- merge(VMS2011t2016_freq, VMS2011t2016_int, by = "FID_HELCOM")

#Surface
FishingPressure$Table2_surface[FishingPressure$freq_cat=="persistent" & FishingPressure$SurfInt=="high"] <- "High"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="persistent" & FishingPressure$SurfInt=="moderate"] <- "High"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="persistent" & FishingPressure$SurfInt=="low"] <- "Moderate"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="persistent" & FishingPressure$SurfInt=="very low"] <- "Moderate"

FishingPressure$Table2_surface[FishingPressure$freq_cat=="frequent" & FishingPressure$SurfInt=="high"] <- "High"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="frequent" & FishingPressure$SurfInt=="moderate"] <- "Moderate"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="frequent" & FishingPressure$SurfInt=="low"] <- "Moderate"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="frequent" & FishingPressure$SurfInt=="very low"] <- "Low"

FishingPressure$Table2_surface[FishingPressure$freq_cat=="regular" & FishingPressure$SurfInt=="high"] <- "Moderate"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="regular" & FishingPressure$SurfInt=="moderate"] <- "Moderate"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="regular" & FishingPressure$SurfInt=="low"] <- "Low"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="regular" & FishingPressure$SurfInt=="very low"] <- "Very low"

FishingPressure$Table2_surface[FishingPressure$freq_cat=="occasional" & FishingPressure$SurfInt=="high"] <- "Moderate"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="occasional" & FishingPressure$SurfInt=="moderate"] <- "Low"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="occasional" & FishingPressure$SurfInt=="low"] <- "Very low"
FishingPressure$Table2_surface[FishingPressure$freq_cat=="occasional" & FishingPressure$SurfInt=="very low"] <- "Very low"

#Subsurface
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="persistent" & FishingPressure$SubInt=="high"] <- "High"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="persistent" & FishingPressure$SubInt=="moderate"] <- "High"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="persistent" & FishingPressure$SubInt=="low"] <- "Moderate"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="persistent" & FishingPressure$SubInt=="very low"] <- "Moderate"

FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="frequent" & FishingPressure$SubInt=="high"] <- "High"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="frequent" & FishingPressure$SubInt=="moderate"] <- "Moderate"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="frequent" & FishingPressure$SubInt=="low"] <- "Moderate"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="frequent" & FishingPressure$SubInt=="very low"] <- "Low"

FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="regular" & FishingPressure$SubInt=="high"] <- "Moderate"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="regular" & FishingPressure$SubInt=="moderate"] <- "Moderate"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="regular" & FishingPressure$SubInt=="low"] <- "Low"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="regular" & FishingPressure$SubInt=="very low"] <- "Very low"

FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="occasional" & FishingPressure$SubInt=="high"] <- "Moderate"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="occasional" & FishingPressure$SubInt=="moderate"] <- "Low"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="occasional" & FishingPressure$SubInt=="low"] <- "Very low"
FishingPressure$Table2_subsurface[FishingPressure$freq_cat=="occasional" & FishingPressure$SubInt=="very low"] <- "Very low"

st_write(FishingPressure,"Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Fishing_intensity\\VMS2011t2016_freq_int.shp", delete_dsn = T)

#Habitat

pathdir <- "Q:/dfad/users/joeg/home/VMS/200921_HELCOM_Cumi/Data_layers/Habitat"
setwd(pathdir)
habitat <- readOGR(dsn = paste(pathdir,"Broad_scale_habitats.shp",sep="/"))
habitat_sf <- st_as_sf(habitat)
habitat_sf$habitatComb[!is.na(habitat_sf$Habitat)] <- "Circalittoral hard substrate"
habitat_sf$habitatComb[!is.na(habitat_sf$Habitat_1)] <- "Circalittoral mud"
habitat_sf$habitatComb[!is.na(habitat_sf$Habitat_12)] <- "Circalittoral sand"
habitat_sf$habitatComb[!is.na(habitat_sf$Habitat_13)] <- "Infralittoral hard substrate"
habitat_sf$habitatComb[!is.na(habitat_sf$Habitat_14)] <- "Infralittoral mixed sediments"
habitat_sf$habitatComb[!is.na(habitat_sf$Habitat_15)] <- "Infralittoral mud"
habitat_sf$habitatComb[!is.na(habitat_sf$Habitat_16)] <- "Infralittoral sand"
habitat_sf$habitatComb[!is.na(habitat_sf$Habitat_17)] <- "Circalittoral mixed sediments"

habitat_sf1 <- subset(habitat_sf, select= -c(Habitat, Habitat_1, Habitat_12, Habitat_13, Habitat_14, Habitat_15, Habitat_16, Habitat_17,GRIDCODE, GRIDCODE_1, GRIDCODE_2,GRIDCODE_3, GRIDCODE_4, GRIDCODE_5, GRIDCODE_6, GRIDCODE_7,
                                             Salinity, FID_circa1, ID_1, Salinity_1, FID_circ_1, ID_12, Salinity_2, FID_infr_h, ID_12_13, 
                                             Salinity_3, FID_infr_m, ID_12_1_14, Salinity_4, FID_infr_1, ID_12_1_15, Salinity_5, FID_infr_s, 
                                             ID_12_1_16, Salinity_6, FID_circ_2, ID_12_1_17, Salinity_7))
#Circalittoral hard substrate
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral hard substrate" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral hard substrate" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral hard substrate" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral hard substrate" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Low"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral hard substrate" & habitat_sf1$SalinityCo=="Polyhaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral hard substrate" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral hard substrate" & habitat_sf1$SalinityCo=="Euhaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral hard substrate" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"

#Infralittoral hard substrate
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral hard substrate" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral hard substrate" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral hard substrate" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral hard substrate" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Low"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral hard substrate" & habitat_sf1$SalinityCo=="Polyhaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral hard substrate" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral hard substrate" & habitat_sf1$SalinityCo=="Euhaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral hard substrate" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"

#Circalittoral sand
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral sand" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral sand" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral sand" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral sand" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral sand" & habitat_sf1$SalinityCo=="Polyhaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral sand" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral sand" & habitat_sf1$SalinityCo=="Euhaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral sand" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"

#Infralittoral sand
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral sand" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral sand" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral sand" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral sand" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral sand" & habitat_sf1$SalinityCo=="Polyhaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral sand" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral sand" & habitat_sf1$SalinityCo=="Euhaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral sand" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"

#Circalittoral mud
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral mud" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral mud" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral mud" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral mud" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral mud" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral mud" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral mud" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral mud" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"

#Infralittoral mud
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral mud" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral mud" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral mud" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral mud" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral mud" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral mud" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral mud" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral mud" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"

#Circalittoral mixed sediments
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral mixed sediments" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral mixed sediments" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral mixed sediments" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral mixed sediments" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Low"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral mixed sediments" & habitat_sf1$SalinityCo=="Polyhaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral mixed sediments" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Circalittoral mixed sediments" & habitat_sf1$SalinityCo=="Euhaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Circalittoral mixed sediments" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"

#Infralittoral mixed sediments
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral mixed sediments" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral mixed sediments" & habitat_sf1$SalinityCo=="Oligohaline"] <- "Moderate"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral mixed sediments" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral mixed sediments" & habitat_sf1$SalinityCo=="Mesohaline"] <- "Low"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral mixed sediments" & habitat_sf1$SalinityCo=="Polyhaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral mixed sediments" & habitat_sf1$SalinityCo=="Polyhaline"] <- "High"
habitat_sf1$sensitivity_subsurface[habitat_sf1$habitatComb=="Infralittoral mixed sediments" & habitat_sf1$SalinityCo=="Euhaline"] <- "Low"
habitat_sf1$sensitivity_Surface[habitat_sf1$habitatComb=="Infralittoral mixed sediments" & habitat_sf1$SalinityCo=="Euhaline"] <- "High"

#Dansk EEZ
habitat_crs <- st_transform(habitat_sf1,crs(EEZ_crs))
habitat_crs1 <- st_make_valid(habitat_crs)
EEZ_crs1 <- st_make_valid(EEZ_crs)
habitat_DK <-st_intersection(habitat_crs1,EEZ_crs1)

st_write(habitat_DK,"Q:/dfad/users/joeg/home/VMS/200921_HELCOM_Cumi/Results/CumI_Habitat_sensitivity_DK.shp", delete_dsn = T)
st_write(habitat_sf1,"Q:/dfad/users/joeg/home/VMS/200921_HELCOM_Cumi/Results/CumI_Habitat_sensitivity.shp", delete_dsn = T)


#Join Habitat og FishingPressure
FishingPressure1 <- st_make_valid(FishingPressure)
habitat_DK1 <- st_make_valid(habitat_DK)
FishingPressure_Habitat <- st_join(FishingPressure1,habitat_DK1, left=TRUE, largest=TRUE)

#Impact - table 4 surface
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="High" & FishingPressure_Habitat$sensitivity_Surface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="High" & FishingPressure_Habitat$sensitivity_Surface=="Moderate"] <- "High"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="High" & FishingPressure_Habitat$sensitivity_Surface=="Low"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="High" & FishingPressure_Habitat$sensitivity_Surface=="Very low"] <- "m1"

FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Moderate" & FishingPressure_Habitat$sensitivity_Surface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Moderate" & FishingPressure_Habitat$sensitivity_Surface=="Moderate"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Moderate" & FishingPressure_Habitat$sensitivity_Surface=="Low"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Moderate" & FishingPressure_Habitat$sensitivity_Surface=="Very low"] <- "Low"

FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Low" & FishingPressure_Habitat$sensitivity_Surface=="High"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Low" & FishingPressure_Habitat$sensitivity_Surface=="Moderate"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Low" & FishingPressure_Habitat$sensitivity_Surface=="Low"] <- "Low"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Low" & FishingPressure_Habitat$sensitivity_Surface=="Very low"] <- "Very low"

FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Very low" & FishingPressure_Habitat$sensitivity_Surface=="High"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Very low" & FishingPressure_Habitat$sensitivity_Surface=="Moderate"] <- "Low"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Very low" & FishingPressure_Habitat$sensitivity_Surface=="Low"] <- "Very low"
FishingPressure_Habitat$Table4_FishingImpactSurface[FishingPressure_Habitat$Table2_surface=="Very low" & FishingPressure_Habitat$sensitivity_Surface=="Very low"] <- "Very low"


#Impact - table 4 subsurface
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface =="High" & FishingPressure_Habitat$sensitivity_subsurface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="High" & FishingPressure_Habitat$sensitivity_subsurface=="Moderate"] <- "High"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="High" & FishingPressure_Habitat$sensitivity_subsurface=="Low"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="High" & FishingPressure_Habitat$sensitivity_subsurface=="Very low"] <- "m1"

FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Moderate" & FishingPressure_Habitat$sensitivity_subsurface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Moderate" & FishingPressure_Habitat$sensitivity_subsurface=="Moderate"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Moderate" & FishingPressure_Habitat$sensitivity_subsurface=="Low"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Moderate" & FishingPressure_Habitat$sensitivity_subsurface=="Very low"] <- "Low"

FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Low" & FishingPressure_Habitat$sensitivity_subsurface=="High"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Low" & FishingPressure_Habitat$sensitivity_subsurface=="Moderate"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Low" & FishingPressure_Habitat$sensitivity_subsurface=="Low"] <- "Low"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Low" & FishingPressure_Habitat$sensitivity_subsurface=="Very low"] <- "Very low"

FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Very low" & FishingPressure_Habitat$sensitivity_subsurface=="High"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Very low" & FishingPressure_Habitat$sensitivity_subsurface=="Moderate"] <- "Low"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Very low" & FishingPressure_Habitat$sensitivity_subsurface=="Low"] <- "Very low"
FishingPressure_Habitat$Table4_FishingImpactSubSurface[FishingPressure_Habitat$Table2_subsurface=="Very low" & FishingPressure_Habitat$sensitivity_subsurface=="Very low"] <- "Very low"

#Surface and subsurface impact - take maximum
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Very low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Very low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m3"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Very low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m2"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Very low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m1"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Very low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Low"] <- "Low"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Very low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Very low"] <- "Very low"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Very low" & is.na(FishingPressure_Habitat$Table4_FishingImpactSubSurface)] <- "Very low"

FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m3"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m2"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m1"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Low"] <- "Low"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Low" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Very low"] <- "Low"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="Low" & is.na(FishingPressure_Habitat$Table4_FishingImpactSubSurface)] <- "Low"

FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m1" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m1" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m3"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m1" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m2"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m1" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m1"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m1" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Low"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m1" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Very low"] <- "m1"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m1" & is.na(FishingPressure_Habitat$Table4_FishingImpactSubSurface)] <- "m1"

FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m2" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m2" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m3"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m2" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m2"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m2" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m1"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m2" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Low"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m2" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Very low"] <- "m2"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m2" & is.na(FishingPressure_Habitat$Table4_FishingImpactSubSurface)] <- "m2"

FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m3" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m3" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m3"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m3" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m2"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m3" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m1"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m3" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Low"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m3" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Very low"] <- "m3"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="m3" & is.na(FishingPressure_Habitat$Table4_FishingImpactSubSurface)] <- "m3"

FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="High" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="High"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="High" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m3"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="High" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m2"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="High" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="m1"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="High" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Low"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="High" & FishingPressure_Habitat$Table4_FishingImpactSubSurface=="Very low"] <- "High"
FishingPressure_Habitat$Table4_FishingImpact[FishingPressure_Habitat$Table4_FishingImpactSurface=="High" & is.na(FishingPressure_Habitat$Table4_FishingImpactSubSurface)] <- "High"

st_write(FishingPressure_Habitat,"Q:/dfad/users/joeg/home/VMS/200921_HELCOM_Cumi/Results/CumI_FishingImpact.shp", delete_dsn = T)

#Joiner with HELCOM Subdivisions. 
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\HELCOM_subbasins_2018"
setwd(pathdir)
Subbasins <- readOGR(dsn = paste(pathdir,"HELCOM_subbasins_2018.shp",sep="/"))
Subbasins_sf <- st_as_sf(Subbasins)
Subbasins_crs <- st_transform(Subbasins_sf,crs(FishingPressure_Habitat))
Subbasins_crs1 <- st_make_valid(Subbasins_crs)

FishingImpact_subbasins <- st_join(FishingPressure_Habitat,Subbasins_crs1, left=TRUE, largest=TRUE)

FishingImpact_subbasins$area_km2 <- st_area(FishingImpact_subbasins)/1000000

FishingImpact_subbasins_sum <- FishingImpact_subbasins %>%
  group_by(level_2, habitatComb, Table4_FishingImpact) %>%
  summarise(area_km2 = sum(area_km2))
st_geometry(FishingImpact_subbasins_sum) <- NULL

write.csv(FishingImpact_subbasins_sum, "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Grafer og tabeller\\FishingImpact_subbasins_sum.csv")

#Other pressure layers
#Cables under construction
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Cables"
setwd(pathdir)
cables_buffer <- readOGR(dsn = paste(pathdir,"cables_under_construction_buffer.shp",sep="/"))

cables_buffer_sf <- st_as_sf(cables_buffer)
cables_buffer_crs <- st_transform(cables_buffer_sf,crs(FishingImpact_subbasins))

FishingImpact_p2 <- st_join(FishingImpact_subbasins,cables_buffer_crs, left=TRUE, largest=TRUE)

#Coastal defence under construction
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Coastal_defense"
setwd(pathdir)
Coastal_defense_buffer <- readOGR(dsn = paste(pathdir,"Coastal_defense_under_construction_buffer.shp",sep="/"))

Coastal_defense_sf <- st_as_sf(Coastal_defense_buffer)
Coastal_defense_crs <- st_transform(Coastal_defense_sf,crs(FishingImpact_p2))

FishingImpact_p2 <- st_join(FishingImpact_p2,Coastal_defense_crs, left=TRUE, largest=TRUE)

#Depositing areas 06-18
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Depositing_areas_06_18"
setwd(pathdir)
Depositing_areas_buffer <- readOGR(dsn = paste(pathdir,"Depositing_areas_06_18_buffer.shp",sep="/"))

Depositing_areas_sf <- st_as_sf(Depositing_areas_buffer)
Depositing_areas_crs <- st_transform(Depositing_areas_sf,crs(FishingImpact_p2))

FishingImpact_p3 <- st_join(FishingImpact_p2,Depositing_areas_crs, left=TRUE, largest=TRUE)

#Depositing points 06-18
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Depositing_points_06_18"
setwd(pathdir)
Depositing_points_buffer <- readOGR(dsn = paste(pathdir,"Depositing_points_06_18_buffer.shp",sep="/"))

Depositing_points_sf <- st_as_sf(Depositing_points_buffer)
Depositing_points_crs <- st_transform(Depositing_points_sf,crs(FishingImpact_p3))

FishingImpact_p4 <- st_join(FishingImpact_p3,Depositing_points_crs, left=TRUE, largest=TRUE)

#Extraction of sand and gravel
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Extraction_of_sand_and_gravel"
setwd(pathdir)
Extraction_buffer <- readOGR(dsn = paste(pathdir,"Extraction_of_sand_and_gravel_buffer.shp",sep="/"))

Extraction_sf <- st_as_sf(Extraction_buffer)
Extraction_crs <- st_transform(Extraction_sf,crs(FishingImpact_p4))

FishingImpact_p5 <- st_join(FishingImpact_p4,Extraction_crs, left=TRUE, largest=TRUE)

#Finfish mariculture
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Finfish_mariculture"
setwd(pathdir)
Finfish_mariculture <- readOGR(dsn = paste(pathdir,"Finfish_mariculture_buffer.shp",sep="/"))

Finfish_mariculture_sf <- st_as_sf(Finfish_mariculture)
Finfish_mariculture_crs <- st_transform(Finfish_mariculture_sf,crs(FishingImpact_p5))

FishingImpact_p6 <- st_join(FishingImpact_p5,Finfish_mariculture_crs, left=TRUE, largest=TRUE)

#Pipelines
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Pipelines"
setwd(pathdir)
Pipelines <- readOGR(dsn = paste(pathdir,"Pipelines_in_operation_buffer.shp",sep="/"))

Pipelines_sf <- st_as_sf(Pipelines)
Pipelines_crs <- st_transform(Pipelines_sf,crs(FishingImpact_p6))

FishingImpact_p7 <- st_join(FishingImpact_p6,Pipelines_crs, left=TRUE, largest=TRUE)

#wind_turbines - wind farms under construction
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\wind_turbines"
setwd(pathdir)
Windfarms_construction <- readOGR(dsn = paste(pathdir,"wind_farms_under_construction_buffer.shp",sep="/"))

Windfarms_construction_sf <- st_as_sf(Windfarms_construction)
Windfarms_construction_crs <- st_transform(Windfarms_construction_sf,crs(FishingImpact_p7))

FishingImpact_p8 <- st_join(FishingImpact_p7,Windfarms_construction_crs, left=TRUE, largest=TRUE)

#wind_turbines - wind farms in operation
pathdir <- "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\wind_turbines"
setwd(pathdir)
Windfarms_operation <- readOGR(dsn = paste(pathdir,"wind_farms_in_operation_buffer.shp",sep="/"))

Windfarms_operation_sf <- st_as_sf(Windfarms_operation)
Windfarms_operation_crs <- st_transform(Windfarms_operation_sf,crs(FishingImpact_p8))

FishingImpact_p9 <- st_join(FishingImpact_p8,Windfarms_operation_crs, left=TRUE, largest=TRUE)

#Depth
Dybde <- raster("Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Depthreliefmap\\Depth_Total_Dynocs1.tif")
dybde_crs <- projectRaster(Dybde, crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
#AIS
AIS <- raster("Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Shipping (AIS-density)\\2011-2015_all_ship_types_MEAN\\2011-2015_all_ship_types_MEAN.tif")
AIS_crs <- projectRaster(AIS, crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
AIS_poly <- rasterToPolygons(AIS_crs)
AIS_poly_sf <- st_as_sf(AIS_poly)
AIS_DK <-st_intersection(AIS_poly_sf,EEZ_crs1)
st_write(AIS_DK,"Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Shipping (AIS-density)\\AIS_DK.shp", delete_dsn = T) 

AIS_dybde <- extract(dybde_crs, AIS_DK, fun=mean, sp=TRUE)

AIS_dybde$DybdeRev <- AIS_dybde$Depth_Total_Dynocs1*-1
AIS_dybde$DybdeRev[is.na(AIS_dybde$DybdeRev)] <- 0
AIS_dybde$shippingIntensity <- 0
AIS_dybde$shippingIntensity <- ifelse(AIS_dybde$DybdeRev>0 & AIS_dybde$DybdeRev <=10,AIS_dybde$X2011_2*1,AIS_dybde$shippingIntensity)
AIS_dybde$shippingIntensity <- ifelse(AIS_dybde$DybdeRev>10 & AIS_dybde$DybdeRev <=15,AIS_dybde$X2011_2*0.5,AIS_dybde$shippingIntensity)
AIS_dybde$shippingIntensity <- ifelse(AIS_dybde$DybdeRev>15 & AIS_dybde$DybdeRev <=20,AIS_dybde$X2011_2*0.25,AIS_dybde$shippingIntensity)
AIS_dybde$shippingIntensity <- ifelse(AIS_dybde$DybdeRev>20 & AIS_dybde$DybdeRev <=25,AIS_dybde$X2011_2*0.1,AIS_dybde$shippingIntensity)

AIS_dybde$test <- AIS_dybde$X2011.2015_all_ship_types_MEAN*0.5

test <- AIS_dybde[AIS_dybde$DybdeRev>10 & AIS_dybde$DybdeRev <=15,]
test_sf <- st_as_sf(test)

ranges <- 1
ranges$max <- max(AIS_dybde$shippingIntensity, na.rm=TRUE)
ranges$pct_75 <- ranges$max*0.75
ranges$pct_50 <- ranges$max*0.5
ranges$pct_25 <- ranges$max*0.25

AIS_dybde$shippingCategory[AIS_dybde$shippingIntensity>= 10000 ] <- "high"
AIS_dybde$shippingCategory[AIS_dybde$shippingIntensity>= 5000 & AIS_dybde$shippingIntensity< 10000] <- "moderate"
AIS_dybde$shippingCategory[AIS_dybde$shippingIntensity>= 100 & AIS_dybde$shippingIntensity< 5000] <- "low"
AIS_dybde$shippingCategory[AIS_dybde$shippingIntensity>= 0 & AIS_dybde$shippingIntensity< 100] <- "very low"
AIS_dybde_sf <- st_as_sf(AIS_dybde)

AIS_dybde_crs <- st_transform(AIS_dybde_sf,crs(FishingImpact_p9))
AIS_dybde_crs1 <- st_make_valid(AIS_dybde_crs)

FishingImpact_p10 <- st_join(FishingImpact_p9,AIS_dybde_crs1, left=TRUE, largest=TRUE)

###############################################
#Joining impacts from layers
#Fishing impact + Cables under construction
Cum_impacts1 <- FishingImpact_p10
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="High" & Cum_impacts1$cc_impact=="very low"] <- "High"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="High" & Cum_impacts1$cc_impact=="low"] <- "High"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="High" & Cum_impacts1$cc_impact=="moderate"] <- "High"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m3" & Cum_impacts1$cc_impact=="very low"] <- "m3"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m3" & Cum_impacts1$cc_impact=="low"] <- "m3"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m3" & Cum_impacts1$cc_impact=="moderate"] <- "m3"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m2" & Cum_impacts1$cc_impact=="very low"] <- "m2"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m2" & Cum_impacts1$cc_impact=="low"] <- "m2"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m2" & Cum_impacts1$cc_impact=="moderate"] <- "m2"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m1" & Cum_impacts1$cc_impact=="very low"] <- "m1"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m1" & Cum_impacts1$cc_impact=="low"] <- "m1"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m1" & Cum_impacts1$cc_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="Low" & Cum_impacts1$cc_impact=="very low"] <- "Low"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="Low" & Cum_impacts1$cc_impact=="low"] <- "Low"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="Low" & Cum_impacts1$cc_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="Very low" & Cum_impacts1$cc_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="Very low" & Cum_impacts1$cc_impact=="low"] <- "Low"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="Very low" & Cum_impacts1$cc_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="NA" & Cum_impacts1$cc_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="NA" & Cum_impacts1$cc_impact=="low"] <- "Low"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="NA" & Cum_impacts1$cc_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact1[is.na(Cum_impacts1$Table4_FishingImpact) & Cum_impacts1$cc_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact1[is.na(Cum_impacts1$Table4_FishingImpact) & Cum_impacts1$cc_impact=="low"] <- "Low"
Cum_impacts1$cumImpact1[is.na(Cum_impacts1$Table4_FishingImpact) & Cum_impacts1$cc_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="High" & is.na(Cum_impacts1$cc_impact)] <- "High"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m3" & is.na(Cum_impacts1$cc_impact)] <- "m3"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m2" & is.na(Cum_impacts1$cc_impact)] <- "m2"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="m1" & is.na(Cum_impacts1$cc_impact)] <- "m1"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="Low" & is.na(Cum_impacts1$cc_impact)] <- "Low"
Cum_impacts1$cumImpact1[Cum_impacts1$Table4_FishingImpact=="Very low" & is.na(Cum_impacts1$cc_impact)] <- "Very low"

Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="High" & Cum_impacts1$cd_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="High" & Cum_impacts1$cd_impact=="moderate"] <- "High"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="High" & Cum_impacts1$cd_impact=="low"] <- "High"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="High" & Cum_impacts1$cd_impact=="very low"] <- "High"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m3" & Cum_impacts1$cd_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m3" & Cum_impacts1$cd_impact=="moderate"] <- "m3"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m3" & Cum_impacts1$cd_impact=="low"] <- "m3"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m3" & Cum_impacts1$cd_impact=="very low"] <- "m3"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m2" & Cum_impacts1$cd_impact=="high"] <- "High"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m2" & Cum_impacts1$cd_impact=="moderate"] <- "m2"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m2" & Cum_impacts1$cd_impact=="low"] <- "m2"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m2" & Cum_impacts1$cd_impact=="very low"] <- "m2"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m1" & Cum_impacts1$cd_impact=="high"] <- "High"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m1" & Cum_impacts1$cd_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m1" & Cum_impacts1$cd_impact=="low"] <- "m1"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m1" & Cum_impacts1$cd_impact=="very low"] <- "m1"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Low" & Cum_impacts1$cd_impact=="high"] <- "High"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Low" & Cum_impacts1$cd_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Low" & Cum_impacts1$cd_impact=="low"] <- "Low"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Low" & Cum_impacts1$cd_impact=="very low"] <- "Low"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Very low" & Cum_impacts1$cd_impact=="high"] <- "High"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Very low" & Cum_impacts1$cd_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Very low" & Cum_impacts1$cd_impact=="low"] <- "Low"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Very low" & Cum_impacts1$cd_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="NA" & Cum_impacts1$cd_impact=="high"] <- "High"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="NA" & Cum_impacts1$cd_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="NA" & Cum_impacts1$cd_impact=="low"] <- "Low"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="NA" & Cum_impacts1$cd_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact2[is.na(Cum_impacts1$cumImpact1) & Cum_impacts1$cd_impact=="high"] <- "High"
Cum_impacts1$cumImpact2[is.na(Cum_impacts1$cumImpact1) & Cum_impacts1$cd_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact2[is.na(Cum_impacts1$cumImpact1) & Cum_impacts1$cd_impact=="low"] <- "Low"
Cum_impacts1$cumImpact2[is.na(Cum_impacts1$cumImpact1) & Cum_impacts1$cd_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="High" & is.na(Cum_impacts1$cd_impact)] <- "High"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m3" & is.na(Cum_impacts1$cd_impact)] <- "m3"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m2" & is.na(Cum_impacts1$cd_impact)] <- "m2"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="m1" & is.na(Cum_impacts1$cd_impact)] <- "m1"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Low" & is.na(Cum_impacts1$cd_impact)] <- "Low"
Cum_impacts1$cumImpact2[Cum_impacts1$cumImpact1=="Very low" & is.na(Cum_impacts1$cd_impact)] <- "Very low"

Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="High" & Cum_impacts1$da_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="High" & Cum_impacts1$da_impact=="moderate"] <- "High"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="High" & Cum_impacts1$da_impact=="low"] <- "High"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="High" & Cum_impacts1$da_impact=="very low"] <- "High"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m3" & Cum_impacts1$da_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m3" & Cum_impacts1$da_impact=="moderate"] <- "m3"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m3" & Cum_impacts1$da_impact=="low"] <- "m3"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m3" & Cum_impacts1$da_impact=="very low"] <- "m3"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m2" & Cum_impacts1$da_impact=="high"] <- "High"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m2" & Cum_impacts1$da_impact=="moderate"] <- "m2"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m2" & Cum_impacts1$da_impact=="low"] <- "m2"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m2" & Cum_impacts1$da_impact=="very low"] <- "m2"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m1" & Cum_impacts1$da_impact=="high"] <- "High"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m1" & Cum_impacts1$da_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m1" & Cum_impacts1$da_impact=="low"] <- "m1"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m1" & Cum_impacts1$da_impact=="very low"] <- "m1"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Low" & Cum_impacts1$da_impact=="high"] <- "High"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Low" & Cum_impacts1$da_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Low" & Cum_impacts1$da_impact=="low"] <- "Low"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Low" & Cum_impacts1$da_impact=="very low"] <- "Low"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Very low" & Cum_impacts1$da_impact=="high"] <- "High"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Very low" & Cum_impacts1$da_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Very low" & Cum_impacts1$da_impact=="low"] <- "Low"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Very low" & Cum_impacts1$da_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="NA" & Cum_impacts1$da_impact=="high"] <- "High"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="NA" & Cum_impacts1$da_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="NA" & Cum_impacts1$da_impact=="low"] <- "Low"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="NA" & Cum_impacts1$da_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact3[is.na(Cum_impacts1$cumImpact2) & Cum_impacts1$da_impact=="high"] <- "High"
Cum_impacts1$cumImpact3[is.na(Cum_impacts1$cumImpact2) & Cum_impacts1$da_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact3[is.na(Cum_impacts1$cumImpact2) & Cum_impacts1$da_impact=="low"] <- "Low"
Cum_impacts1$cumImpact3[is.na(Cum_impacts1$cumImpact2) & Cum_impacts1$da_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="High" & is.na(Cum_impacts1$da_impact)] <- "High"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m3" & is.na(Cum_impacts1$da_impact)] <- "m3"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m2" & is.na(Cum_impacts1$da_impact)] <- "m2"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="m1" & is.na(Cum_impacts1$da_impact)] <- "m1"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Low" & is.na(Cum_impacts1$da_impact)] <- "Low"
Cum_impacts1$cumImpact3[Cum_impacts1$cumImpact2=="Very low" & is.na(Cum_impacts1$da_impact)] <- "Very low"

Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="High" & Cum_impacts1$dp_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="High" & Cum_impacts1$dp_impact=="moderate"] <- "High"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="High" & Cum_impacts1$dp_impact=="low"] <- "High"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="High" & Cum_impacts1$dp_impact=="very low"] <- "High"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m3" & Cum_impacts1$dp_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m3" & Cum_impacts1$dp_impact=="moderate"] <- "m3"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m3" & Cum_impacts1$dp_impact=="low"] <- "m3"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m3" & Cum_impacts1$dp_impact=="very low"] <- "m3"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m2" & Cum_impacts1$dp_impact=="high"] <- "High"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m2" & Cum_impacts1$dp_impact=="moderate"] <- "m2"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m2" & Cum_impacts1$dp_impact=="low"] <- "m2"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m2" & Cum_impacts1$dp_impact=="very low"] <- "m2"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m1" & Cum_impacts1$dp_impact=="high"] <- "High"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m1" & Cum_impacts1$dp_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m1" & Cum_impacts1$dp_impact=="low"] <- "m1"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m1" & Cum_impacts1$dp_impact=="very low"] <- "m1"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Low" & Cum_impacts1$dp_impact=="high"] <- "High"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Low" & Cum_impacts1$dp_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Low" & Cum_impacts1$dp_impact=="low"] <- "Low"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Low" & Cum_impacts1$dp_impact=="very low"] <- "Low"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Very low" & Cum_impacts1$dp_impact=="high"] <- "High"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Very low" & Cum_impacts1$dp_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Very low" & Cum_impacts1$dp_impact=="low"] <- "Low"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Very low" & Cum_impacts1$dp_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="NA" & Cum_impacts1$dp_impact=="high"] <- "High"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="NA" & Cum_impacts1$dp_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="NA" & Cum_impacts1$dp_impact=="low"] <- "Low"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="NA" & Cum_impacts1$dp_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact4[is.na(Cum_impacts1$cumImpact3) & Cum_impacts1$dp_impact=="high"] <- "High"
Cum_impacts1$cumImpact4[is.na(Cum_impacts1$cumImpact3) & Cum_impacts1$dp_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact4[is.na(Cum_impacts1$cumImpact3) & Cum_impacts1$dp_impact=="low"] <- "Low"
Cum_impacts1$cumImpact4[is.na(Cum_impacts1$cumImpact3) & Cum_impacts1$dp_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="High" & is.na(Cum_impacts1$dp_impact)] <- "High"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m3" & is.na(Cum_impacts1$dp_impact)] <- "m3"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m2" & is.na(Cum_impacts1$dp_impact)] <- "m2"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="m1" & is.na(Cum_impacts1$dp_impact)] <- "m1"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Low" & is.na(Cum_impacts1$dp_impact)] <- "Low"
Cum_impacts1$cumImpact4[Cum_impacts1$cumImpact3=="Very low" & is.na(Cum_impacts1$dp_impact)] <- "Very low"

Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="High" & Cum_impacts1$sg_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="High" & Cum_impacts1$sg_impact=="moderate"] <- "High"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="High" & Cum_impacts1$sg_impact=="low"] <- "High"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="High" & Cum_impacts1$sg_impact=="very low"] <- "High"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m3" & Cum_impacts1$sg_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m3" & Cum_impacts1$sg_impact=="moderate"] <- "m3"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m3" & Cum_impacts1$sg_impact=="low"] <- "m3"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m3" & Cum_impacts1$sg_impact=="very low"] <- "m3"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m2" & Cum_impacts1$sg_impact=="high"] <- "High"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m2" & Cum_impacts1$sg_impact=="moderate"] <- "m2"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m2" & Cum_impacts1$sg_impact=="low"] <- "m2"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m2" & Cum_impacts1$sg_impact=="very low"] <- "m2"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m1" & Cum_impacts1$sg_impact=="high"] <- "High"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m1" & Cum_impacts1$sg_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m1" & Cum_impacts1$sg_impact=="low"] <- "m1"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m1" & Cum_impacts1$sg_impact=="very low"] <- "m1"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Low" & Cum_impacts1$sg_impact=="high"] <- "High"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Low" & Cum_impacts1$sg_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Low" & Cum_impacts1$sg_impact=="low"] <- "Low"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Low" & Cum_impacts1$sg_impact=="very low"] <- "Low"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Very low" & Cum_impacts1$sg_impact=="high"] <- "High"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Very low" & Cum_impacts1$sg_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Very low" & Cum_impacts1$sg_impact=="low"] <- "Low"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Very low" & Cum_impacts1$sg_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="NA" & Cum_impacts1$sg_impact=="high"] <- "High"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="NA" & Cum_impacts1$sg_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="NA" & Cum_impacts1$sg_impact=="low"] <- "Low"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="NA" & Cum_impacts1$sg_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact5[is.na(Cum_impacts1$cumImpact4) & Cum_impacts1$sg_impact=="high"] <- "High"
Cum_impacts1$cumImpact5[is.na(Cum_impacts1$cumImpact4) & Cum_impacts1$sg_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact5[is.na(Cum_impacts1$cumImpact4) & Cum_impacts1$sg_impact=="low"] <- "Low"
Cum_impacts1$cumImpact5[is.na(Cum_impacts1$cumImpact4) & Cum_impacts1$sg_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="High" & is.na(Cum_impacts1$sg_impact)] <- "High"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m3" & is.na(Cum_impacts1$sg_impact)] <- "m3"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m2" & is.na(Cum_impacts1$sg_impact)] <- "m2"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="m1" & is.na(Cum_impacts1$sg_impact)] <- "m1"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Low" & is.na(Cum_impacts1$sg_impact)] <- "Low"
Cum_impacts1$cumImpact5[Cum_impacts1$cumImpact4=="Very low" & is.na(Cum_impacts1$sg_impact)] <- "Very low"

Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="High" & Cum_impacts1$fm_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="High" & Cum_impacts1$fm_impact=="moderate"] <- "High"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="High" & Cum_impacts1$fm_impact=="low"] <- "High"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="High" & Cum_impacts1$fm_impact=="very low"] <- "High"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m3" & Cum_impacts1$fm_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m3" & Cum_impacts1$fm_impact=="moderate"] <- "m3"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m3" & Cum_impacts1$fm_impact=="low"] <- "m3"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m3" & Cum_impacts1$fm_impact=="very low"] <- "m3"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m2" & Cum_impacts1$fm_impact=="high"] <- "High"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m2" & Cum_impacts1$fm_impact=="moderate"] <- "m2"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m2" & Cum_impacts1$fm_impact=="low"] <- "m2"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m2" & Cum_impacts1$fm_impact=="very low"] <- "m2"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m1" & Cum_impacts1$fm_impact=="high"] <- "High"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m1" & Cum_impacts1$fm_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m1" & Cum_impacts1$fm_impact=="low"] <- "m1"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m1" & Cum_impacts1$fm_impact=="very low"] <- "m1"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Low" & Cum_impacts1$fm_impact=="high"] <- "High"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Low" & Cum_impacts1$fm_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Low" & Cum_impacts1$fm_impact=="low"] <- "Low"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Low" & Cum_impacts1$fm_impact=="very low"] <- "Low"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Very low" & Cum_impacts1$fm_impact=="high"] <- "High"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Very low" & Cum_impacts1$fm_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Very low" & Cum_impacts1$fm_impact=="low"] <- "Low"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Very low" & Cum_impacts1$fm_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="NA" & Cum_impacts1$fm_impact=="high"] <- "High"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="NA" & Cum_impacts1$fm_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="NA" & Cum_impacts1$fm_impact=="low"] <- "Low"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="NA" & Cum_impacts1$fm_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact6[is.na(Cum_impacts1$cumImpact5) & Cum_impacts1$fm_impact=="high"] <- "High"
Cum_impacts1$cumImpact6[is.na(Cum_impacts1$cumImpact5) & Cum_impacts1$fm_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact6[is.na(Cum_impacts1$cumImpact5) & Cum_impacts1$fm_impact=="low"] <- "Low"
Cum_impacts1$cumImpact6[is.na(Cum_impacts1$cumImpact5) & Cum_impacts1$fm_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="High" & is.na(Cum_impacts1$fm_impact)] <- "High"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m3" & is.na(Cum_impacts1$fm_impact)] <- "m3"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m2" & is.na(Cum_impacts1$fm_impact)] <- "m2"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="m1" & is.na(Cum_impacts1$fm_impact)] <- "m1"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Low" & is.na(Cum_impacts1$fm_impact)] <- "Low"
Cum_impacts1$cumImpact6[Cum_impacts1$cumImpact5=="Very low" & is.na(Cum_impacts1$fm_impact)] <- "Very low"

Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="High" & Cum_impacts1$po_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="High" & Cum_impacts1$po_impact=="moderate"] <- "High"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="High" & Cum_impacts1$po_impact=="low"] <- "High"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="High" & Cum_impacts1$po_impact=="very low"] <- "High"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m3" & Cum_impacts1$po_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m3" & Cum_impacts1$po_impact=="moderate"] <- "m3"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m3" & Cum_impacts1$po_impact=="low"] <- "m3"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m3" & Cum_impacts1$po_impact=="very low"] <- "m3"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m2" & Cum_impacts1$po_impact=="high"] <- "High"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m2" & Cum_impacts1$po_impact=="moderate"] <- "m2"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m2" & Cum_impacts1$po_impact=="low"] <- "m2"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m2" & Cum_impacts1$po_impact=="very low"] <- "m2"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m1" & Cum_impacts1$po_impact=="high"] <- "High"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m1" & Cum_impacts1$po_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m1" & Cum_impacts1$po_impact=="low"] <- "m1"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m1" & Cum_impacts1$po_impact=="very low"] <- "m1"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Low" & Cum_impacts1$po_impact=="high"] <- "High"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Low" & Cum_impacts1$po_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Low" & Cum_impacts1$po_impact=="low"] <- "Low"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Low" & Cum_impacts1$po_impact=="very low"] <- "Low"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Very low" & Cum_impacts1$po_impact=="high"] <- "High"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Very low" & Cum_impacts1$po_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Very low" & Cum_impacts1$po_impact=="low"] <- "Low"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Very low" & Cum_impacts1$po_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="NA" & Cum_impacts1$po_impact=="high"] <- "High"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="NA" & Cum_impacts1$po_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="NA" & Cum_impacts1$po_impact=="low"] <- "Low"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="NA" & Cum_impacts1$po_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact7[is.na(Cum_impacts1$cumImpact6) & Cum_impacts1$po_impact=="high"] <- "High"
Cum_impacts1$cumImpact7[is.na(Cum_impacts1$cumImpact6) & Cum_impacts1$po_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact7[is.na(Cum_impacts1$cumImpact6) & Cum_impacts1$po_impact=="low"] <- "Low"
Cum_impacts1$cumImpact7[is.na(Cum_impacts1$cumImpact6) & Cum_impacts1$po_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="High" & is.na(Cum_impacts1$po_impact)] <- "High"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m3" & is.na(Cum_impacts1$po_impact)] <- "m3"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m2" & is.na(Cum_impacts1$po_impact)] <- "m2"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="m1" & is.na(Cum_impacts1$po_impact)] <- "m1"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Low" & is.na(Cum_impacts1$po_impact)] <- "Low"
Cum_impacts1$cumImpact7[Cum_impacts1$cumImpact6=="Very low" & is.na(Cum_impacts1$po_impact)] <- "Very low"

Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="High" & Cum_impacts1$WC_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="High" & Cum_impacts1$WC_impact=="moderate"] <- "High"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="High" & Cum_impacts1$WC_impact=="low"] <- "High"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="High" & Cum_impacts1$WC_impact=="very low"] <- "High"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m3" & Cum_impacts1$WC_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m3" & Cum_impacts1$WC_impact=="moderate"] <- "m3"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m3" & Cum_impacts1$WC_impact=="low"] <- "m3"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m3" & Cum_impacts1$WC_impact=="very low"] <- "m3"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m2" & Cum_impacts1$WC_impact=="high"] <- "High"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m2" & Cum_impacts1$WC_impact=="moderate"] <- "m2"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m2" & Cum_impacts1$WC_impact=="low"] <- "m2"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m2" & Cum_impacts1$WC_impact=="very low"] <- "m2"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m1" & Cum_impacts1$WC_impact=="high"] <- "High"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m1" & Cum_impacts1$WC_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m1" & Cum_impacts1$WC_impact=="low"] <- "m1"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m1" & Cum_impacts1$WC_impact=="very low"] <- "m1"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Low" & Cum_impacts1$WC_impact=="high"] <- "High"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Low" & Cum_impacts1$WC_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Low" & Cum_impacts1$WC_impact=="low"] <- "Low"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Low" & Cum_impacts1$WC_impact=="very low"] <- "Low"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Very low" & Cum_impacts1$WC_impact=="high"] <- "High"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Very low" & Cum_impacts1$WC_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Very low" & Cum_impacts1$WC_impact=="low"] <- "Low"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Very low" & Cum_impacts1$WC_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="NA" & Cum_impacts1$WC_impact=="high"] <- "High"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="NA" & Cum_impacts1$WC_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="NA" & Cum_impacts1$WC_impact=="low"] <- "Low"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="NA" & Cum_impacts1$WC_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact8[is.na(Cum_impacts1$cumImpact7) & Cum_impacts1$WC_impact=="high"] <- "High"
Cum_impacts1$cumImpact8[is.na(Cum_impacts1$cumImpact7) & Cum_impacts1$WC_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact8[is.na(Cum_impacts1$cumImpact7) & Cum_impacts1$WC_impact=="low"] <- "Low"
Cum_impacts1$cumImpact8[is.na(Cum_impacts1$cumImpact7) & Cum_impacts1$WC_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="High" & is.na(Cum_impacts1$WC_impact)] <- "High"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m3" & is.na(Cum_impacts1$WC_impact)] <- "m3"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m2" & is.na(Cum_impacts1$WC_impact)] <- "m2"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="m1" & is.na(Cum_impacts1$WC_impact)] <- "m1"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Low" & is.na(Cum_impacts1$WC_impact)] <- "Low"
Cum_impacts1$cumImpact8[Cum_impacts1$cumImpact7=="Very low" & is.na(Cum_impacts1$WC_impact)] <- "Very low"

Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="High" & Cum_impacts1$wo_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="High" & Cum_impacts1$wo_impact=="moderate"] <- "High"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="High" & Cum_impacts1$wo_impact=="low"] <- "High"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="High" & Cum_impacts1$wo_impact=="very low"] <- "High"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m3" & Cum_impacts1$wo_impact=="high"] <- "Very high"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m3" & Cum_impacts1$wo_impact=="moderate"] <- "m3"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m3" & Cum_impacts1$wo_impact=="low"] <- "m3"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m3" & Cum_impacts1$wo_impact=="very low"] <- "m3"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m2" & Cum_impacts1$wo_impact=="high"] <- "High"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m2" & Cum_impacts1$wo_impact=="moderate"] <- "m2"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m2" & Cum_impacts1$wo_impact=="low"] <- "m2"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m2" & Cum_impacts1$wo_impact=="very low"] <- "m2"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m1" & Cum_impacts1$wo_impact=="high"] <- "High"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m1" & Cum_impacts1$wo_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m1" & Cum_impacts1$wo_impact=="low"] <- "m1"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m1" & Cum_impacts1$wo_impact=="very low"] <- "m1"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Low" & Cum_impacts1$wo_impact=="high"] <- "High"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Low" & Cum_impacts1$wo_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Low" & Cum_impacts1$wo_impact=="low"] <- "Low"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Low" & Cum_impacts1$wo_impact=="very low"] <- "Low"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Very low" & Cum_impacts1$wo_impact=="high"] <- "High"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Very low" & Cum_impacts1$wo_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Very low" & Cum_impacts1$wo_impact=="low"] <- "Low"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Very low" & Cum_impacts1$wo_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="NA" & Cum_impacts1$wo_impact=="high"] <- "High"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="NA" & Cum_impacts1$wo_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="NA" & Cum_impacts1$wo_impact=="low"] <- "Low"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="NA" & Cum_impacts1$wo_impact=="very low"] <- "Very low"
Cum_impacts1$cumImpact9[is.na(Cum_impacts1$cumImpact8) & Cum_impacts1$wo_impact=="high"] <- "High"
Cum_impacts1$cumImpact9[is.na(Cum_impacts1$cumImpact8) & Cum_impacts1$wo_impact=="moderate"] <- "m1"
Cum_impacts1$cumImpact9[is.na(Cum_impacts1$cumImpact8) & Cum_impacts1$wo_impact=="low"] <- "Low"
Cum_impacts1$cumImpact9[is.na(Cum_impacts1$cumImpact8) & Cum_impacts1$wo_impact=="very low"] <- "very low"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="High" & is.na(Cum_impacts1$wo_impact)] <- "High"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m3" & is.na(Cum_impacts1$wo_impact)] <- "m3"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m2" & is.na(Cum_impacts1$wo_impact)] <- "m2"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="m1" & is.na(Cum_impacts1$wo_impact)] <- "m1"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Low" & is.na(Cum_impacts1$wo_impact)] <- "Low"
Cum_impacts1$cumImpact9[Cum_impacts1$cumImpact8=="Very low" & is.na(Cum_impacts1$wo_impact)] <- "Very low"

Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="High" & Cum_impacts1$shippingCategory=="high"] <- "Very high"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="High" & Cum_impacts1$shippingCategory=="moderate"] <- "High"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="High" & Cum_impacts1$shippingCategory=="low"] <- "High"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="High" & Cum_impacts1$shippingCategory=="very low"] <- "High"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m3" & Cum_impacts1$shippingCategory=="high"] <- "Very high"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m3" & Cum_impacts1$shippingCategory=="moderate"] <- "m3"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m3" & Cum_impacts1$shippingCategory=="low"] <- "m3"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m3" & Cum_impacts1$shippingCategory=="very low"] <- "m3"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m2" & Cum_impacts1$shippingCategory=="high"] <- "High"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m2" & Cum_impacts1$shippingCategory=="moderate"] <- "m2"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m2" & Cum_impacts1$shippingCategory=="low"] <- "m2"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m2" & Cum_impacts1$shippingCategory=="very low"] <- "m2"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m1" & Cum_impacts1$shippingCategory=="high"] <- "High"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m1" & Cum_impacts1$shippingCategory=="moderate"] <- "m1"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m1" & Cum_impacts1$shippingCategory=="low"] <- "m1"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m1" & Cum_impacts1$shippingCategory=="very low"] <- "m1"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Low" & Cum_impacts1$shippingCategory=="high"] <- "High"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Low" & Cum_impacts1$shippingCategory=="moderate"] <- "m1"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Low" & Cum_impacts1$shippingCategory=="low"] <- "Low"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Low" & Cum_impacts1$shippingCategory=="very low"] <- "Low"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Very low" & Cum_impacts1$shippingCategory=="high"] <- "High"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Very low" & Cum_impacts1$shippingCategory=="moderate"] <- "m1"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Very low" & Cum_impacts1$shippingCategory=="low"] <- "Low"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Very low" & Cum_impacts1$shippingCategory=="very low"] <- "Very low"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="NA" & Cum_impacts1$shippingCategory=="high"] <- "High"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="NA" & Cum_impacts1$shippingCategory=="moderate"] <- "m1"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="NA" & Cum_impacts1$shippingCategory=="low"] <- "Low"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="NA" & Cum_impacts1$shippingCategory=="very low"] <- "Very low"
Cum_impacts1$cumImpact10[is.na(Cum_impacts1$cumImpact9) & Cum_impacts1$shippingCategory=="high"] <- "High"
Cum_impacts1$cumImpact10[is.na(Cum_impacts1$cumImpact9) & Cum_impacts1$shippingCategory=="moderate"] <- "m1"
Cum_impacts1$cumImpact10[is.na(Cum_impacts1$cumImpact9) & Cum_impacts1$shippingCategory=="low"] <- "Low"
Cum_impacts1$cumImpact10[is.na(Cum_impacts1$cumImpact9) & Cum_impacts1$shippingCategory=="very low"] <- "Very low"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="High" & is.na(Cum_impacts1$shippingCategory)] <- "High"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m3" & is.na(Cum_impacts1$shippingCategory)] <- "m3"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m2" & is.na(Cum_impacts1$shippingCategory)] <- "m2"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="m1" & is.na(Cum_impacts1$shippingCategory)] <- "m1"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Low" & is.na(Cum_impacts1$shippingCategory)] <- "Low"
Cum_impacts1$cumImpact10[Cum_impacts1$cumImpact9=="Very low" & is.na(Cum_impacts1$shippingCategory)] <- "Very low"

Cum_impacts1$cell_area_km2 <- st_area(Cum_impacts1)/1000000
st_write(Cum_impacts1,"Q:/dfad/users/joeg/home/VMS/200921_HELCOM_Cumi/Results/CumI_CumulativeImpacts3.shp", delete_dsn = T)

#####################
#Loss
loss_total <- st_read("Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Data_layers\\Loss\\Loss_total_DK.shp")
loss_total_crs <- st_transform(loss_total,crs(Cum_impacts1))
loss_total_valid <- st_make_valid(loss_total_crs)
loss_total_valid$loss <- "Loss"
Cum_impacts1_valid <- st_make_valid(Cum_impacts1)

Cum_impacts2 <- st_join(Cum_impacts1_valid,loss_total_valid, left=TRUE, largest=TRUE)


Cum_impacts2$cumImpact_loss <-Cum_impacts2$cumImpact10 
Cum_impacts2$cumImpact_loss[Cum_impacts2$loss=="Loss"] <-"Loss" 

Cum_impacts2$cell_area_km2 <- st_area(Cum_impacts2)/1000000

st_write(Cum_impacts2,"Q:/dfad/users/joeg/home/VMS/200921_HELCOM_Cumi/Results/CumI_CumulativeImpacts_Loss2.shp", delete_dsn = T)

Cum_impacts2_loss <- Cum_impacts2

Cum_impacts2_loss$area_km2 <- st_area(Cum_impacts2_loss)/1000000

cumulative_impact_loss <- Cum_impacts2_loss %>%
  group_by(level_2, habitatComb, cumImpact_loss) %>%
  summarise(area_km2 = sum(area_km2))
st_geometry(cumulative_impact_loss) <- NULL

write.csv(cumulative_impact_loss, "Q:\\dfad\\users\\joeg\\home\\VMS\\200921_HELCOM_Cumi\\Grafer og tabeller\\cumulative_impact_loss.csv")


