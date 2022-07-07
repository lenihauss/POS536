#### R Script POS536 zooplankton distribution plot
# Author: H.Hauss (hhauss@geomar.de)
# Calculate zooplankton abundance and biomass from Ecotaxa tsv export
#POS536. Mesozooplankton was collected using a 
#Bongo net (200µm mesh size?), oblique hauls 
#samples were scanned using an Epson V750 pro scanner and imported into Ecotaxa

#load libraries (install first)
library(dplyr)
library(tidyr)
library(ggplot2)
library(maps)
library(mapdata)
library(scatterpie)

## set working directory -- change this to where you saved your .tsv data
setwd("V:/Daten/Cruises/POS536_SO")  

## list all ecotaxa in a folder
file_list <- list.files(pattern="*.tsv") # create list of all .tsv files in that folder
print(file_list)
##read in all ecotaxa tsvs into one dataframe and select variables to keep
data_raw <- do.call(rbind,lapply(file_list,read.csv, header=TRUE, quote = '"', sep = "\t"))[ ,c('sample_id','object_annotation_category','object_annotation_hierarchy', 'object_date', 'object_time','object_lat', 'object_lon','object_depth_min', 'object_depth_max', 'object_area', 'sample_volconc')]

## remove detritus and other non-living organisms from dataframe
data <-data_raw[-grep("not-living", data_raw$object_annotation_hierarchy),] 
rm(data_raw)

## check format of variables (numeric, character, factor, etc.)
str(data)

## rename some variables for easier/shorter names
names(data) = gsub(pattern = "object_", replacement = "", x = names(data))
data$spec_id <- data$annotation_category
##date format
data$date <- as.Date(as.character(data$date), format =  "%Y%m%d")
##time format
data$time <- sprintf("%06d", data$time)
data$time <- format(strptime(data$time, format="%H%M%S"), format = "%H:%M")
data$sample_volconc <- data$sample_volconc/1000 # to get to m3
metadata<- data %>% 
  distinct(sample_id, lat, lon, date, time, depth_min, depth_max,sample_volconc, .keep_all = F)
#seems lat is missing for 2019_pos536_bongo12l_lrg and volumes seem super large...check!!!

##calculate image area from pixels (pixel size = 10.6µm)
data$area_mm2 <- data$area*0.00011236

##convert area to biomass for different taxa -- these relationships are from Lehette & Hernandez-Leon 2006 and
data$biomass_ug <- with(data, 
                        ifelse(grepl("Calanoida", annotation_hierarchy),45.25*data$area_mm2^1.59,
                               ifelse(grepl("Annelida", annotation_hierarchy),43.38*data $area_mm2^1.54,
                                      ifelse(grepl("Chaetognatha", annotation_hierarchy),23.45*data $area_mm2^1.19,
                                             ifelse(grepl("Gastropoda", annotation_hierarchy),43.38*data $area_mm2^1.54,
                                                    44.78*data$area_mm2^1.56)))))
#from here, it is aggregation
#total abundance
samples_abund <- aggregate(biomass_ug ~ (sample_id+sample_volconc), data, length)
samples_abund$total_abundance_m3 <- samples_abund$biomass_ug/samples_abund$sample_volconc
samples_abund <-samples_abund[c("sample_id","total_abundance_m3")]
#same for individual taxonomic categories
species_abund <- aggregate(biomass_ug ~ (spec_id+sample_id + sample_volconc), data, length)
species_abund$abundance_m3 <- species_abund$biomass_ug/species_abund$sample_volconc
species_abund <-species_abund[c("sample_id","spec_id","abundance_m3")]
#spread to wide format
species_abund <-species_abund %>% spread(spec_id, abundance_m3)
#merge to metadata and total abundance
abundance <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "sample_id", all.x = TRUE),
       list(metadata, samples_abund, species_abund))
#replace NAs with zeros
abundance[is.na(abundance)] <- 0

#total biomass
samples_biomass <- aggregate(biomass_ug ~ (sample_id+sample_volconc), data, sum)
samples_biomass$total_biomass_ug_m3 <- samples_biomass$biomass_ug/samples_biomass$sample_volconc
samples_biomass <-samples_biomass[c("sample_id","total_biomass_ug_m3")]
#same for individual taxonomic categories
species_biomass <- aggregate(biomass_ug ~ (spec_id+sample_id + sample_volconc), data, sum)
species_biomass$biomass_ug_m3 <- species_biomass$biomass_ug/species_biomass$sample_volconc
species_biomass <-species_biomass[c("sample_id","spec_id","biomass_ug_m3")]
#spread to wide format
species_biomass <-species_biomass %>% spread(spec_id, biomass_ug_m3)
#merge to metadata and total abundance
biomass <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "sample_id", all.x = TRUE),
                    list(metadata, samples_biomass, species_biomass))
#replace NAs with zeros
biomass[is.na(biomass)] <- 0
#remove intermediate dataframes
rm(data, samples_biomass, species_biomass, species_abund, samples_abund)
##export summary data in case this shall be published (e.g. on Pangaea)
write.table(biomass, file = "POS536_biomass.txt")
write.table(abundance, file = "POS536_abundance.txt")

##calculate fraction from total abundance for calanoids, chaetognaths, krill and others
abundance$calanoid_fraction = abundance$Calanoida/abundance$total_abundance_m3
abundance$harpacticoid_fraction = abundance$Harpacticoida/abundance$total_abundance_m3
abundance$salp_fraction = abundance$Thaliacea/abundance$total_abundance_m3
abundance$fish_fraction = abundance$Actinopterygii/abundance$total_abundance_m3
abundance$krill_fraction = abundance$Euphausiacea/abundance$total_abundance_m3
abundance$chaetognath_fraction = abundance$Chaetognatha/abundance$total_abundance_m3
abundance$other_fraction = 1 - abundance$calanoid_fraction - abundance$chaetognath_fraction - abundance$krill_fraction - abundance$fish_fraction - abundance$salp_fraction  - abundance$harpacticoid_fraction
fractions <-abundance[c("sample_id","lat", "lon", "total_abundance_m3", "calanoid_fraction", "krill_fraction","chaetognath_fraction", "fish_fraction", "salp_fraction", "harpacticoid_fraction", "other_fraction")]
write.table(fractions, file = "POS536_abundance_fractions.txt")


##calculate fraction from total biomass for calanoids, chaetognaths, krill and others
biomass$calanoid_fraction = biomass$Calanoida/biomass$total_biomass_ug_m3
biomass$krill_fraction = biomass$Euphausiacea/biomass$total_biomass_ug_m3
biomass$chaetognath_fraction = biomass$Chaetognatha/biomass$total_biomass_ug_m3
biomass$other_fraction = 1 - biomass$calanoid_fraction - biomass$chaetognath_fraction - biomass$krill_fraction
fractions <-biomass[c("sample_id","lat", "lon", "total_biomass_ug_m3", "calanoid_fraction", "krill_fraction","chaetognath_fraction", "other_fraction")]
write.table(fractions, file = "POS536_fractions.txt")

################plot on map
###load coastline data for Europe, using the map_data function
Area <- map_data("world")


#plot map of total biomass as bubble plot
ggplot() +
  theme_bw()+
  geom_polygon(fill="darkgoldenrod1", color = "black", data = Area, aes(x=long, y = lat, group = group)) + ##this adds the coastlines
  geom_point(data = biomass, aes(x = lon, y = lat, size = total_biomass_ug_m3), color = "red", alpha=0.7) +
  coord_quickmap(xlim=c(-40,-4), ylim=c(25,50))  + ##map boundaries
  labs(title = "Total Biomass",x = "Longitude °W", y = "Latitude °N", size = "ug / m3")                                         ##axis labels

ggsave("V:/Daten/Cruises/POS536_SO/total_biomass_bubble.png", width = 7, height = 5, bg = "transparent")

#plot map with piecharts to show fraction of different groups
#something is still wrong here, check...
ggplot() +
  theme_bw()+
  geom_polygon(fill="darkgoldenrod1", color = "black", data = Area, aes(x=long, y = lat, group = group)) +  ## coastlines                                                       ##map boundaries
  labs(title = "POS536 & SO579",x = "Longitude °W", y = "Latitude °N", size = "Biomass (mg / m²)")  +                          ##axis labels
  #geom_point(data = fractions, aes(x = lon, y = lat, size = total_biomass_ug_m3), pch=21, fill = "transparent", color= "black") +
  geom_scatterpie(aes(x=lon, y=lat, group = sample_id, r = 3), 
                  data = fractions, cols = colnames(fractions[,c(5:11)])) +
  #scale_fill_maname="Group",labels=c("Calanoida","Krill","Chaetognatha", "Fish", "Salp", "Harpacticoid", "other"), values=c("green1","coral1","darkslateblue", "yellow")) +
  coord_fixed(xlim=c(-40,-4), ylim=c(25,50))

ggsave("V:/Daten/Cruises/POS536_SO/biomass_pie.png", width = 7, height = 5, bg = "transparent")