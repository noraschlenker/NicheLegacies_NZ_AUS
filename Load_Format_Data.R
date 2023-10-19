#Upload GBIF data, filter to only wanted points, combine with climate data, and create one list for each family with desired data for analysis
#October 2023

# *** NOTE: This code filters species occurrence data and creates the file Q1_data. 
# Q1_data is included in the data folder, you do not need to run this code to run the analysis code ****

library(tidyverse)
library(raster)
library(sp)
library(ggplot2)
library(maps)
library(rgdal)
library(gridExtra)
library(grid)

#Create file paths and load data ------------
#dir.create(path = "data")
#dir.create(path = "output")
load("data/tablespec.Rdata")
bioclim.data <- getData(name = "worldclim", var = "bio", res = 2.5, path = "data/")
load("./data/clim_data.RData")
bioclim.data <- clim_data
countries <- raster("./data/CountryRaster_40km.tif")

#Create Load Family Function ----------

#This function loads each family, filters to only wanted points, separates out individual genera, determines if points are in New Zealand or Australia combines data into a list for further analysis
Load_Family <- function(fam, gen, clim_data, country_raster, table_spec, LatMin = -50, LatMax = 0, LongMin = 110, LongMax = 180) {
  if(length(fam) > 1 | class(fam) != "character") {
    stop("Family must be a character of length 1", call. = FALSE)
  }
  if(length(gen) > 6 | length(gen) < 2 | class(gen) != "character") {
    stop("Genera must be a character vector of length 2, 3, 4, 5, or 6", call. = FALSE)
  }
  
  #get plant occurrence data from hard drive and name it "fam_unfiltered", row.names = NULL because original data does not have row numbers
  fam_unfiltered <- readr::read_tsv(file = paste0("./data/", fam,"/Data.csv"), col_types = table_spec) 
  
  convert_to_raster <- function(x, gnus){
    #filters out genera points
    pnts <- x %>%
      dplyr::filter(genus == gnus) %>%
      dplyr::select(decimalLatitude, decimalLongitude) %>%
      dplyr::mutate(presence = 1)
    #converts points to spatial points
    coordinates(pnts) <- ~decimalLongitude + decimalLatitude
    #rasterizes points
    rstr <- rasterize(pnts, clim_data, 'presence')
    names(rstr)<- gnus
    rstr
  }
  
  #Filter and organize plant occurrence data into desired long format
  #creates a table from fam_unfiltered table which is filtered by Lat/Long, Human Observation, and Genus
  #then removes observations that are not identified to species level and which have a coordinate uncertainty of more than 2250 meters
  
  #filters fam_1 table to include only the area, observation, and genera in study 
  
  if(length(gen) == 2){
    obs_list <- fam_unfiltered %>%
      dplyr::select(gbifID, family, genus, species,decimalLatitude, decimalLongitude, basisOfRecord,
                    countryCode,coordinateUncertaintyInMeters)%>%
      dplyr::filter(decimalLatitude <= LatMax, decimalLatitude >= LatMin, decimalLongitude >= LongMin, decimalLongitude <= LongMax, 
                    basisOfRecord == "HUMAN_OBSERVATION" | basisOfRecord == "PRESERVED_SPECIMEN",  
                    genus == paste(gen[1])| genus == paste(gen[2]), species != "") %>%
      dplyr::filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 2250) %>%
      dplyr::arrange(genus, species)
    rm(fam_unfiltered)
    g1 <- convert_to_raster(obs_list, gnus = paste(gen[1]))  
    g2 <- convert_to_raster(obs_list, gnus = paste(gen[2]))
    final_raster <- stack(clim_data, g1, g2, country_raster)
    rm(g1, g2, clim_data, country_raster, table_spec)
    
  } else {
    if(length(gen) == 3){
      obs_list <- fam_unfiltered %>%
        dplyr::select(gbifID, family, genus, species,decimalLatitude, decimalLongitude, basisOfRecord,
                      countryCode,coordinateUncertaintyInMeters)%>%
        dplyr::filter(decimalLatitude <= LatMax, decimalLatitude >= LatMin, decimalLongitude >= LongMin, decimalLongitude <= LongMax, 
                      basisOfRecord == "HUMAN_OBSERVATION" | basisOfRecord == "PRESERVED_SPECIMEN",  
                      genus == paste(gen[1])| genus == paste(gen[2])| genus == paste(gen[3]), species != "") %>%
        dplyr::filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 2250) %>%
        dplyr::arrange(genus, species)
      rm(fam_unfiltered)
      g1 <- convert_to_raster(obs_list, gnus = paste(gen[1]))  
      g2 <- convert_to_raster(obs_list, gnus = paste(gen[2]))
      g3 <- convert_to_raster(obs_list, gnus = paste(gen[3]))
      final_raster <- stack(clim_data, g1, g2, g3, country_raster)
      rm(g1, g2, g3, clim_data, country_raster, table_spec)
    } else {
      if(length(gen) == 4){
        obs_list <- fam_unfiltered %>%
          dplyr::select(gbifID, family, genus, species,decimalLatitude, decimalLongitude, basisOfRecord,
                        countryCode,coordinateUncertaintyInMeters)%>%
          dplyr::filter(decimalLatitude <= LatMax, decimalLatitude >= LatMin, decimalLongitude >= LongMin, decimalLongitude <= LongMax, 
                        basisOfRecord == "HUMAN_OBSERVATION" | basisOfRecord == "PRESERVED_SPECIMEN",  
                        genus == paste(gen[1])| genus == paste(gen[2])| genus == paste(gen[3]) | genus == paste(gen[4]), species != "") %>%
          dplyr::filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 2250) %>%
          dplyr::arrange(genus, species)
        rm(fam_unfiltered)
        g1 <- convert_to_raster(obs_list, gnus = paste(gen[1]))  
        g2 <- convert_to_raster(obs_list, gnus = paste(gen[2]))
        g3 <- convert_to_raster(obs_list, gnus = paste(gen[3]))
        g4 <- convert_to_raster(obs_list, gnus = paste(gen[4]))
        final_raster <- stack(clim_data, g1, g2, g3, g4, country_raster)
        rm(g1, g2, g3, g4, clim_data, country_raster, table_spec)
      } else {
        if(length(gen) == 5){
          obs_list <- fam_unfiltered %>%
            dplyr::select(gbifID, family, genus, species,decimalLatitude, decimalLongitude, basisOfRecord,
                          countryCode,coordinateUncertaintyInMeters)%>%
            dplyr::filter(decimalLatitude <= LatMax, decimalLatitude >= LatMin, decimalLongitude >= LongMin, decimalLongitude <= LongMax, 
                          basisOfRecord == "HUMAN_OBSERVATION" | basisOfRecord == "PRESERVED_SPECIMEN",  
                          genus == paste(gen[1])| genus == paste(gen[2])| genus == paste(gen[3]) | genus == paste(gen[4])| genus == paste(gen[5]), species != "") %>%
            dplyr::filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 2250) %>%
            dplyr::arrange(genus, species)
          rm(fam_unfiltered)
          g1 <- convert_to_raster(obs_list, gnus = paste(gen[1]))  
          g2 <- convert_to_raster(obs_list, gnus = paste(gen[2]))
          g3 <- convert_to_raster(obs_list, gnus = paste(gen[3]))
          g4 <- convert_to_raster(obs_list, gnus = paste(gen[4]))
          g5 <- convert_to_raster(obs_list, gnus = paste(gen[5]))
          final_raster <- stack(clim_data, g1, g2, g3, g4, g5, country_raster)
          rm(g1, g2, g3, g4, g5, clim_data, country_raster, table_spec)
        } else {
          if(length(gen) == 6){
            obs_list <- fam_unfiltered %>%
              dplyr::select(gbifID, family, genus, species,decimalLatitude, decimalLongitude, basisOfRecord,
                            countryCode,coordinateUncertaintyInMeters)%>%
              dplyr::filter(decimalLatitude <= LatMax, decimalLatitude >= LatMin, decimalLongitude >= LongMin, decimalLongitude <= LongMax, 
                            basisOfRecord == "HUMAN_OBSERVATION" | basisOfRecord == "PRESERVED_SPECIMEN",  
                            genus == paste(gen[1])| genus == paste(gen[2])| genus == paste(gen[3]) | genus == paste(gen[4])| genus == paste(gen[5]) | genus == paste(gen[6]), species != "") %>%
              dplyr::filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 2250) %>%
              dplyr::arrange(genus, species)
            rm(fam_unfiltered)
            g1 <- convert_to_raster(obs_list, gnus = paste(gen[1]))  
            g2 <- convert_to_raster(obs_list, gnus = paste(gen[2]))
            g3 <- convert_to_raster(obs_list, gnus = paste(gen[3]))
            g4 <- convert_to_raster(obs_list, gnus = paste(gen[4]))
            g5 <- convert_to_raster(obs_list, gnus = paste(gen[5]))
            g6 <- convert_to_raster(obs_list, gnus = paste(gen[6]))
            final_raster <- stack(clim_data, g1, g2, g3, g4, g5, g6, country_raster)
            rm(g1, g2, g3, g4, g5, g6, clim_data, country_raster, table_spec)
            
          } else{
            stop("Incorrect number of genera entered", call. = FALSE)
          }
        }
      }
    }
  }
  
  #Convert final_raster to data.frame
  final_df_long <- as.data.frame(final_raster, row.names=NULL, optional=FALSE, xy=TRUE, na.rm=FALSE, long=FALSE)
  colnames(final_df_long) <- c(colnames(final_df_long)[1:length(colnames(final_df_long))-1], "Country") #rename final column
  rm(final_raster) 
  
  #Filter final_df_long to include only the study area
  final_df_short <- dplyr::filter(final_df_long, y <= LatMax, y >= LatMin, x >= LongMin, x <= LongMax)
  
  #Remove row with no climate variables (i.e. Oceans)
  final_df <- final_df_short[complete.cases(final_df_short[ , 3:21]),]
  
  Country_Abreviations <- data.frame(code = c(26, 93, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221),
                                     country_abrv = c("FJI",  "NZL", 	"IDN", 	"TLS",  "AUS", "NRU",  "NCL", "NFK", 	"PNG", 	"SLB",  "TUV",  "VUT"))
  final_df <- left_join(final_df, Country_Abreviations, by = c("Country" = "code"))
  
  #Put data into list and save 
  list(family = fam, genera = gen, obs = obs_list, points = final_df)
}


Summary_Tables <- function(x, fam){
  s1 <- x %>%
    dplyr::group_by(genus) %>%
    dplyr::summarise(record_count = n(),
                     species_count = n_distinct(species)) %>%
    dplyr::arrange(desc(record_count))
  s2<-  rbind(s1, data.frame(genus ="Total", record_count=sum(s1$record_count), species_count = sum(s1$species_count)))
  jpeg(paste0("output/", fam, "GeneraSummary.jpeg"))
  grid.table(s2)
  dev.off()  
  s3 <- x %>%
    dplyr::group_by(species, genus, countryCode) %>%
    dplyr::summarise(record_count = n()) %>%
    dplyr::arrange(genus, species, countryCode, desc(record_count))
  write.csv(s3, paste0("output/", fam, "SpeciesSummary.csv"))
  
  p1 <- map_data("world")
  p2 <- ggplot() + geom_polygon(data = p1, aes(x = long, y = lat, group = group), fill = "grey70", col =
                                  "black") + xlim(110, 180) + ylim(-50, 0)
                                p3 <- p2 +
                                  geom_point(data = x, aes(x=decimalLongitude, y=decimalLatitude, col = genus)) +
                                  labs(y= "Latitude", x = "Longitude", colour = "Genus") +
                                  ggtitle(paste(fam, "Occurrence Map"))
                                ggsave(paste0("output/", fam, "OccurrenceMap.jpeg"),p3, width = 6, height = 4, units =
                                         "in")
                                rm(s1, s2, s3, p1, p2, p3)
                                
}

# Create family data ------------
# First download GBIF data using the references in Appendix 1
# NOTE: GBIF data must be formatted in a folder with the folder name equal to the family and the spreadsheet named "Data"

Family <- "Araucariaceae"
Genera <- c("Araucaria", "Agathis")
Araucariaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Araucariaceae[["obs"]], fam = Family)
save(Araucariaceae, file = "data_PS/Araucariaceae.Rdata")
rm(Araucariaceae)


Family <- "Arecaceae"
Genera <- c("Calamus", "Cocos", "Nypa", "Rhopalostylis")
Arecaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Arecaceae[["obs"]], fam = Family)
save(Arecaceae, file = "data_PS/Arecaceae.Rdata")
rm(Arecaceae)

Family <- "Argophyllaceae"
Genera <- c("Argophyllum", "Corokia")
Argophyllaceae  <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Argophyllaceae[["obs"]], fam = Family)
save(Argophyllaceae, file = "data_PS/Argophyllaceae.Rdata")
rm(Argophyllaceae)

Family <- "Elaeocarpaceae"
Genera <- c("Sloanea", "Elaeocarpus", "Aristotelia")
Elaeocarpaceae  <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Elaeocarpaceae[["obs"]], fam = Family)
save(Elaeocarpaceae, file = "data_PS/Elaeocarpaceae.Rdata")
rm(Elaeocarpaceae)

Family <- "Euphorbiaceae"
Genera <- c("Mallotus", "Euphorbia")
Euphorbiaceae   <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Euphorbiaceae[["obs"]], fam = Family)
save(Euphorbiaceae , file = "data_PS/Euphorbiaceae.Rdata")
rm(Euphorbiaceae)

Family <- "Fabaceae"
Genera <- c("Acacia", "Caesalpinia", "Sophora", "Carmichaelia", "Clianthus", "Canavalia")
Fabaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Fabaceae[["obs"]], fam = Family)
save(Fabaceae, file = "data_PS/Fabaceae.Rdata")
rm(Fabaceae)

Family <- "Lauraceae"
Genera <- c("Cryptocarya", "Beilschmiedia", "Litsea")
Lauraceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Lauraceae[["obs"]], fam = Family)
save(Lauraceae, file = "data_PS/Lauraceae.Rdata")
rm(Lauraceae)

Family <- "Malvaceae"
Genera <- c("Brachychiton", "Bombax", "Ungeria", "Hoheria", "Plagianthus", "Entelea")
Malvaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Malvaceae[["obs"]], fam = Family)
save(Malvaceae, file = "data_PS/Malvaceae.Rdata")
rm(Malvaceae)

Family <- "Onagraceae"
Genera <- c("Ludwigia", "Fuchsia", "Epilobium")
Onagraceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Onagraceae[["obs"]], fam = Family)
save(Onagraceae, file = "data_PS/Onagraceae.Rdata")
rm(Onagraceae)

Family <- "Paracryphiaceae"
Genera <- c("Paracryphia", "Quintinia")
Paracryphiaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Paracryphiaceae[["obs"]], fam = Family)
save(Paracryphiaceae, file = "data_PS/Paracryphiaceae.Rdata")
rm(Paracryphiaceae)

Family <- "Podocarpaceae"
Genera <- c("Microcachrys", "Podocarpus", "Dacrydium")
Podocarpaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Podocarpaceae[["obs"]], fam = Family)
save(Podocarpaceae, file = "data_PS/Podocarpaceae.Rdata")
rm(Podocarpaceae)

Family <- "Proteaceae"
Genera <- c("Banksia", "Beauprea", "Alloxylon", "Stenocarpus", "Toronia", "Knightia")
Proteaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Proteaceae[["obs"]], fam = Family)
save(Proteaceae, file = "data_PS/Proteaceae.Rdata")
rm(Proteaceae)

Family <- "Santalaceae"
Genera <- c("Amphorogyne", "Exocarpos", "Korthalsella")
Santalaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Santalaceae[["obs"]], fam = Family)
save(Santalaceae, file = "data_PS/Santalaceae.Rdata")
rm(Santalaceae)

Family <- "Sapindaceae"
Genera <- c("Cupaniopsis", "Mischocarpus", "Alectryon")
Sapindaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Sapindaceae[["obs"]], fam = Family)
save(Sapindaceae, file = "data_PS/Sapindaceae.Rdata")
rm(Sapindaceae)

Family <- "Strasburgeriaceae"
Genera <- c("Strasburgeria", "Ixerba")
Strasburgeriaceae <- Load_Family(fam = Family, gen = Genera, clim_data = bioclim.data, country_raster = countries, table_spec = table_specs)
Summary_Tables(Strasburgeriaceae[["obs"]], fam = Family)
save(Strasburgeriaceae, file = "data_PS/Strasburgeriaceae.Rdata")
rm(Strasburgeriaceae)


Points_long <- cbind(Araucariaceae[[4]], 
                     Arecaceae[[4]][,22:(21+length(Arecaceae[[4]])-22)],
                     Argophyllaceae[[4]][,22:(21+length(Argophyllaceae[[4]])-22)],
                     Elaeocarpaceae[[4]][,22:(21+length(Elaeocarpaceae[[4]])-22)],
                     Euphorbiaceae[[4]][,22:(21+length(Euphorbiaceae[[4]])-22)],
                     Fabaceae[[4]][,22:(21+length(Fabaceae[[4]])-22)],
                     Lauraceae[[4]][,22:(21+length(Lauraceae[[4]])-22)],
                     Malvaceae[[4]][,22:(21+length(Malvaceae[[4]])-22)],
                     Onagraceae[[4]][,22:(21+length(Onagraceae[[4]])-22)],
                     Paracryphiaceae[[4]][,22:(21+length(Paracryphiaceae[[4]])-22)],
                     Podocarpaceae[[4]][,22:(21+length(Podocarpaceae[[4]])-22)],
                     Proteaceae[[4]][,22:(21+length(Proteaceae[[4]])-22)],
                     Santalaceae[[4]][,22:(21+length(Santalaceae[[4]])-22)],
                     Sapindaceae[[4]][,22:(21+length(Sapindaceae[[4]])-22)],
                     Strasburgeriaceae[[4]][,22:(21+length(Strasburgeriaceae[[4]])-22)]
                     )

# Create filtered master table for question1 named "Q1_data" ----------
Australia <- Points_long %>%
  dplyr::filter(Country == 214) %>%
  dplyr::select(x:bio19)

Araucaria <- Points_long %>%
  dplyr::filter(Araucaria == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Agathis <- Points_long %>%
  dplyr::filter(Agathis == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Argophyllum <- Points_long %>%
  dplyr::filter(Argophyllum == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Corokia<- Points_long %>%
  dplyr::filter(Corokia == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Sloanea<- Points_long %>%
  dplyr::filter(Sloanea == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Elaeocarpus <- Points_long %>%
  dplyr::filter(Elaeocarpus == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Aristotelia <- Points_long %>%
  dplyr::filter(Aristotelia == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Mallotus <- Points_long %>%
  dplyr::filter(Mallotus == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Euphorbia <- Points_long %>%
  dplyr::filter(Euphorbia == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Acacia <- Points_long %>%
  dplyr::filter(Acacia == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Caesalpinia <- Points_long %>%
  dplyr::filter(Caesalpinia == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Sophora <- Points_long %>%
  dplyr::filter(Sophora == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Cryptocarya <- Points_long %>%
  dplyr::filter(Cryptocarya == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Beilschmiedia <- Points_long %>%
  dplyr::filter(Beilschmiedia == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Litsea <- Points_long %>%
  dplyr::filter(Litsea == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Ludwigia <- Points_long %>%
  dplyr::filter(Ludwigia == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Epilobium <- Points_long %>%
  dplyr::filter(Epilobium == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Microcachrys <- Points_long %>%
  dplyr::filter(Microcachrys == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Podocarpus <- Points_long %>%
  dplyr::filter(Podocarpus == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Cupaniopsis <- Points_long %>%
  dplyr::filter(Cupaniopsis == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Mischocarpus <- Points_long %>%
  dplyr::filter(Mischocarpus == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Alectryon <- Points_long %>%
  dplyr::filter(Alectryon == 1 & Country == 214) %>%
  dplyr::select(x:bio19)

Q1_data <- list("Australia" = Australia, "Araucaria" = Araucaria, "Agathis" = Agathis, "Argophyllum" = Argophyllum, 
                "Corokia" = Corokia, "Sloanea" = Sloanea, "Elaeocarpus"= Elaeocarpus, "Aristotelia" = Aristotelia,
                "Mallotus" = Mallotus, "Euphorbia" = Euphorbia, "Acacia" = Acacia, "Caesalpinia" = Caesalpinia,
                "Sophora" = Sophora, "Cryptocarya" = Cryptocarya, "Beilschmiedia" = Beilschmiedia, "Litsea" = Litsea,
                "Ludwigia" = Ludwigia, "Epilobium" = Epilobium, "Microcachrys" = Microcachrys, "Podocarpus" = Podocarpus,
                "Cupaniopsis" = Cupaniopsis, "Mischocarpus" = Mischocarpus, "Alectryon" = Alectryon)
save(Q1_data, file = "output/Q1_data.Rdata")

#Analysis continued in "Niche_Comparison_Analysis.R" file
