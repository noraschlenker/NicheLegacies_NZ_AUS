# Using data from "Load_Format_Data.R" final data Q1_data - complete niche analysis 
# This code completes the niche analysis and creates figures

library(ecospat)
library(tidyverse)
library(ade4)
library(ggrepel)
library(ggtext)
library(tidytext)
library(overlapping)

# Loadt Data -----------------
# Create climate subsets and clips for eco-regions --------
#load polygond for outlines of tropical and temperate Australia, full Australia, and New Zealand
proj_data <- sp::CRS(proj4string(clim_data))
clip_poly <- readOGR("./data/shpfiles/Beck_KG_0083_poly_TropTemp_1km_shp.shp")
trop_poly <- readOGR("./data/shpfiles/Beck_KG_Trop_1m_shp.shp")
temp_poly <- readOGR("./data/NZ_AUS_NicheComparison/Beck_KG_Temp_1m_shp.shp")
aus_poly <- readOGR("./data/shpfiles/AUS_shp.shp")
nz_poly <- readOGR("./data/shpfiles/NZ_outline_shp.shp")

#Clip each respective region with the climate data to make climate data tables for each region
clim_clip <- raster::extract(clim_data, clip_poly)
clim_clip_tbl <- clim_clip[[1]]
clim_clip_tbl <- clim_clip_tbl[complete.cases(clim_clip_tbl),]
class(clim_clip_tbl) <- "numeric"

clim_trop <- raster::extract(clim_data, trop_poly)
clim_trop_tbl <- clim_trop[[1]]
clim_trop_tbl <- clim_trop_tbl[complete.cases(clim_trop_tbl),]
class(clim_trop_tbl) <- "numeric"

clim_temp <- raster::extract(clim_data, temp_poly)
clim_temp_tbl <- clim_temp[[1]]
clim_temp_tbl <- clim_temp_tbl[complete.cases(clim_temp_tbl),]
class(clim_temp_tbl) <- "numeric"

clim_aus <- raster::extract(clim_data, aus_poly)
clim_aus_tbl <- clim_aus[[1]]
clim_aus_tbl <- clim_aus_tbl[complete.cases(clim_aus_tbl),]
class(clim_aus_tbl) <- "numeric"

clim_nz <- raster::extract(clim_data, nz_poly)
clim_nz_tbl <- clim_nz[[1]]
clim_nz_tbl <- clim_nz_tbl[complete.cases(clim_nz_tbl),]
class(clim_nz_tbl) <- "numeric"

nz_in_aus <- inner_join(round(data.frame(clim_aus_tbl),2), round(data.frame(clim_nz_tbl), 2))
aus <- data.frame(clim_aus_tbl)
nz_in_aus <- data.frame(clim_nz_tbl) %>% dplyr::filter(bio1 > min(aus$bio1) & bio1 < max(aus$bio1) &
                                                         bio2 > min(aus$bio2) & bio2 < max(aus$bio2) &
                                                         bio3 > min(aus$bio3) & bio3 < max(aus$bio3) &
                                                         bio4 > min(aus$bio4) & bio4 < max(aus$bio4) &
                                                         bio5 > min(aus$bio5) & bio5 < max(aus$bio5) &
                                                         bio6 > min(aus$bio6) & bio6 < max(aus$bio6) &
                                                         bio7 > min(aus$bio7) & bio7 < max(aus$bio7) &
                                                         bio8 > min(aus$bio8) & bio8 < max(aus$bio8) &
                                                         bio9 > min(aus$bio9) & bio9 < max(aus$bio9) &
                                                         bio10 > min(aus$bio10) & bio10 < max(aus$bio10) &
                                                         bio11 > min(aus$bio11) & bio11 < max(aus$bio11) &
                                                         bio12 > min(aus$bio12) & bio12 < max(aus$bio12) &
                                                         bio13 > min(aus$bio13) & bio13 < max(aus$bio13) &
                                                         bio14 > min(aus$bio14) & bio14 < max(aus$bio14) &
                                                         bio15 > min(aus$bio15) & bio15 < max(aus$bio15) &
                                                         bio16 > min(aus$bio16) & bio16 < max(aus$bio16) &
                                                         bio17 > min(aus$bio17) & bio17 < max(aus$bio17) &
                                                         bio18 > min(aus$bio18) & bio18 < max(aus$bio18) &
                                                         bio19 > min(aus$bio19) & bio19 < max(aus$bio19))

nz_in_aus <-as.matrix(nz_in_aus)



# Download genera data ------------
names(Q1_data)
#extract climate data from genera for only area in tropical and temperate Australia 
extract_genus <- function(genus, #character of genus (capitol first letter)
                          clim = clim_data, #climate raster stack
                          pr = proj_data, #projection data
                          cl = clip_poly #clip data
){
  x <- Q1_data[[genus]][,1:2]
  x_sites <- SpatialPoints(x, proj4string = pr)
  x_clim <- raster::extract(clim, x_sites[cl,])
  x_clim
}

Araucaria <- extract_genus("Araucaria")
Agathis <- extract_genus("Agathis")
Euphorbia <- extract_genus("Euphorbia")
Mallotus <- extract_genus("Mallotus")
Argophyllum <- extract_genus("Argophyllum")
Corokia <- extract_genus("Corokia")
Sloanea <- extract_genus("Sloanea")
Elaeocarpus <- extract_genus("Elaeocarpus")
Aristotelia <- extract_genus("Aristotelia")
Acacia <- extract_genus("Acacia")
Sophora <- extract_genus("Sophora")
Caesalpinia <- extract_genus("Caesalpinia")
Cryptocarya <- extract_genus("Cryptocarya")
Beilschmiedia <- extract_genus("Beilschmiedia")
Litsea <- extract_genus("Litsea")
Ludwigia <- extract_genus("Ludwigia")
Epilobium <- extract_genus("Epilobium")
Microcachrys <- extract_genus("Microcachrys")
Podocarpus <- extract_genus("Podocarpus")
Cupaniopsis <- extract_genus("Cupaniopsis")
Mischocarpus <- extract_genus("Mischocarpus")
Alectryon <- extract_genus("Alectryon")

# Create climate space PCA (Australia tropical/temperate) -----------
bioclims <- c(1:19)

pca_env <- dudi.pca(clim_clip_tbl[,bioclims], scannf = FALSE, nf = 2)

jpeg("./figures/pca_contrib.jpeg")
ecospat.plot.contrib(contrib=pca_env$co, eigen=pca_env$eig)
dev.off()

scores_globclim <- pca_env$li

#Get PCA score for different subsets of Australian climate space
scores_trop <- suprow(pca_env,clim_trop_tbl[,bioclims])$li
scores_temp <- suprow(pca_env,clim_temp_tbl[,bioclims])$li
scores_AUS <- suprow(pca_env,clim_aus_tbl[,bioclims])$li
#Get PCA score for NZ climate space
scores_NZ <- suprow(pca_env,clim_nz_tbl[,bioclims])$li

#Get PCA scores from AUS that are not tropical or temperate
clim_aus_des_tbl <- anti_join(data.frame(clim_aus_tbl), data.frame(clim_clip_tbl))
scores_AUS_des <- suprow(pca_env,clim_aus_des_tbl[,bioclims])$li

#Determine areas of NZ that overlap with Australian climate space (Used for exploratory analysis not in publication)
nz_clim_scores <- data.frame(clim_nz_tbl, scores_NZ, round(scores_NZ,1))
colnames(nz_clim_scores) <- c(colnames(nz_clim_scores)[1:21], "Axis1_round", "Axis2_round")
nz_aus <- inner_join(round(scores_AUS,1), nz_clim_scores, by = c("Axis1" = "Axis1_round", "Axis2" = "Axis2_round"))
nz_in_aus <- nz_aus[3:21]
scores_NZ <- suprow(pca_env,nz_in_aus[,bioclims])$li

clim_nz_aus <- rbind(clim_aus_tbl, nz_in_aus)
clim_nz_aus_pca <- dudi.pca(clim_nz_aus[,bioclims], scannf = FALSE, nf = 2)
clim_nz_pca <- dudi.pca(nz_in_aus[,bioclims], scannf = FALSE, nf = 2)

plot(scores_globclim)
plot(scores_trop)
plot(scores_temp)

climate_space <- ggplot() + 
  geom_point(data = scores_AUS, mapping = aes(Axis1, Axis2), col = "grey", fill = NA, alpha = 0.1, shape = 16, size = 1) +
  geom_point(data = scores_temp, mapping = aes(Axis1, Axis2), col = "#D81B60", fill = NA, alpha = 0.2, shape = 16, size = 1) +
  geom_point(data = scores_trop, mapping = aes(Axis1, Axis2), col = "#62E89D", fill = NA, alpha = 0.2, shape = 16, size = 1) +
  geom_point(data = scores_NZ, mapping = aes(Axis1, Axis2), col = "black", fill = NA, alpha = 0.2, shape = 16, size = 1) +
  theme_minimal()

ggsave(file = "./figures/climate_space.jpeg", climate_space, dpi = 600)

##Climate Space Figure ---------
clim_space_all <- ggplot() +
  geom_vline(xintercept = 0, col = "black")+
  geom_hline(yintercept = 0, col = "black")+
  stat_density_2d(data = scores_AUS_des, mapping = aes(Axis1, Axis2), fill = "#848484", alpha = 0.1, geom = "polygon", breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 1)) +
  stat_density_2d(data = scores_AUS_des, mapping = aes(Axis1, Axis2), col = "#848484",  breaks = c(0.001, 1), size = 0.3, linetype = "dashed") +
  stat_density_2d(data = scores_NZ, mapping = aes(Axis1, Axis2), fill = "#464646", alpha = 0.1, geom = "polygon", breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 1)) +
  stat_density_2d(data = scores_NZ, mapping = aes(Axis1, Axis2), col = "#464646",  breaks = c(0.001, 1), size = 0.3, linetype = "dashed") +
  stat_density_2d(data = scores_temp, mapping = aes(Axis1, Axis2), fill = "#C285A3", alpha = 0.3, geom = "polygon", breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 1)) +
  stat_density_2d(data = scores_temp, mapping = aes(Axis1, Axis2), col = "#C285A3",  breaks = c(0.001, 1), size = 0.3, linetype = "dashed") +
  stat_density_2d(data = scores_trop, mapping = aes(Axis1, Axis2), fill = "#92D050", alpha = 0.3, geom = "polygon", breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 1)) +
  stat_density_2d(data = scores_trop, mapping = aes(Axis1, Axis2), col = "#92D050",  breaks = c(0.001, 1), size = 0.3, linetype = "dashed") +
  stat_density_2d(data = scores_globclim, mapping = aes(Axis1, Axis2), col = "black",  breaks = c(0.0001, 1), size = 0.35) +
  xlim(-7.5, 7.5)+
  ylim(-10.5, 6.5)+
  xlab("PCA 1")+
  ylab("PCA 2")+
  theme_minimal()


#ggplot(scores_globclim, aes(Axis1, Axis2)) +stat_density_2d() +geom_point()

#, col = "#A05078"
#, col = "#69A12B"
ggsave("./figures/clim_space_all.jpeg", clim_space_all, dpi = 600,  width = 5, height = 4.5, units = "in")

#PCA Contribution Figure
pcacontrib_AUS <- data.frame(bios = c("bio1", "bio2","bio3","bio4","bio5","bio6","bio7","bio8",
                                      "bio9", "bio10","bio11","bio12","bio13","bio14","bio15","bio16",
                                      "bio17", "bio18","bio19"),
                             Comp1 = pca_env$co[,1],
                             Comp2 = pca_env$co[,2])
plot_contrib <- ggplot(pcacontrib_AUS)+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  geom_segment(aes(x = 0, y = 0, xend = Comp1, yend = Comp2), size = 1, lineend = "round", arrow = arrow(length = unit(2, "mm")))+
  geom_label_repel(aes(x = Comp1, y = Comp2, label = bios), size = 3.5, color = "black", force = 2, min.segment.length = 0.2) +
  xlab("Contribution to PCA 1 = 63.7%")+
  ylab("Contribution to PCA 2 = 20.3%")+
  theme_minimal()
ggsave(filename = "./figures/PCAContrib.jpeg", plot_contrib, height = 5, width = 5, units = "in", dpi = 400)

## Schoener's D and niche plots ---------------- 
#create function that compares two niches and saves the niche overlap figure
niche_plot <- function(species_1, 
                       species_2, 
                       biocl = bioclims, #list of bioclims to include
                       clim_pca = pca_env #climate data
) {
  scores_clim <- clim_pca$li
  scores_1 <- suprow(clim_pca,species_1[,biocl])$li
  scores_2 <- suprow(clim_pca,species_2[,biocl])$li
  
  
  grid_1 <- ecospat.grid.clim.dyn(glob=scores_clim,
                                  glob1=scores_clim,
                                  sp=scores_1, R=100,
                                  th.sp=0)
  
  grid_2 <- ecospat.grid.clim.dyn(glob=scores_clim,
                                  glob1=scores_clim,
                                  sp=scores_2, R=100,
                                  th.sp=0)
  
  #niche_dyn <- ecospat.niche.dyn.index(grid_agathis, grid_araucaria, intersection = 0.1)
  
  nplot <- ecospat.plot.niche.dyn(grid_1, grid_2, quant=0, interest=2,
                                  title= paste("Niche Overlap (Schoener's D = ", round(ecospat.niche.overlap(grid_1, grid_2, cor = TRUE)$D,2), ")"), name.axis1="PC1",
                                  name.axis2="PC2"
                                  , col.unf = "#FB8150"
                                    , col.exp = "#4E87F8"
                                    # , col.unf = "#D65F02"
                                  # , col.exp = "#056C9A"
                                  , col.stab = "#788D90" #"#807f7e"
                                    , colZ1 = "black"
                                    , colZ2 = "black"
                                    , transparency = 30
  )
  nplot
}



jpeg("./figures/Araucaria_Agathis.jpeg")
niche_plot(Araucaria, Agathis) 
dev.off()

jpeg("./figures/Argophyllum_Corokia.jpeg")
niche_plot(Argophyllum, Corokia)
dev.off()

jpeg("./figures/Sloanea_Elaeocarpus.jpeg")
niche_plot(Sloanea, Elaeocarpus)
dev.off()

jpeg("./figures/Sloanea_Aristotelia.jpeg")
niche_plot(Sloanea, Aristotelia)
dev.off()

jpeg("./figures/Mallotus_Euphorbia.jpeg")
niche_plot(Mallotus, Euphorbia)
dev.off()

jpeg("./figures/Acacia_Sophora.jpeg")
niche_plot(Acacia, Sophora)
dev.off()

jpeg("./figures/Caesalpinia_Sophora.jpeg")
niche_plot(Caesalpinia, Sophora)
dev.off()

jpeg("./figures/Cryptocarya_Beilschmiedia.jpeg")
niche_plot(Cryptocarya, Beilschmiedia)
dev.off()

jpeg("./figures/Cryptocarya_Litsea.jpeg")
niche_plot(Cryptocarya, Litsea)
dev.off()

jpeg("./figures/Ludwigia_Epilobium.jpeg")
niche_plot(Ludwigia, Epilobium)
dev.off()

jpeg("./figures/Microcachrys_Podocarpus.jpeg")
niche_plot(Microcachrys, Podocarpus)
dev.off()

jpeg("./figures/Cupaniopsis_Alectryon.jpeg")
niche_plot(Cupaniopsis, Alectryon)
dev.off()

jpeg("./figures/Mischocarpus_Alectryon.jpeg")
niche_plot(Mischocarpus, Alectryon)
dev.off()



#ecospat.shift.centroids(scores_araucaria, scores_agathis, scores_globclim, scores_globclim)

## Niche Volume ---------------
niche_D <- function(species_1, 
                    species_2, 
                    biocl = bioclims, #list of bioclims to include
                    clim_pca = pca_env, #climate data
                    r_num = 2 #amount of rounding
) {
  scores_clim <- clim_pca$li
  scores_1 <- suprow(clim_pca,species_1[,biocl])$li
  scores_2 <- suprow(clim_pca,species_2[,biocl])$li
  
  
  grid_1 <- ecospat.grid.clim.dyn(glob=scores_clim,
                                  glob1=scores_clim,
                                  sp=scores_1, R=100,
                                  th.sp=0)
  
  grid_2 <- ecospat.grid.clim.dyn(glob=scores_clim,
                                  glob1=scores_clim,
                                  sp=scores_2, R=100,
                                  th.sp=0)
  
  #niche_dyn <- ecospat.niche.dyn.index(grid_agathis, grid_araucaria, intersection = 0.1)
  
  round(ecospat.niche.overlap(grid_1, grid_2, cor = TRUE)$D,r_num)
}

aus_troptemp <- clim_clip_tbl

volume_AUS_D <- rbind(niche_D(Araucaria, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Agathis, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Argophyllum, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Corokia, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Sloanea, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Elaeocarpus, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Aristotelia, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Mallotus, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Euphorbia, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Acacia, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Sophora, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Caesalpinia, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Cryptocarya, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Beilschmiedia, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Litsea, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Ludwigia, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Epilobium, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Microcachrys, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Podocarpus, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Cupaniopsis, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Alectryon, aus_troptemp, clim_pca = clim_nz_aus_pca),
                      niche_D(Mischocarpus, aus_troptemp, clim_pca = clim_nz_aus_pca)
)



niche_vol <- data.frame(
  taxa = c("Araucaria", "Agathis", 
           "Argophyllum", "Corokia", 
           "Sloanea", "Elaeocarpus", "Aristotelia", 
           "Mallotus", "Euphorbia", 
           "Acacia", "Sophora", "Caesalpinia", 
           "Cryptocarya", "Beilschmiedia", "Litsea", 
           "Ludwigia", "Epilobium", 
           "Microcachrys", "Podocarpus", 
           "Cupaniopsis", "Alectryon", "Mischocarpus"),
  status = c("extinct", "extant",
             "extinct", "extant",
             "extinct", "extant", "extant",
             "extinct", "extant",
             "extinct",  "extant", "extinct",
             "extinct", "extant","extant",
             "extinct", "extant",
             "extinct", "extant",
             "extinct", "extant", "extinct"),
  family = c("Araucariaceae", "Araucariaceae", 
             "Argophyllaceae", "Argophyllaceae",
             "Elaeocarpaceae", "Elaeocarpaceae","Elaeocarpaceae",
             "Euphorbiaceae", "Euphorbiaceae", 
             "Fabaceae", "Fabaceae","Fabaceae",
             "Lauraceae", "Lauraceae", "Lauraceae", 
             "Onagraceae", "Onagraceae",
             "Podocarpaceae", "Podocarpaceae", 
             "Sapindaceae", "Sapindaceae", "Sapindaceae"),
  volume_AUS = volume_AUS_D,
  volume_NZ = volume_NZ_D
  
)

niche_vol %>% ggplot(aes(status, volume_NZ_D, col = status)) +geom_boxplot()
niche_vol %>% ggplot(aes(status, volume_AUS_D, col = status)) +geom_boxplot()
niche_vol %>% pivot_longer(cols = c(-taxa, -status, -family)) %>% ggplot(aes(value, fill = status)) +geom_histogram()+facet_wrap(status~name)

vol_plot <- ggplot(niche_vol, aes(x = taxa, y = volume, fill = status, col = status)) + 
  geom_bar(stat = "identity") + 
  facet_grid(. ~ family, scale= "free_x", space = "free_x") +
  #geom_text(aes(label=Species), vjust = -0.5, position = position_dodge(0.9), size = 3) +
  theme_minimal()+ 
  scale_fill_manual(values = c("#056C9A", "#D65F02"), labels = c("Extant", "Extinct"))+
  scale_color_manual(values = c("#056C9A", "#D65F02"), labels = c("Extant", "Extinct"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(size = 8),
        legend.position="bottom", legend.title = element_blank()) +
  ylab("Niche Volume")
ggsave(filename = "./figures/NicheVolBar.jpeg", vol_plot, width = 10, height = 5, units = "in") 

NicheVolEvEBar <- ggplot(niche_vol, aes(x=status, y= volume_AUS_D, fill = status)) + geom_boxplot() + 
  theme_minimal()+
  scale_fill_manual(values = c("#4E87F8", "#FB8150"), labels = c("Extant", "Extinct"))+
  scale_color_manual(values = c("#4E87F8", "#FB8150"), labels = c("Extant", "Extinct"))+
  # scale_fill_manual(values = c("#056C9A", "#D65F02"), labels = c("Extant", "Extinct"))+
  # scale_color_manual(values = c("#056C9A", "#D65F02"), labels = c("Extant", "Extinct"))+
  #annotate(geom="text", x = 0, y=0.57, label="p-value = 0.47", color="black", hjust = 0, size = 3) +
  theme(legend.position = "none")+
  xlab("") +
  ylab("Niche Volume in Australia") +
  scale_x_discrete(labels=c("extinct" = "Extinct", "extant" = "Extant"))
ggsave(filename = "figures/NicheVolEvEBar.jpeg", NicheVolEvEBar, width = 3, height = 3, units = "in") 

NicheVolEvEBar_nz <- ggplot(niche_vol, aes(x=status, y= volume_NZ_D, fill = status)) + geom_boxplot() + 
  theme_minimal()+
  scale_fill_manual(values = c("#056C9A", "#D65F02"), labels = c("Extant", "Extinct"))+
  scale_color_manual(values = c("#056C9A", "#D65F02"), labels = c("Extant", "Extinct"))+
  annotate(geom="text", x = 0, y=0.33, label="p-value = 0.26", color="black", hjust = 0, size = 3) +
  theme(legend.position = "none")+
  xlab("") +
  ylab("Niche Volume in New Zealand Overlap") +
  scale_x_discrete(labels=c("extinct" = "Extinct", "extant" = "Extant"))
ggsave(filename = "figures/NicheVolEvEBar_nz.jpeg", NicheVolEvEBar_nz, width = 3, height = 3, units = "in") 

NicheVolExtinct <-niche_vol %>% dplyr::filter(status == "extinct") 
NicheVolExtant <- niche_vol %>% dplyr::filter(status == "extant") 
t.test(NicheVolExtinct$volume_AUS, NicheVolExtant$volume_AUS, alternative = "two.sided")
t.test(NicheVolExtinct$volume_NZ, NicheVolExtant$volume_NZ, alternative = "two.sided")


vol_difference <- data.frame(
  pair = taxa_pairs,
  extinct = c(
    niche_vol$volume_AUS[which(niche_vol$taxa == "Araucaria")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Argophyllum")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Sloanea")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Sloanea")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Mallotus")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Acacia")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Caesalpinia")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Cryptocarya")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Cryptocarya")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Ludwigia")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Microcachrys")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Cupaniopsis")],
    niche_vol$volume_AUS[which(niche_vol$taxa == "Mischocarpus")]),
  extant = c(
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Agathis")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Corokia")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Elaeocarpus")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Aristotelia")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Euphorbia")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Sophora")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Sophora")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Beilschmiedia")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Litsea")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Epilobium")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Podocarpus")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Alectryon")]),
    (niche_vol$volume_AUS[which(niche_vol$taxa == "Alectryon")])
  )
)

vol_diff_points_aus <- ggplot(vol_difference, aes(extinct, extant))+geom_point()+ 
  geom_smooth(method='lm', formula= y~x, se = FALSE, col = "black") +
  xlab("Extinct Niche Volume") +
  ylab("Extant Niche Volume") +
  annotate(geom="text", x =0.4, y=0.2, label="y = 0.1x + 0.3", color="black", hjust = 0, size = 3) +
  annotate(geom="text", x=0.4, y=0.17, label="p-value = 0.72", color="black", hjust = 0, size = 3) +
  annotate(geom="text", x=0.4, y=0.14, label= paste("r^2 == ", 0.01), color="black", parse = TRUE, hjust = 0, size = 3) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(), 
        panel.background = element_blank())

ggsave("./figures/vol_diff_points_aus.jpeg", vol_diff_points_aus, width = 3, height = 3, units = "in")

summary(lm(vol_difference$extinct~vol_difference$extant))

# Niche similarity ----------------
#Calculate niche similarity between genera pairs
niche_similarity <- function(species_1, #extant
                             species_2, #extinct
                             biocl = bioclims, #list of bioclims to include
                             clim_pca = pca_env #climate data
) {
  scores_clim <- clim_pca$li
  scores_1 <- suprow(clim_pca,species_1[,biocl])$li
  scores_2 <- suprow(clim_pca,species_2[,biocl])$li
  
  
  grid_1 <- ecospat.grid.clim.dyn(glob=scores_clim,
                                  glob1=scores_clim,
                                  sp=scores_1, R=100,
                                  th.sp=0)
  
  grid_2 <- ecospat.grid.clim.dyn(glob=scores_clim,
                                  glob1=scores_clim,
                                  sp=scores_2, R=100,
                                  th.sp=0)
  
  extant_extinct <- ecospat.niche.similarity.test(grid_1, grid_2,rep=1000,
                                                  overlap.alternative = "higher",
                                                  expansion.alternative = "lower",
                                                  stability.alternative = "higher",
                                                  unfilling.alternative = "lower",
                                                  intersection = NA,
                                                  rand.type=2)
  
  extinct_extant<- ecospat.niche.similarity.test(grid_2, grid_1,rep=1000,
                                                 overlap.alternative = "higher",
                                                 expansion.alternative = "lower",
                                                 stability.alternative = "higher",
                                                 unfilling.alternative = "lower",
                                                 intersection = NA,
                                                 rand.type=2)
  
  sim <- c(extant_extinct$p.D, extinct_extant$p.D)
  sim
}


taxa_pairs <- c("Araucaria-Agathis", "Argophyllum-Corokia", "Sloanea-Elaeocarpus", "Sloanea-Aristotelia", 
                "Mallotus-Euphorbia", "Acacia-Sophora", "Caesalpinia-Sophora", "Cryptocarya-Beilschmiedia",
                "Cryptocarya-Litsea", "Ludwigia-Epilobium", "Microcachrys-Podocarpus", 
                "Cupaniopsis-Alectryon", "Mischocarpus-Alectryon")

sim_tests <- rbind(niche_similarity(Araucaria, Agathis),  
                   niche_similarity(Argophyllum, Corokia),
                   niche_similarity(Sloanea, Elaeocarpus),
                   niche_similarity(Sloanea, Aristotelia),
                   niche_similarity(Mallotus, Euphorbia),
                   niche_similarity(Acacia, Sophora),
                   niche_similarity(Caesalpinia, Sophora),
                   niche_similarity(Cryptocarya, Beilschmiedia),
                   niche_similarity(Cryptocarya, Litsea),
                   niche_similarity(Ludwigia, Epilobium),
                   niche_similarity(Microcachrys, Podocarpus),
                   niche_similarity(Cupaniopsis, Alectryon),
                   niche_similarity(Mischocarpus, Alectryon)
)

niche_similarity <- data.frame(
  pair = taxa_pairs,
  extant_extinct = sim_tests[,1],
  extinct_extant = sim_tests[,2]
  
)

# Individual bioclim comparison ------------------
a <- mutate(data.frame(Araucaria), Genus = "Araucaria", Family = "Araucariaceae", Status = "extinct", Pair = "Araucaria vs.Agathis")
b <- mutate(data.frame(Agathis), Genus = "Agathis", Family = "Araucariaceae", Status = "extant", Pair = "Araucaria vs.Agathis")
c <- mutate(data.frame(Argophyllum), Genus = "Argophyllum", Family = "Argophyllaceae", Status = "extinct", Pair = "Argophyllum vs. Corokia")
d <- mutate(data.frame(Corokia), Genus = "Corokia", Family = "Argophyllaceae", Status = "extant", Pair = "Argophyllum vs. Corokia")
e <- mutate(data.frame(Sloanea), Genus = "Sloanea", Family = "Elaeocarpaceae", Status = "extinct", Pair = "Sloanea vs. Elaeocarpus")
f <- mutate(data.frame(Elaeocarpus), Genus = "Elaeocarpus", Family = "Elaeocarpaceae", Status = "extant", Pair = "Sloanea vs. Elaeocarpus")
g <- mutate(data.frame(Sloanea), Genus = "Sloanea", Family = "Elaeocarpaceae", Status = "extinct", Pair = "Sloanea vs. Aristotelia")
h <- mutate(data.frame(Aristotelia), Genus = "Aristotelia", Family = "Elaeocarpaceae", Status = "extant", Pair = "Sloanea vs. Aristotelia")
i <- mutate(data.frame(Mallotus), Genus = "Mallotus", Family = "Euphorbiaceae", Status = "extinct", Pair = "Mallotus vs. Euphorbia")
j <- mutate(data.frame(Euphorbia), Genus = "Euphorbia", Family = "Euphorbiaceae", Status = "extant", Pair = "Mallotus vs. Euphorbia")
k <- mutate(data.frame(Acacia), Genus = "Acacia", Family = "Fabaceae", Status = "extinct", Pair = "Acacia vs. Sophora")
l <- mutate(data.frame(Sophora), Genus = "Sophora", Family = "Fabaceae", Status = "extant", Pair = "Acacia vs. Sophora")
m <- mutate(data.frame(Caesalpinia), Genus = "Caesalpinia", Family = "Fabaceae", Status = "extinct", Pair = "Caesalpinia vs. Sophora")
n <- mutate(data.frame(Sophora), Genus = "Sophora", Family = "Fabaceae", Status = "extant", Pair = "Caesalpinia vs. Sophora")
o <- mutate(data.frame(Cryptocarya), Genus = "Cryptocarya", Family = "Lauraceae", Status = "extinct", Pair = "Cryptocarya vs. Beilschmiedia")
p <- mutate(data.frame(Beilschmiedia), Genus = "Beilschmiedia", Family = "Lauraceae", Status = "extant", Pair = "Cryptocarya vs. Beilschmiedia")
q <- mutate(data.frame(Cryptocarya), Genus = "Cryptocarya", Family = "Lauraceae", Status = "extinct", Pair = "Cryptocarya vs. Litsea")
r <- mutate(data.frame(Litsea), Genus = "Litsea", Family = "Lauraceae", Status = "extant", Pair = "Cryptocarya vs. Litsea")
s <- mutate(data.frame(Ludwigia), Genus = "Ludwigia", Family = "Onagraceae", Status = "extinct", Pair = "Ludwigia vs. Epilobium")
t <- mutate(data.frame(Epilobium), Genus = "Epilobium", Family = "Onagraceae", Status = "extant", Pair = "Ludwigia vs. Epilobium")
u <- mutate(data.frame(Microcachrys), Genus = "Microcachrys", Family = "Podocarpaceae", Status = "extinct", Pair = "Microcachrys vs. Podocarpus")
v <- mutate(data.frame(Podocarpus), Genus = "Podocarpus", Family = "Podocarpaceae", Status = "extant", Pair = "Microcachrys vs. Podocarpus")
w <- mutate(data.frame(Cupaniopsis), Genus = "Cupaniopsis", Family = "Sapindaceae", Status = "extinct", Pair = "Cupaniopsis vs. Alectryon")
x <- mutate(data.frame(Alectryon), Genus = "Alectryon", Family = "Sapindaceae", Status = "extant", Pair = "Cupaniopsis vs. Alectryon")
y <- mutate(data.frame(Mischocarpus), Genus = "Mischocarpus", Family = "Sapindaceae", Status = "extinct", Pair = "Mischocarpus vs. Alectryon")
z <- mutate(data.frame(Alectryon), Genus = "Alectryon", Family = "Sapindaceae", Status = "extant", Pair = "Mischocarpus vs. Alectryon")

Q1_dataLong_NicheComp <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,n,m,o,p,q,r,s,t,u,v,w,x,y,z)
rm(a,b,c,d,e,f,g,h,i,j,k,l,n,m,o,p,q,r,s,t,u,v,w,x,y,z)
Q1_dataLong_NicheComp$Pair <- factor(Q1_dataLong_NicheComp$Pair, levels = c("Araucaria vs.Agathis", "Argophyllum vs. Corokia", "Sloanea vs. Elaeocarpus", 
                                                                            "Sloanea vs. Aristotelia", "Mallotus vs. Euphorbia", "Acacia vs. Sophora", 
                                                                            "Caesalpinia vs. Sophora", "Cryptocarya vs. Beilschmiedia", "Cryptocarya vs. Litsea", 
                                                                            "Ludwigia vs. Epilobium", "Microcachrys vs. Podocarpus", "Cupaniopsis vs. Alectryon", 
                                                                            "Mischocarpus vs. Alectryon"))
Q1_dataLong_NicheComp$Genus <- factor(Q1_dataLong_NicheComp$Genus, levels = c("Araucaria", "Agathis", "Argophyllum", "Corokia", "Sloanea", "Elaeocarpus", 
                                                                              "Aristotelia", "Mallotus", "Euphorbia", "Acacia" , 
                                                                              "Caesalpinia", "Sophora", "Cryptocarya", "Beilschmiedia", "Litsea", 
                                                                              "Ludwigia", "Epilobium", "Microcachrys", "Podocarpus", "Cupaniopsis", 
                                                                              "Mischocarpus", "Alectryon"))
#Density Plot
#comparison of individual climate variable density 
dens_plot <- Q1_dataLong_NicheComp %>%
  filter(Pair == "Araucaria vs.Agathis" | Pair == "Sloanea vs. Aristotelia" | Pair ==  "Mallotus vs. Euphorbia" | Pair == "Ludwigia vs. Epilobium") %>%
  mutate(bio12cm = round(bio12/10,0), bio1c = (bio1/10), bio15dec = (bio15/100)) %>%
  gather(key = "bios", value = "value", -Genus, -Family, -Status, -Pair) %>%
  filter(bios == "bio1c" | bios  == "bio7" | bios == "bio12cm" | bios == "bio15dec") %>%
  mutate(bios = factor(bios, levels = c("bio1c", "bio7", "bio12cm", "bio15dec"))) %>%
  mutate(Pair = factor(Pair, levels = c("Sloanea vs. Aristotelia", "Ludwigia vs. Epilobium", "Mallotus vs. Euphorbia", "Araucaria vs.Agathis"))) %>%
  ggplot(aes(x = value, fill = Status, col = Status)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#4E87F8", "#FB8150"), labels = c("Extant", "Extinct"))+
  scale_color_manual(values = c("#4E87F8", "#FB8150"), labels = c("Extant", "Extinct"))+
  facet_wrap(Pair~bios, scales = "free", labeller = label_wrap_gen(width = 15)) +
  theme_minimal()+
  theme(strip.text = element_blank(), axis.title = element_blank(), legend.position = "none", axis.text = element_text(size = 5))

ggsave(filename = "./figures/dens_plot.jpeg", dens_plot, width = 6, height = 5, units = "in", dpi = 600) 

#Calculate kernal density overlap
c(overlap(list(Araucaria[,1], Agathis[,1]), type = "2"),
  overlap(list(Araucaria[,7], Agathis[,7]), type = "2"),
  overlap(list(Araucaria[,12], Agathis[,12]), type = "2"),
  overlap(list(Araucaria[,15], Agathis[,15]), type = "2"))


c(overlap(list(Mallotus[,1], Euphorbia[,1]), type = "2"),
  overlap(list(Mallotus[,7], Euphorbia[,7]), type = "2"),
  overlap(list(Mallotus[,12], Euphorbia[,12]), type = "2"),
  overlap(list(Mallotus[,15], Euphorbia[,15]), type = "2"))

c(overlap(list(Ludwigia[,1], Epilobium[,1]), type = "2"),
  overlap(list(Ludwigia[,7], Epilobium[,7]), type = "2"),
  overlap(list(Ludwigia[,12], Epilobium[,12]), type = "2"),
  overlap(list(Ludwigia[,15], Epilobium[,15]), type = "2"))

c(overlap(list(Sloanea[,1], Aristotelia[,1]), type = "2"),
  overlap(list(Sloanea[,7], Aristotelia[,7]), type = "2"),
  overlap(list(Sloanea[,12], Aristotelia[,12]), type = "2"),
  overlap(list(Sloanea[,15], Aristotelia[,15]), type = "2"))


#Overlap vs. extinction age
age_overlap <- data.frame(extinct_age = c(2.5,16.5,16.5,16.5,2.5,1,28,16.5,16.5,24,1,2.5,2.5), 
                          overlap = c(0.33,0.37,0.75,0.09,0.44,0.32,0.72,0.83,0.85,0.06,0.24,0.64,0.52))
ggplot(age_overlap, aes(extinct_age, overlap))+ geom_point() + geom_smooth(method = "lm")

summary(lm(age_overlap$overlap~age_overlap$extinct_age))

