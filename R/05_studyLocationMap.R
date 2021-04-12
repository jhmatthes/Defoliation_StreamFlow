# Map a map of the watersheds with defoliation in 2015-2017
library(sf)
#library(spData)
#library(maps)
library(USAboundaries)

# Read gages point data
wshed_gages <- st_read("data/Gages_WatershedID.shp")
staid_use <- unique(FDC_data$STAID)
gages_staid <- dplyr::filter(wshed_gages, STAID %in% staid_use) %>%
  dplyr::arrange(STAID)

# Read watershed shapefile
wshed_shp <- st_read("../../MANUSCRIPT/Gage_Data/QGIS_Data/Watershed_Data/NHD_01_NorthEast/WBDHU12.shp") %>%
  dplyr::filter(HUC12 %in% gages_staid$HUC12)

wshed_join <- st_join(wshed_shp, gages_staid) %>%
  dplyr::mutate(defol_15 = X_15mean*-1, defol_16 = X_16mean*-1, defol_17 = X_17mean*-1)

# US state boundaries
data("us_states", package = "spData")
us_states_crs = st_transform(us_states, crs = raster::crs(gages_staid))

contemporary_ne <- us_states(states = c("Massachusetts", "Rhode Island",
                                        "Connecticut"))
contemporary_ne <- st_transform(contemporary_ne, crs = raster::crs(gages_staid))


# Get bounding box
wshed_bb = st_as_sfc(st_bbox(contemporary_ne))

# Make U.S. map with study area in red box 
ggm1 <- ggplot() + 
  geom_sf(data = us_states_crs, fill = "white") + 
  geom_sf(data = wshed_bb, fill = NA, color = "red", size = 0.75) +
  labs(subtitle = "A) Study Location in U.S.") +
  coord_sf(crs=2163) + 
  ggspatial::annotation_scale(height = unit(0.1,'cm'),aes(location="br")) +
  theme_void() +
  theme(legend.position = "none", plot.subtitle = element_text(face = "bold"))

# Defoliation by watershed for 2015, 2016, 2017
defol_range <- range(c(wshed_join$defol_15,wshed_join$defol_16,wshed_join$defol_17))

ggm2 <- ggplot() + 
  geom_sf(data = wshed_join, aes(fill = defol_15), size = 0.25) +
  geom_sf(data = contemporary_ne, fill = NA) +
  geom_sf(data = gages_staid, aes(shape = CLASS), size = 0.5, fill = "black") +
  scale_shape_manual(name = "Gage class",values = c(16, 25)) +
  rcartocolor::scale_fill_carto_c(limits = c(defol_range[1],defol_range[2]),
                     name = "Defol score", palette = "BrwnYl") +
  labs(subtitle = "B) 2015 Defoliation") +
  theme_void() +
  theme(legend.position = "none", plot.subtitle = element_text(face = "bold"))

ggm3 <- ggplot() + 
  geom_sf(data = wshed_join, aes(fill = defol_16), size = 0.25) +
  geom_sf(data = contemporary_ne, fill = NA) +
  geom_sf(data = gages_staid, aes(shape = CLASS), size = 0.5, fill = "black") +
  scale_shape_manual(name = "Gage class",values = c(16, 25)) +
  rcartocolor::scale_fill_carto_c(limits = c(defol_range[1],defol_range[2]),
                     name = "Defol score", palette = "BrwnYl") +
  labs(subtitle = "C) 2016 Defoliation") +
  theme_void() +
  theme(legend.position = "none", plot.subtitle = element_text(face = "bold"))

ggm4 <- ggplot() + 
  geom_sf(data = wshed_join, aes(fill = defol_17), size = 0.25) +
  geom_sf(data = contemporary_ne, fill = NA) +
  geom_sf(data = gages_staid, aes(shape = CLASS), size = 0.5, fill = "black") +
  scale_shape_manual(name = "Gage class",values = c(16, 25)) +
  rcartocolor::scale_fill_carto_c(limits = c(defol_range[1],defol_range[2]),
                     name = "Defoliation \n score", palette = "BrwnYl") +
  labs(subtitle = "D) 2017 Defoliation") +
  ggspatial::annotation_scale(height = unit(0.1,'cm'), aes(location="br")) +
  theme_void() +
  guides(shape = guide_legend(override.aes = list(size = 2))) +  
  theme(legend.spacing.x = unit(0.5, 'cm'),
        legend.position = "none", 
        plot.subtitle = element_text(face = "bold"))

map_leg <- cowplot::get_legend(ggm4 + 
                                 theme(legend.box.just = "left", legend.position="bottom"))
map_plots <- cowplot::plot_grid(ggm1, ggm2, ggm3, ggm4, nrow = 2, align = "vh") 
                              
#png("figures/map_v3.png", width = 3500, height = 2500, res = 600)
cowplot::plot_grid(map_plots, map_leg, 
                   nrow = 2, rel_heights = c(1.1, 0.2))
#dev.off()
