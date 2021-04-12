# 01_loadData.R
# 1. Load stream gage information and watershed defoliation data processed by GIS.
# 2. Download 15-minute stream gage data (if needed)
# 3. Download Daymet precipitation data at stream gage locations

library(magrittr)

# Function to download high-res 15-min gage data
source("R/Functions/download_15min_USGSgages.R") 

# Load stream gage dataset with metadata
gages_watershedID <- sf::read_sf("data/Gages_WatershedID.shp")
sf::st_geometry(gages_watershedID) <- NULL

# Find duplicate stream gages in the same subwatershed
duplicate_HUC12 <- gages_watershedID %>%
  dplyr::group_by(HUC12) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  dplyr::mutate(dups_HUC12 = ifelse(count > 1, 1, 0)) %>%
  dplyr::filter(dups_HUC12 == 1)

# Get the STAID of smaller drainage areas to remove (keeping max drainage area per subwatershed)
for(w in 1:nrow(duplicate_HUC12)) {
  HUC12_match <- gages_watershedID[which(duplicate_HUC12$HUC12[w] == gages_watershedID$HUC12),]
  HUC12_remove <- HUC12_match[which(max(HUC12_match$DRAIN_SQKM) != HUC12_match$DRAIN_SQKM),]$STAID
  
  if(w == 1){
    HUC12_remove_dups <- HUC12_remove
  } else {
    HUC12_remove_dups <- c(HUC12_remove_dups, HUC12_remove)
  }
}

# Filter the duplicate gages in the same subwatershed to keep the one with the largest drainage
gages_nondups <- gages_watershedID %>%
  dplyr::filter(!(STAID %in% HUC12_remove_dups))

# Only keep defoliation columns with the mean 
stats_cols <- c(grep("mean", names(gages_nondups)))
gages_nondups <- gages_nondups[,c(1:3,5,stats_cols)]
name_cols <- which(grepl("mean", names(gages_nondups)))
names(gages_nondups)[name_cols] <- paste0(as.character(substr(names(gages_nondups)[name_cols], 3 , 10)))

# Gather and clean the columns into year, mean, stdev with rows by gage station
gages_defol <- tidyr::gather(gages_nondups, key = "year", value = "defol_mean", `2015mean`:`2017mean`) 
gages_defol$year <- as.numeric(substr(gages_defol$year,1,4))
gages_defol <- dplyr::rename(gages_defol, ref_gage = CLASS) 

# If needed: Download 15-min USGS gage data for 1995-2017 to local data/Gage_Data directory
# This may take awhile (1-2 hours) depending on download speeds and for this analysis
# will take up about ~2GB of space. 
#download_15min_USGSgages(gages_goodDischargeDuration)

# Load table with info for the stream gages that had good 15-min data 
gages_goodDischargeDuration <- readr::read_csv("data/gagelocations_v2.csv") %>%
  dplyr::mutate(STAID = stringr::str_pad(STAID, width = 8, side = "left", pad = "0"))

# Filter all subwatersheds to those with good 15-min gage data
gages_defol <- dplyr::filter(gages_defol, 
                             STAID %in% gages_goodDischargeDuration$STAID)

gages_STAID <- data.frame(STAID = unique(gages_defol$STAID), 
                          DRAIN_SQKM = unique(gages_defol$DRAIN_SQKM)) 

gages_goodDischargeDuration <- gages_goodDischargeDuration %>%
  dplyr::left_join(gages_STAID)

# Download Daymet data at the USGS stream gage locaitons
precipitationdata <- daymetr::download_daymet_batch(file_location = "data/gagelocations_v2.csv",
                                                    start = 1995, end = 2017, internal = T, silent = F)

# Remove intermediate variables that aren't used later on
# need gages_STAID$STAID gages_defol
rm(HUC12_remove, HUC12_remove_dups, HUC12_match, duplicate_HUC12, gages_nondups)

# Next up: 02_calculateBudgetsFDC.R

