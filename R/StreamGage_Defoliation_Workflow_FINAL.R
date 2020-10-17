# Code for flow duration curves and stats
#library(dataRetrieval) # for downloading stream gage data

# Functions from external libraries are called by functions where
# needed, except for magrittr, lme4, and lmerTest
#library(dplyr)
#library(tidyr)
#library(lubridate)
#library(cowplot)
library(magrittr)
library(daymetr)
library(sf)
library(lme4)
library(lmerTest)

# Load the flow duration curve functions
source("R/Functions/FDC_calc.R") # Calculate FDC curves
source("R/Functions/FDC_percentiles.R") # Find percentiles on FDC curves

# Stream gage dataset with metadata
gages_watershedID <- sf::read_sf("data/Gages_WatershedID.shp")
sf::st_geometry(gages_watershedID) <- NULL

# Find duplicate stream gages in the same subwatershed
duplicate_HUC12 <- gages_watershedID %>%
  dplyr::group_by(HUC12) %>%
  dplyr::summarize(count = dplyr::n()) %>%
  dplyr::mutate(dups_HUC12 = ifelse(count > 1, 1, 0)) %>%
  dplyr::filter(dups_HUC12 == 1)

# Get the STAID of smaller drainage areas to remove (keeping max drainage area gage per subwatershed)
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

# Only keep defoliation columns with mean 
stats_cols <- c(grep("mean", names(gages_nondups)))
gages_nondups <- gages_nondups[,c(1:3,5,stats_cols)]
name_cols <- which(grepl("mean", names(gages_nondups)))
names(gages_nondups)[name_cols] <- paste0(as.character(substr(names(gages_nondups)[name_cols], 3 , 10)))

# Gather and clean the columns into year, mean, stdev with rows by gage station
gages_defol <- tidyr::gather(gages_nondups, key = "year", value = "defol_mean", `2015mean`:`2017mean`) 
gages_defol$year <- as.numeric(substr(gages_defol$year,1,4))
gages_defol <- dplyr::rename(gages_defol, ref_gage = CLASS) 
gages_STAID <- data.frame(STAID = unique(gages_defol$STAID), 
                          DRAIN_SQKM = unique(gages_defol$DRAIN_SQKM)) %>%
  dplyr::arrange(STAID)


# #### Code used to download the 15-min instantaneous stream gage data
#  # Download 15-minute instantaneous discharge USGS stream gage data (just run once!)
# start_date <- "05-31"
# end_date <- "11-01"
# 
# for(g in 1:nrow(gages_dischargeDuration)){
#    years <- gages_dischargeDuration$yr_start[g]:gages_dischargeDuration$yr_end[g]
#   
#    for(y in 1:length(years)){
# 
#      tmp_yr <- readNWISuv(gages_dischargeDuration$STAID[g],
#                           parameterCd =  "00060",
#                           startDate = paste0(years[y],"-",start_date),
#                           endDate = paste0(years[y],"-",end_date))
#      
#      if(nrow(tmp_yr)!=0){
#        if(y == 1){
#          site_yrs <- tmp_yr
#        } else {
#          site_yrs <- rbind(site_yrs, tmp_yr)
#        }
#      }
#    }
#    write_csv(site_yrs, paste0("data/discharge_15min/",
#                               gages_dischargeDuration$STAID[g],"_15min.csv"))
# }

# Plot and save each gage's FDC curve in a PDF file? (NEED TO SET FILE NAME BELOW)
plot_FDC <- FALSE

# Window for calculating growing season metrics
start_month <- 6
end_month <- 9 

# Load and process 15-minute instantaneous stream gage data (previously downloaded locally):
# Aggregate data for seasonal budgets & calculate statistics for flow duration curves
for(g in 1:length(gages_STAID$STAID)){
  
  print(paste0("Working on ",g))
  
  # Read 15-min gage discharge data
  dat <- read.csv(paste0("../Gage_Data/discharge_15min/",gages_STAID$STAID[g],"_15min.csv"), header=T) 
  dat <- tidyr::separate(dat, dateTime, into = c("date","time"), sep = "T", remove = F) 
  dat$time <- substr(dat$time,1,8)
  dat <- dplyr::mutate(dat, date = as.Date(date, format = "%Y-%m-%d"),
                time = lubridate::hms(time),
                month = lubridate::month(date), 
                year = lubridate::year(date))
  
  # Calculate seasonal budgets of discharge for each stream gage x year
  watershed_area <- gages_STAID$DRAIN_SQKM[g]
  gages_monthlyDischarge_hires_site <- dat %>% 
    dplyr::mutate(time = lubridate::hms(lubridate::as_datetime(dateTime))) %>%
    dplyr::filter(month >= start_month & month <= end_month) %>%
    dplyr::mutate(month_year = paste0(month,"-",year)) %>%
    dplyr::group_by(site_no, month_year, month, year) %>%
    dplyr::summarize(discharge_monthly = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE), #ft3/s to m3/s to m3/15min per subwatershed
              yield_monthly = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE)/(mean(watershed_area)*1000^2), # yield in m (m3/m2) 
              sum_NA = sum(is.na(X_00060_00000)),
              sum_good = length(grep("A", X_00060_00000_cd))) %>%
    dplyr::mutate(days_good = sum_good/(4*24),
           days_good_frac = sum_good/(4*24)/lubridate::days_in_month(month))
  
  # Calculate baseline (pre-2015) flow duration curve & pull percentile flow stats 
  dat_baseline <- dplyr::filter(dat, date < "2015-01-01")
  baseline <- FDC_calc(dat_baseline$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  baseline_stats <- FDC_pertiles(baseline, gages_STAID$STAID[g], "baseline")
   
  # Add 2015, 2016, 2017, 2018 flow duration curves & pull percentile flow stats
  dat_2015 <- dplyr::filter(dat, date >= "2015-01-01" & date <= "2015-12-31")
  FDC_2015 <- FDC_calc(dat_2015$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  F2015_stats <- FDC_pertiles(FDC_2015, gages_STAID$STAID[g], "2015")
  
  dat_2016 <- dplyr::filter(dat, date >= "2016-01-01" & date <= "2016-12-31")
  FDC_2016 <- FDC_calc(dat_2016$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  F2016_stats <- FDC_pertiles(FDC_2016, gages_STAID$STAID[g], "2016")
  
  dat_2017 <- dplyr::filter(dat, date >= "2017-01-01" & date <= "2017-12-31")
  FDC_2017 <- FDC_calc(dat_2017$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  F2017_stats <- FDC_pertiles(FDC_2017, gages_STAID$STAID[g], "2017")
  
  #dat_2018 <- filter(dat, date >= "2018-01-01" & date <= "2018-12-31")
  #FDC_2018 <- FDC_calc(dat_2018$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  #F2018_stats <- FDC_pertiles(FDC_2018, gages_dischargeDuration$STAID[g], "2018")
  
  # Combine all FDC stats for this site
  FDC_stats <- rbind(baseline_stats, F2015_stats, F2016_stats, F2017_stats)
  
  if(g == 1){
    FDC_stats_sites <- FDC_stats 
    gages_monthlyDischarge_hires <- gages_monthlyDischarge_hires_site
  } else {
    FDC_stats_sites <- rbind(FDC_stats_sites, FDC_stats)
    gages_monthlyDischarge_hires <- rbind(gages_monthlyDischarge_hires, gages_monthlyDischarge_hires_site)
  }
  
  if(plot_FDC == TRUE){
    pdf(paste0("FDC_plot_",dat$site_no[g],".pdf"))
    p1 <- ggplot2::ggplot(baseline) +
      ggplot2::geom_line(ggplot2::aes(x = exceed, y = q_dis)) +
      ggplot2::geom_hline(yintercept = median(baseline$exceed), lty=2) +
      ggplot2::geom_vline(xintercept = 0.5, lty=2) +
      ggplot2::geom_line(data = FDC_2015, ggplot2::aes(x = exceed, y = q_dis), color = "#E69F00") +
      ggplot2::geom_line(data = FDC_2016, ggplot2::aes(x = exceed, y = q_dis), color = "#56B4E9") +
      ggplot2::geom_line(data = FDC_2017, ggplot2::aes(x = exceed, y = q_dis), color = "#009E73") +
      ggplot2::scale_y_log10() +
      ggplot2::labs(y = expression(paste("Discharge (",m^{3},s^{-1},")")), 
           x = "Probability of exceeding discharge rate", 
           title = paste("Discharge Rate:",dat_2017$site_no[1]))
    print(p1)
    dev.off()
  }
}

# Seasonal Budgets --------------------------------------------------------
# Make months with more than 90% missing data NA
gages_monthlyDischarge_hires$days_good_frac[gages_monthlyDischarge_hires$days_good_frac < 0.90] <- NA

# Summarize monthly to seasonal yearly discharge, keep years with > 95% non-NA data
gages_seasonalDischarge <- gages_monthlyDischarge_hires %>%
  dplyr::mutate(STAID = stringr::str_pad(as.factor(site_no),side="left",pad="0",width=8)) %>%
  dplyr::group_by(STAID, year)  %>% 
  dplyr::summarize(discharge_seasonal = sum(discharge_monthly), #m3/subwatershed
            yield_seasonal = sum(yield_monthly), #m at gage
            sum_NA = sum(is.na(discharge_monthly)),
            prop_good = sum(days_good)/(30+31+31+30), groups = "keep") %>%
  dplyr::filter(prop_good > 0.95) 

# Make table of coverage stats for min/max year and number of years per gage
gages_coverage <- gages_seasonalDischarge %>%
  dplyr::group_by(STAID) %>%
  dplyr::summarize(yr_start = min(year),
            yr_end = max(year),
            duration = max(year)-min(year) +1,
            sum_NA = sum(is.na(discharge_seasonal)),
            prop_yrs_missing = sum_NA/duration) 

# Filter gages with fewer than 10 years of baseline data
gages_seasonalDischarge <- dplyr::left_join(gages_seasonalDischarge, gages_coverage)  %>%
  dplyr::ungroup() %>% 
  dplyr::filter(duration > 14) %>%
  dplyr::mutate(STAID = stringr::str_pad(STAID, width = 8, side = "left", pad = "0"))

# Load Daymet Precipitation Data ------------------------------------------
# Load set of stream gage locations (lat, lon) 
gage_locations <- read.csv("../Gage_Data/gagelocations.csv", header = T) %>%
  dplyr::mutate(STAID_string = stringr::str_pad(STAID, width = 8, side = "left", pad = "0")) %>%
  dplyr::filter(STAID_string %in% gages_defol$STAID) 

readr::write_csv(gage_locations,"../Gage_Data/gagelocations_v2.csv")

# Download Daymet data at the USGS stream gage locaitons
precipitationdata <- daymetr::download_daymet_batch(file_location = "../Gage_Data/gagelocations_v2.csv",
                                                start = 1995, end = 2017, internal = T, silent = F)

# Calculate seasonal Daymet precipitation at each stream gage
for(g in 1:nrow(gage_locations)){
  gage_seasonalPrecip <- precipitationdata[[g]]$data %>%
    dplyr::mutate(date = (as.Date(paste0(precipitationdata[[g]]$data$year,'-01-01')) + 
                     precipitationdata[[g]]$data$yday),
           month = lubridate::month(date)) %>%
    dplyr::filter(month >= start_month & month <= end_month) %>%
    dplyr::group_by(year) %>%
    dplyr::summarize(precip_ann = sum(prcp..mm.day.)/1000, #convert precip mm to m
              num_na = sum(is.na(prcp..mm.day.))) %>%
    dplyr::mutate(STAID = stringr::str_pad(precipitationdata[[g]]$site, width = 8, side = "left", pad = "0"),
           precip_gage = precip_ann*(gages_STAID$DRAIN_SQKM[g]*1000^2)) #m * m2(drainage area) = total m3 precip/gage
  
  # gage_monthlyPrecip <- precipitationdata[[g]]$data %>%
  #   dplyr::mutate(date = (as.Date(paste0(precipitationdata[[g]]$data$year,'-01-01')) + 
  #                    precipitationdata[[g]]$data$yday),
  #          month = lubridate::month(date)) %>%
  #   dplyr::filter(month >= start_month & month <= end_month) %>%
  #   dplyr::group_by(month, year) %>%
  #   dplyr::summarize(precip_sum = sum(prcp..mm.day.)/1000, #convert mm to m
  #             num_na = sum(is.na(prcp..mm.day.))) %>%
  #   dplyr::mutate(STAID = stringr::str_pad(precipitationdata[[g]]$site, width = 8, side = "left", pad = "0"),
  #          precip_gage = precip_sum*(gages_STAID$DRAIN_SQKM[g]*1000^2)) #m * m2(drainage area) = m3 precip/gage
  # 
  # Aggregate seasonal gage precip into one dataframe
  if(g == 1){
    gages_seasonalPrecip <- gage_seasonalPrecip
    #gages_monthlyPrecip <- gage_monthlyPrecip
  } else {
    gages_seasonalPrecip <- rbind(gages_seasonalPrecip, gage_seasonalPrecip)
    #gages_monthlyPrecip <- rbind(gages_monthlyPrecip, gage_monthlyPrecip)
  }
}

# Combine precipitation and discharge data: calculate yield-to-precip ratio
gages_dischargePrecip <- dplyr::full_join(gages_seasonalPrecip, gages_seasonalDischarge, 
                                   by = c("STAID", "year")) %>%
  dplyr::mutate(runoff_ratio = yield_seasonal/precip_ann, na.rm=TRUE)

# Calculate the departures in precipitation & discharge from the 20-year mean (1995-2014)
dischargePrecip_mean <- gages_dischargePrecip %>%
  dplyr::filter(year < 2015) %>%
  dplyr::group_by(STAID) %>%
  dplyr::summarize(precip_mean = mean(precip_ann, na.rm=TRUE),
            precip_gage_mean = mean(precip_gage, na.rm=TRUE),
            discharge_mean = mean(discharge_seasonal, na.rm=TRUE),
            yield_mean = mean(yield_seasonal, na.rm=TRUE),
            runoff_ratio_mean = mean(runoff_ratio, na.rm=TRUE)) %>%
  dplyr::filter(runoff_ratio_mean < 2)

# Yield ratio calculates the departure of runoff ratio 
dischargePrecip <- dplyr::left_join(dischargePrecip_mean, gages_dischargePrecip, by = "STAID") %>%
  dplyr::mutate(precip_anom = precip_ann - precip_mean,
         discharge_anom = discharge_seasonal - discharge_mean,
         yield_anom = yield_seasonal - yield_mean,
         runoff_ratio_anom = runoff_ratio - runoff_ratio_mean)

gages_defol_vals <- gages_defol %>% 
  dplyr::filter(STAID %in% unique(dischargePrecip$STAID)) %>% 
  dplyr::mutate(defol_mean = defol_mean*-1)

dischargePrecip <- dischargePrecip %>% 
  dplyr::left_join(gages_defol_vals, by = c("STAID", "year")) %>% #make defoliation positive 
  dplyr::filter(STAID != '01105880', year >=2015 & year<2018,
                !is.na(runoff_ratio_anom))  #remove stream gage on cape cod, that is not a single basin

# Plots and Linear Models -----------------------------------------------
# FIG 2: Yearly Defoliation Boxplot
dischargePrecip$ref_gage <- as.factor(dischargePrecip$ref_gage)
levels(dischargePrecip$ref_gage) <- c("Non-Ref", "Ref")
#png("figures/FIG2.png", width = 4000, height = 2500, res = 600)
defol_boxplot <- ggplot(dischargePrecip, aes(as.factor(year), defol_mean, 
                            color = as.factor(year), shape = ref_gage, 
                            size = ref_gage, group = year)) + 
  geom_boxplot(alpha = 0.2, outlier.colour = NA, show.legend=FALSE) +
  geom_hline(yintercept = 0, lty = 2) + 
  geom_point( alpha = 0.8, position = "jitter")+
  scale_size_discrete(range = c(2,3.5), name = NULL, labels = NULL,breaks = NULL) + 
  scale_color_manual(name = "Year", values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4"))+
  scale_shape_discrete(name = "Gage class", breaks = c("Non-Ref", "Ref"), labels = c("Non-Reference", "Reference"))+
  labs(x = "Year", y = "Defoliation Metric")+
  cowplot::theme_cowplot()+
  theme(legend.key = element_rect(fill = NA, color = NA))
defol_boxplot
#dev.off()

# SUPP: Yearly Precipitation Boxplot
#png("figures/SuppPrecip.png", width = 4000, height = 2500, res = 600)
precip_boxplot <- ggplot(dplyr::filter(dischargePrecip, year >=2015),
                         aes(as.factor(year),precip_ann*1000, color = as.factor(year)))+
  geom_hline(yintercept = mean(dischargePrecip_mean$precip_mean)*1000, lty = 2) + 
  geom_boxplot(alpha = 0.2, outlier.colour = NA, show.legend=FALSE)+
  geom_point(aes(shape = ref_gage), size = 1.5, alpha = 0.8, position = "jitter")+
  scale_color_manual(name = "Year", values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  scale_shape_discrete(name = "Gage class", breaks = c("Non-Ref", "Ref"), labels = c("Non-Reference", "Reference"))+
  labs(x = "Year", y = "Precipitation (mm)")+
  cowplot::theme_cowplot()
precip_boxplot
#dev.off()

# FIG 3: Runoff Ratio ~ Defoliation, A) All Gages and B) Reference gages only
YPmod_all <- lmer(runoff_ratio_anom ~ defol_mean + (1 | year),
                  data = dischargePrecip, control = lmerControl(optimizer ="Nelder_Mead"))
summary(YPmod_all)
#YPmod_all_residuals = residuals(YPmod_all)
dischargePrecip$YPmod_all<- predict(YPmod_all) 

RunoffRatio_All <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = runoff_ratio_anom, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = YPmod_all, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Runoff ratio anomaly", col = "Year")+
  ylim(c(-0.25, 0.5))+
  xlim(c(-1, 2))+
  cowplot::theme_cowplot()+
  theme(legend.position="none")

# Runoff Ratio ~ defoliation; References gages only
ref_dischargePrecip <- dplyr::filter(dischargePrecip, ref_gage == "Ref")
YPmod_ref <- lmer(runoff_ratio_anom ~ defol_mean + (1 | year),
                  data = ref_dischargePrecip, 
                  control = lmerControl(optimizer ="Nelder_Mead"))
ref_dischargePrecip$YPmod_ref<- predict(YPmod_ref) 
summary(YPmod_ref)

RunoffRatio_Ref <- ggplot(ref_dischargePrecip) +
  geom_point(aes(x = defol_mean, y = runoff_ratio_anom, color = as.factor(year)), pch = 17) +
  geom_line(aes(x = defol_mean, y = YPmod_ref, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Runoff ratio anomaly", col = 'Year')+
  cowplot::theme_cowplot()+
  ylim(c(-0.25, 0.5)) +
  xlim(c(-1, 2)) +
  theme(legend.position="none")

RunoffRatio_leg <- cowplot::get_legend(RunoffRatio_All + theme(legend.position="bottom"))
RunoffRatio_Plot <- cowplot::plot_grid(RunoffRatio_All, RunoffRatio_Ref, labels = c("A", "B"))

RunoffRatio_Plot_Leg <- cowplot::plot_grid(RunoffRatio_Plot, RunoffRatio_leg, ncol = 1, rel_heights = c(1.1, 0.2))
#cowplot::ggsave2(filename = paste0(getwd(), "/figures/FIG3.png"), RunoffRatio_Plot_Leg,  width = 7.5, height = 4, dpi = 600)
#dev.off()

# Bind together LMM models outputs in table format
RRAnoms_YReffects <- data.frame(data = c(rep("All", 3),rep("Ref Only", 3)),
                           response = rep("RR Anom",6),
                           year = rep(2015:2017,2),
                           intercept = signif(c(coef(YPmod_all)$year[,1],
                                               coef(YPmod_ref)$year[,1]),2),
                           std_err_intercept = signif(c(rep(sqrt(diag(vcov(YPmod_all)))[1],3), 
                                             rep(sqrt(diag(vcov(YPmod_ref)))[1],3)),2),
                           slope = signif(c(coef(YPmod_all)$year[,2],
                                           coef(YPmod_ref)$year[,2]),2), 
                           std_err_slope = signif(c(rep(sqrt(diag(vcov(YPmod_all)))[2],3), 
                                                   rep(sqrt(diag(vcov(YPmod_ref)))[2],3)),2))

#write.csv(RRAnoms_YReffects, file = "figures/year_effects_inteceptmodel.csv")

# FLOW DURATION CURVE STATISTICS ------------------------------------------
# Clean up FDC stats by site
FDC_stats_sites <- dplyr::mutate(FDC_stats_sites, year = datasub)
#FDC_stats_sites <- as_tibble(FDC_stats_sites)

# Calculate FDC defoliation year departure from baseline stats
FDC_discharge <- dplyr::filter(FDC_stats_sites, FDC_type == "discharge")
sites <- unique(FDC_stats_sites$STAID)

# Get stats column indices for difference-from-baseline calculation
stats_columns <- grep("flow_", colnames(FDC_discharge))

# Calculate the FDC percentile departures in 2015-2017 from baseline
for(g in 1:length(sites)){
  for(s in 1:length(stats_columns)){
    
    # Filter FDC stats by site 
    discharge_site <- dplyr::filter(FDC_discharge, STAID == sites[g])
    
    # Baseline stat
    baseline <- as.numeric(dplyr::filter(discharge_site, year == "baseline")[,stats_columns[s]])
    
    # Calculate differences between 2015-2018 and baseline FDC stats
    diff_15 <- (as.numeric(dplyr::filter(discharge_site, year == "2015")[,stats_columns[s]])-
                  baseline)/baseline
    
    diff_16 <- (as.numeric(dplyr::filter(discharge_site, year == "2016")[,stats_columns[s]])-
                  baseline)/baseline
    
    diff_17 <- (as.numeric(dplyr::filter(discharge_site, year == "2017")[,stats_columns[s]])-
                  baseline)/baseline
    
    # Combine site, stat, values for year diffs
    tmp <- data.frame(STAID = as.character(sites[g]), 
                  stat_type = colnames(FDC_discharge)[stats_columns[s]],
                  year = as.character(c(2015, 2016, 2017)),
                  flow_diff = c(diff_15, diff_16, diff_17))
    
    if(g == 1 & s == 1){
      FDC_departures <- tmp
    } else {
      FDC_departures <- rbind(FDC_departures, tmp)
    }
  }
}

# Join defoliation with FDC stats difference data
gages_defol$year <- as.character(gages_defol$year)
flow_diff_defol <- dplyr::left_join(FDC_departures, gages_defol, 
                             by = c("STAID", "year")) %>%
  dplyr::filter(year >= 2015 & year < 2018) %>%
  dplyr::mutate(defol_mean = defol_mean*-1,
    flow_diff = flow_diff * 100)

# FDC percentile departures random effects models  --------------------------------------------
FDC_data <- tidyr::pivot_wider(flow_diff_defol, names_from = stat_type, values_from = flow_diff)

# All Gages 50% exceedence (medium probability of flow value) Model
ALL_FDC50_int <- lmer(flow_50 ~ defol_mean + (1 | year), 
                      data = FDC_data,
                      control = lmerControl(optimizer ="Nelder_Mead"))
summary(ALL_FDC50_int)
FDC_data$flow_50_model <- predict(ALL_FDC50_int)
#MuMIn::r.squaredGLMM(ALL_FDC50_int)

FDC_REFdata <- dplyr::filter(FDC_data, ref_gage == "Ref")
REF_FDC50_int <- lmer(flow_50 ~ defol_mean + (1 | year), 
                       data = FDC_REFdata, 
                      control = lmerControl(optimizer ="Nelder_Mead"))
summary(REF_FDC50_int)
FDC_REFdata$flow_50_model <- predict(REF_FDC50_int)

# All Gages 25% exceedence (low probability of flow value) Model
ALL_FDC25_int <- lmer(flow_25 ~ defol_mean + (1 | year), 
                      data = FDC_data)
FDC_data$flow_25_model <- predict(ALL_FDC25_int)

# Reference gages 25% exceedence (low probability of flow value) Model
REF_FDC25_int <- lmer(flow_25 ~ defol_mean + (1 | year), 
                      data = FDC_REFdata)
FDC_REFdata$flow_25_model <- predict(REF_FDC25_int)
summary(REF_FDC25_int)

# All Gages 75% exceedence (high probability of flow value) Model
ALL_FDC75_int <- lmer(flow_75 ~ defol_mean + (1 | year), 
                      data = FDC_data)
FDC_data$flow_75_model <- predict(ALL_FDC75_int)

# Reference Gages 75% exceedence (high probability of flow value) Model
REF_FDC75_int <- lmer(flow_75 ~ defol_mean + (1 | year), 
                      data = FDC_REFdata)
FDC_REFdata$flow_75_model <- predict(REF_FDC75_int)
summary(REF_FDC75_int)

FDC_data <- mutate(FDC_data, Year = year)
#### Plot defoliation ~ FDC percentile changes
# Plot changes in 25%tile against defoliation, by year
all_plot_25 <- ggplot(FDC_data) + 
  geom_point(aes( x = defol_mean, y = flow_25, color = year)) +
  geom_line(aes(x = defol_mean, y = flow_25_model, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 25% prob."))) +
  cowplot::theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 50%tile against defoliation, by year
all_plot_50 <- ggplot(FDC_data)  + 
  geom_point(aes( x = defol_mean, y = flow_50, color = year)) +
  geom_line(aes(x = defol_mean, y = flow_50_model, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 50% prob."))) +
  cowplot::theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 75%tile against defoliation, by year
all_plot_75 <- ggplot(FDC_data) + 
  geom_point(aes( x = defol_mean, y = flow_75, color = Year)) +
  geom_line(aes(x = defol_mean, y = flow_75_model, color = Year), linetype = "dashed")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 75% prob."))) +
  cowplot::theme_cowplot()+
  theme(legend.position="none")

all_FDC_leg <- cowplot::get_legend(all_plot_75 + theme(legend.position="bottom"))
top_row <- cowplot::plot_grid(all_plot_25, all_plot_50, all_plot_75, labels = c("A", "B", "C"), nrow = 1)

#png("figures/FIG4.png", width = 6000, height = 2500, res = 600)
cowplot::plot_grid(top_row, all_FDC_leg, ncol = 1, rel_heights = c(1.1, 0.2))
#dev.off()

# REFERENCE PLOTS ONLY
ref_plot_25 <- ggplot(FDC_REFdata) + 
  geom_point(aes( x = defol_mean, y = flow_25, color = year)) +
  geom_line(aes(x = defol_mean, y = flow_25_model, color = year))+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  #  facet_wrap(~year) +
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 25% prob."))) +
  cowplot::theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 50%tile against defoliation, by year
ref_plot_50 <- ggplot(FDC_REFdata)  + 
  geom_point(aes( x = defol_mean, y = flow_50, color = year)) +
  geom_line(aes(x = defol_mean, y = flow_50_model, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 50% prob."))) +
  cowplot::theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 75%tile against defoliation, by year
ref_plot_75 <- ggplot(FDC_REFdata) + 
  geom_point(aes( x = defol_mean, y = flow_75, color = year)) +
  geom_line(aes(x = defol_mean, y = flow_75_model, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 75% prob."))) +
  cowplot::theme_cowplot()+
  theme(legend.position="none")

ref_FDC_leg <- cowplot::get_legend(ref_plot_75 + theme(legend.position="bottom"))
top_row <- cowplot::plot_grid(ref_plot_25, ref_plot_50, ref_plot_75, labels = c("A", "B", "C"), nrow = 1)

#png("figures/SUPP_REFFDC.png", width = 6000, height = 2500, res = 600)
cowplot::plot_grid(top_row, ref_FDC_leg, ncol = 1, rel_heights = c(1.1, 0.2))
#dev.off()


# FDC Statistics Tables ---------------------------------------------------

#Bind LMM outputs in data.frame 
FDC_75_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = signif(c(coef(ALL_FDC75_int)$year[,1], 
                                               coef(REF_FDC75_int)$year[,1]),2), 
                           std_err_intercept = signif(c(rep(sqrt(diag(vcov(ALL_FDC75_int)))[1],3), 
                                                       rep(sqrt(diag(vcov(REF_FDC75_int)))[1],3)),2),
                           slope = signif(c(coef(ALL_FDC75_int)$year[,2],
                                         coef(REF_FDC75_int)$year[,2]),2), 
                           std_err_slope = signif(c(rep(sqrt(diag(vcov(ALL_FDC75_int)))[2],3),
                                   rep(sqrt(diag(vcov(REF_FDC75_int)))[2],3)),2),
                           stat_type = rep("75tile",6))
FDC_50_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = signif(c(coef(ALL_FDC50_int)$year[,1], 
                                               coef(REF_FDC50_int)$year[,1]),2), 
                           std_err_intercept = signif(c(rep(sqrt(diag(vcov(ALL_FDC50_int)))[1],3), 
                                                       rep(sqrt(diag(vcov(REF_FDC50_int)))[1],3)),2),
                           slope = signif(c(coef(ALL_FDC50_int)$year[,2],
                                           coef(REF_FDC50_int)$year[,2]),2), 
                           std_err_slope = signif(c(rep(sqrt(diag(vcov(ALL_FDC50_int)))[2],3),
                                                   rep(sqrt(diag(vcov(REF_FDC50_int)))[2],3)),2),
                           stat_type = rep("50tile",6))
FDC_25_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = signif(c(coef(ALL_FDC25_int)$year[,1], 
                                               coef(REF_FDC25_int)$year[,1]),2), 
                           std_err_intercept = signif(c(rep(sqrt(diag(vcov(ALL_FDC25_int)))[1],3), 
                                                       rep(sqrt(diag(vcov(REF_FDC25_int)))[1],3)),2),
                           slope = signif(c(coef(ALL_FDC25_int)$year[,2],
                                           coef(REF_FDC25_int)$year[,2]),2), 
                           std_err_slope = signif(c(rep(sqrt(diag(vcov(ALL_FDC25_int)))[2],3),
                                                   rep(sqrt(diag(vcov(REF_FDC25_int)))[2],3)),2),
                           stat_type = rep("25tile",6))
FDC_lm_stats <- rbind(FDC_25_stats, FDC_50_stats) #drop 75 stats here because the model is not well supported
#write.csv(FDC_lm_stats, file = "figures/FDC_lm_stats_int.csv")


# # SUPP: Yield Anomaly ~ Defoliation, Precip Anomaly ~ Defoliation
dischargePrecip <- filter(dischargePrecip, !is.na(defol_mean))
Ymod_all <- lmer(yield_anom ~ defol_mean + (1 | year),
                  data = dischargePrecip, 
                  control = lmerControl(optimizer ="Nelder_Mead"))
summary(Ymod_all)
dischargePrecip$Ymod_all<- predict(Ymod_all)
 
# Plot water yield anomalies ~ defoliation, w/ year random intercept
all_yieldanom <- ggplot(dischargePrecip) +
   geom_point(aes(x = defol_mean, y = yield_anom, color = as.factor(year))) +
   geom_line(aes(x = defol_mean, y = Ymod_all, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "solid") +
   scale_color_manual(name = "Year",values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4"))+
   labs(x = "Defoliation metric", y = "Water yield anomaly (m)")+
   cowplot::theme_cowplot()+
   theme(legend.position="none")
 
#Precip Anomaly Model
Precipmod_all <- lmer(precip_anom ~ defol_mean + (1 | year),
                       data = dischargePrecip)
summary(Precipmod_all)
dischargePrecip$Precipmod_all<- predict(Precipmod_all) 
 
all_precipanom <- ggplot(dischargePrecip) +
   geom_point(aes(x = defol_mean, y = precip_anom, color = as.factor(year))) +
   geom_line(aes(x = defol_mean, y = Precipmod_all, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "dashed") +
   scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
   labs(x = "Defoliation metric", y = "Precipitation anomaly (m)")+
   cowplot::theme_cowplot()+
   theme(legend.position="none")
 
all_3_leg <- cowplot::get_legend(all_yieldanom + theme(legend.position="bottom"))
allgages_plots <- cowplot::plot_grid(all_yieldanom, all_precipanom, labels = c("A", "B", "C"), nrow = 1)
#png("figures/SUPP_REFYanomPanom.png", width = 4000, height = 2500, res = 600)
cowplot::plot_grid(allgages_plots, all_3_leg, ncol = 1, rel_heights = c(1.1, 0.2))
#dev.off()

# # REF GAGES ONLY: Yield anomaly, precip anomaly, yield:precip anomaly
# # Filter just to reference gages
# ref_dischargePrecip <- filter(dischargePrecip, ref_gage == "Ref")
# 
# # ref GAGES: Yield anomaly, precip anomaly, yield:precip anomaly
# Ymod_ref <- lmer(yield_norm ~ defol_mean + (1 | year),
#                  data = ref_dischargePrecip, control = lmerControl(optimizer ="Nelder_Mead"))
# #ref_dischargePrecip$Ymod_ref<- predict(Ymod_ref) #cannot calculate prediction_refs with both standard errors and random effects
# 
# #Yield anomaly
# ref_1 <- ggplot(ref_dischargePrecip) +
#   geom_point(aes(x = defol_mean, y = yield_norm, color = as.factor(year))) +
#   #  geom_line(aes(x = defol_mean, y = Ymod_ref, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "solid") +
#   scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
#   labs(x = "Defoliation metric", y = "Water yield anomaly (m)")+
#   #facet_wrap(~year)+
#   theme_cowplot()+
#   theme(legend.position="none")
# 
# #Precip anomaly
# Precipmod_ref <- lmer(precip_norm ~ defol_mean + (1 | year),
#                       data = ref_dischargePrecip, control = lmerControl(optimizer ="Nelder_Mead"))
# ref_dischargePrecip$Precipmod_ref<- predict(Precipmod_ref) #cannot calculate prediction_refs with both standard errors and random effects
# 
# ref_2 <- ggplot(ref_dischargePrecip) +
#   geom_point(aes(x = defol_mean, y = precip_norm, color = as.factor(year))) +
#   #  geom_line(aes(x = defol_mean, y = Precipmod_ref, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "dashed") +
#   scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
#   labs(x = "Defoliation metric", y = "Precipitation anomaly (m)")+
#   #facet_wrap(~year)+
#   theme_cowplot()+
#   theme(legend.position="none")
# 
# 
# ref_3_leg <- get_legend(ref_3 + theme(legend.position="bottom"))
# 
# refgages_plots <- plot_grid(ref_1, ref_2, ref_3, labels = c("A", "B", "C"), nrow = 1)
# plot_grid(refgages_plots, ref_3_leg, ncol = 1, rel_heights = c(1, 0.2))

# Finding number of gages with different baseline period durations
gages_uniq <- unique(dischargePrecip$STAID)
tmp2 <- dplyr::filter(dischargePrecip, ref_gage == "Ref")
length(unique(tmp2$STAID))
tmp <- gages_coverage %>%
  dplyr::filter(STAID %in% gages_uniq) %>%
  dplyr::mutate(yrs_cat = dplyr::if_else(duration >= 21, 1, 
                                         dplyr::if_else(duration >= 18 & duration < 21, 2, 3)))

