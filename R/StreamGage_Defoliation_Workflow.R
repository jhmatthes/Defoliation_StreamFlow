# This code is the main workflow for reproducing the analysis in: 
#   Smith-Tripp, S., A. Griffith, V. Pasquarella, J. H. Matthes, 'Impacts of a 
#   regional multi-year insect defoliation on seasonal water yield and instantaneous 
#   streamflow characteristics'.

library(dataRetrieval) 
library(tidyverse) 
library(lubridate)
library(daymetr)
library(cowplot)
library(lme4)

# Load the homegrown functions
source("R/Functions/download_15min_USGSgages.R") # Download USGS gage data
source("R/Functions/FDC_calc.R") # Calculate FDC curves
source("R/Functions/FDC_percentiles.R") # Find percentiles on FDC curves

# Make output directory, if not already there, to hold figures & tables
if(dir.exists("output/")){
  print("Will export figures and table to output/ folder in this working directory.") 
} else{
  dir.create("output/")
  print("Created output/ folder in this working directory to hold figures and tables.") 
}

# Load Landsat forest condition assessment product values 
# for in each subwatershed (HU12 code) that was calculated in QGIS
# Landsat defoliation data product: 10.5281/zenodo.1163679
allgages_defol <- read_csv("data/allgages_defol.csv")

# Load table with info for the stream gages that had 15-min data 
gages_dischargeDuration <- read_csv("data/streamGagesCoverage.csv") %>%
  left_join(filter(allgages_defol, year == 2016)) 

# Download 15-min USGS gage data for 1995-2017 to local data/Gage_Data directory
# This may take awhile (1-2 hours) depending on download speeds. 
download_15min_USGSgages(gages_dischargeDuration)

# Growing season window: June through September
start_month <- 6
end_month <- 9 

# Calcculate Seasonal Water Yield ------------------------------------------
# Plot and save each gage's FDC curve in a PDF file? (IF Y, SET FILE NAME BELOW)
plot_FDC <- FALSE

# Load and process 15-minute instantaneous stream gage data (previously downloaded locally):
# Aggregate data for seasonal budgets & calculate statistics for flow duration curves
for(g in 1:nrow(gages_dischargeDuration)){
  
  print(paste0("Working on ",g))
  
  # Read 15-min gage discharge data
  dat <- read.csv(paste0("data/Gage_Data/",gages_dischargeDuration$STAID[g],"_15min.csv"), header=T) 
  dat <- tidyr::separate(dat, dateTime, into = c("date","time"), sep = " ", remove = F) 
  dat$time <- substr(dat$time,1,8)
  dat <- mutate(dat, date = as.Date(date, format = "%Y-%m-%d"),
                time = hms(time),
                month = month(date), 
                year = year(date))
  
  # Calculate seasonal budgets of discharge for each stream gage x year
  watershed_area <- gages_dischargeDuration[g,]$DRAIN_SQKM
  gages_seasonalDischarge_site <- dat %>% 
    filter(month >= start_month & month <= end_month) %>%
    group_by(site_no, year) %>%
    summarize(discharge_seasonal = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE), #ft3/s to m3/s to m3/15min
              yield_seasonal = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE)/(mean(watershed_area,na.rm=TRUE)*1000^2), # yield in m (m3/m2) 
              sum_NA = sum(is.na(X_00060_00000)),
              sum_good = length(grep("A", X_00060_00000_cd))) %>%
    mutate(prop_good = sum_good/(4*24*(30+31+31+30)))
  
  # Calculate baseline (pre-2015) flow duration curve & pull percentile flow stats 
  dat_baseline <- filter(dat, date < "2015-01-01")
  baseline <- FDC_calc(dat_baseline$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  baseline_stats <- FDC_pertiles(baseline, gages_dischargeDuration$STAID[g], "baseline")
   
  # Add 2015, 2016, 2017, 2018 flow duration curves & pull percentile flow stats
  dat_2015 <- filter(dat, date >= "2015-01-01" & date <= "2015-12-31")
  FDC_2015 <- FDC_calc(dat_2015$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  F2015_stats <- FDC_pertiles(FDC_2015, gages_dischargeDuration$STAID[g], "2015")
  
  dat_2016 <- filter(dat, date >= "2016-01-01" & date <= "2016-12-31")
  FDC_2016 <- FDC_calc(dat_2016$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  F2016_stats <- FDC_pertiles(FDC_2016, gages_dischargeDuration$STAID[g], "2016")
  
  dat_2017 <- filter(dat, date >= "2017-01-01" & date <= "2017-12-31")
  FDC_2017 <- FDC_calc(dat_2017$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  F2017_stats <- FDC_pertiles(FDC_2017, gages_dischargeDuration$STAID[g], "2017")
  
  dat_2018 <- filter(dat, date >= "2018-01-01" & date <= "2018-12-31")
  FDC_2018 <- FDC_calc(dat_2018$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
  F2018_stats <- FDC_pertiles(FDC_2018, gages_dischargeDuration$STAID[g], "2018")
  
  # Combine all FDC stats for this site
  FDC_stats <- rbind(baseline_stats, F2015_stats, F2016_stats, 
                     F2017_stats, F2018_stats)
  
  if(g == 1){
    FDC_stats_sites <- FDC_stats 
    gages_seasonalDischarge <- gages_seasonalDischarge_site
  } else {
    FDC_stats_sites <- rbind(FDC_stats_sites, FDC_stats)
    gages_seasonalDischarge <- rbind(gages_seasonalDischarge, gages_seasonalDischarge_site)
  }
  
  if(plot_FDC == TRUE){
    pdf(paste0("FDC_plot_",dat$site_no[g],".pdf"))
    p1 <- ggplot(baseline) +
      geom_line(aes(x = exceed, y = q_dis)) +
      geom_hline(yintercept = median(baseline$exceed), lty=2) +
      geom_vline(xintercept = 0.5, lty=2) +
      geom_line(data = FDC_2015, aes(x = exceed, y = q_dis), color = "#E69F00") +
      geom_line(data = FDC_2016, aes(x = exceed, y = q_dis), color = "#56B4E9") +
      geom_line(data = FDC_2017, aes(x = exceed, y = q_dis), color = "#009E73") +
      geom_line(data = FDC_2018, aes(x = exceed, y = q_dis), color = "#0072B2") +
      scale_y_log10() +
      labs(y = expression(paste("Discharge (",m^{3},s^{-1},")")), 
           x = "Probability of exceeding discharge rate", 
           title = paste("Discharge Rate:",dat_2018$site_no[1]))
    print(p1)
    dev.off()
  }
}
 
# format gage station ID (STAID) column 
gages_seasonalDischarge <- gages_seasonalDischarge %>%
  mutate(STAID = str_pad(site_no, width = 8, side = "left", pad = "0")) 

# Make years with <95% good data NA
discharge_yrs_NA <- which(gages_seasonalDischarge$prop_good < 0.95)
gages_seasonalDischarge$discharge_seasonal[discharge_yrs_NA] <- NA
gages_seasonalDischarge$yield_seasonal[discharge_yrs_NA] <- NA

# Make table of coverage stats for min/max year and number of years per gage
gages_coverage <-  gages_seasonalDischarge %>%
  group_by(STAID) %>%
  summarize(yr_start = min(year),
            yr_end = max(year),
            duration = max(year)-min(year) +1,
            yrs_NA = sum(is.na(yield_seasonal)),
            good_years = duration - yrs_NA,
            prop_yrs_missing = yrs_NA/duration) 

# Filter gages with fewer than 10 years of baseline data
gages_seasonalDischarge <- left_join(gages_seasonalDischarge, gages_coverage, by = "STAID") %>%
  filter(good_years > 13) %>%
  ungroup(gages_seasonalDischarge) 

# Calculate Seasonal Precipitation ------------------------------------------
# Load set of stream gage locations (lat, lon) 
gage_locations <- read.csv("data/gagelocations.csv", header = T) %>%
  mutate(STAID_string = stringr::str_pad(STAID, width = 8, side = "left", pad = "0")) 

# Download Daymet data at the USGS stream gage locaitons
precipitationdata <- daymetr::download_daymet_batch(file_location = "data/gagelocations.csv",
                                                start = 1995, end = 2018, internal = T, silent = F)

# Calculate seasonal precipitation at each stream gage
for(g in 1:nrow(gage_locations)){
  
  # Match STAID to get correct drainage size
  gage_match <- which(allgages_defol$STAID==str_pad(precipitationdata[[g]]$site, width = 8, side = "left", pad = "0"))
  gage_drainage <- allgages_defol$DRAIN_SQKM[gage_match][1]
  
  gage_seasonalPrecip <- precipitationdata[[g]]$data %>%
    mutate(date = (as.Date(paste0(precipitationdata[[g]]$data$year,'-01-01')) + 
                     precipitationdata[[g]]$data$yday),
           month = lubridate::month(date)) %>%
    filter(month >= start_month & month <= end_month) %>%
    group_by(year) %>%
    summarize(precip_ann = sum(prcp..mm.day.,na.rm=TRUE)/1000, #convert precip mm to m
              num_na = sum(is.na(prcp..mm.day.))) %>%
    mutate(STAID = stringr::str_pad(precipitationdata[[g]]$site, width = 8, side = "left", pad = "0"),
           precip_gage = precip_ann*(gage_drainage*1000^2)) #m * m2(drainage area) = total m3 precip/gage
  
  # Aggregate seasonal gage precip into one dataframe
  if(g == 1){
    gages_seasonalPrecip <- gage_seasonalPrecip
  } else {
    gages_seasonalPrecip <- rbind(gages_seasonalPrecip, gage_seasonalPrecip)
  }
}

# Combine precipitation and discharge data & calculate yield-to-precip ratio
gages_dischargePrecip <- left_join(gages_seasonalDischarge, gages_seasonalPrecip,
                                   by = c("STAID", "year")) %>%
  mutate(yield_seasonal_mm = yield_seasonal*1000,
         precip_ann_mm = precip_ann*1000,
         prop_precip = yield_seasonal_mm/precip_ann_mm)


# Calculate Anomalies from Baseline Sesaonal Yield & Yield:Precip  --------------------------------------------
# Calculate the departures in precipitation & discharge from the 20-year mean (1995-2014)
dischargePrecip_baselineMean <- gages_dischargePrecip %>%
  filter(year < 2015) %>%
  group_by(STAID) %>%
  summarize(precip_mean = mean(precip_ann_mm, na.rm=TRUE),
            yield_mean = mean(yield_seasonal_mm, na.rm=TRUE),
            yieldratio_mean = mean(yield_seasonal_mm/precip_ann_mm, na.rm=TRUE)) %>%
  filter(yieldratio_mean < 2)

dischargePrecip <- left_join(dischargePrecip_baselineMean, gages_dischargePrecip, by = "STAID") %>%
  mutate(precip_norm = precip_ann_mm - precip_mean,
         yield_norm = yield_seasonal_mm - yield_mean,
         yieldratio_norm = yield_seasonal_mm/precip_ann_mm - yieldratio_mean)

# Combine yield anomalies with defoliation from Landsat data product
gages_defol_vals <- allgages_defol %>% 
  filter(STAID %in% unique(dischargePrecip$STAID)) %>% 
  mutate(defol_mean = defol_mean*-1)

dischargePrecip <- dischargePrecip %>% 
  left_join(gages_defol_vals, by = c("STAID", "year")) %>% 
  filter(STAID != '01105880', year >=2015 & year<2018,
         !is.na(yieldratio_norm))  #remove stream gage on cape cod that is not a single basin

# Plots and Linear Models -----------------------------------------------
#FIG 2: Seasonal Defoliation Boxplot by Year
dischargePrecip$ref_gage <- as.factor(dischargePrecip$ref_gage)
levels(dischargePrecip$ref_gage) <- c("Non-Ref", "Ref")
FIG2_defol_boxplot <- ggplot(dischargePrecip, aes(as.factor(year), defol_mean, 
                                             color = as.factor(year), shape = ref_gage, 
                                             size = ref_gage, group = year)) + 
  geom_hline(yintercept = 0, lty = 2) + 
  geom_boxplot(alpha = 0.2, outlier.colour = NA)+ 
  geom_point( alpha = 0.8, position = "jitter")+
  scale_size_discrete(range = c(2,3.5), name = NULL, labels = NULL,breaks = NULL) + 
  scale_color_manual(name = "Year", values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  scale_shape_discrete(name = "Gage class", breaks = c("Non-Ref", "Ref"), labels = c("Non-Reference", "Reference"))+
  #  facet_wrap(~ref_gage) + 
  labs(x = "Year", y = "Defoliation Metric")+
  theme_cowplot()

png("output/FIG2.png", width = 3500, height = 2000, res = 600)
FIG2_defol_boxplot
dev.off()

#FIG S1: Seasonal Precipitation Boxplot
FIGS1_precip_boxplot <- ggplot(filter(dischargePrecip, year >=2015),
                               aes(as.factor(year),precip_ann*1000, color = as.factor(year))) +
  geom_hline(yintercept = mean(filter(dischargePrecip, year < 2015)$precip_mean, na.rm = T), lty = 2) +
  geom_boxplot(alpha = 0.2, outlier.colour = NA)+
  geom_point(aes(shape = ref_gage), size = 1.5, alpha = 0.8, position = "jitter")+
  scale_color_manual(name = "Year", values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  scale_shape_discrete(name = "Class", breaks = c("FALSE", "TRUE"), labels = c("Non-Reference", "Reference"))+
  labs(x = "Year", y = "Total seasonal precipitation (mm)")+
  theme_cowplot()

png("output/FIGS1.png", width = 3500, height = 2000, res = 600)
FIGS1_precip_boxplot
dev.off()

#Yield Anomaly ~ Defoliation model
Ymod_all <- lmer(yield_norm ~ defol_mean + (1 | year),
                 data = dischargePrecip, 
                 control = lmerControl(optimizer ="Nelder_Mead"))
dischargePrecip$Ymod_all<- predict(Ymod_all) 

FIG3A_YAnom_Defol <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = yield_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = Ymod_all, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Water yield anomaly (mm)")+
  theme_cowplot()+
  theme(legend.position="none")

#Precip Anomaly ~ Defoliation Model
Precipmod_all <- lmer(precip_norm ~ defol_mean + (1 | year),
                      data = dischargePrecip)
dischargePrecip$Precipmod_all<- predict(Precipmod_all) 

FIG3B_PAnom_Defol <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = precip_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = Precipmod_all, group = as.factor(year), 
                color = as.factor(year)), size = 1, linetype = "dashed") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Precipitation anomaly (mm)")+
  theme_cowplot()+
  theme(legend.position="none")

#Yield:Precip Ratio ~ Defoliation Model
YPmod_all <- lmer(yieldratio_norm ~ defol_mean + (1 | year),
                  data = dischargePrecip, control = lmerControl(optimizer ="Nelder_Mead"))
dischargePrecip$YPmod_all<- predict(YPmod_all) 

FIG3C_YPRatio_Defol <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = yieldratio_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = YPmod_all, group = as.factor(year), color = as.factor(year)), 
            size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Yield:Precip anomaly (mm)", col = "Year")+
  theme_cowplot()+
  theme(legend.position="none")

FIG3_leg <- get_legend(FIG3C_YPRatio_Defol + theme(legend.position="bottom"))

FIG3_plots <- plot_grid(FIG3A_YAnom_Defol, FIG3B_PAnom_Defol, 
                        FIG3C_YPRatio_Defol, labels = c("A", "B", "C"), nrow = 1)

png("output/FIG3.png", width = 7000, height = 2500, res = 600)
plot_grid(FIG3_plots, FIG3_leg, ncol = 1, rel_heights = c(1.1, 0.2))
dev.off()

# REF GAGES ONLY: Yield anomaly, precip anomaly, yield:precip anomaly
# Filter just to reference gages
dischargePrecip_ref <- filter(dischargePrecip, ref_gage == "Ref")

# ref gages: Yield anomaly, precip anomaly, yield:precip anomaly
Ymod_ref <- lmer(yield_norm ~ defol_mean + (1 | year),
            data = dischargePrecip_ref, control = lmerControl(optimizer ="Nelder_Mead"))
dischargePrecip_ref$Ymod_ref<- predict(Ymod_ref) 

# Ref gages: Yield anomaly ~ Defol
FIGS2A <- ggplot(dischargePrecip_ref) +
  geom_point(aes(x = defol_mean, y = yield_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = Ymod_ref, group = as.factor(year), 
                color = as.factor(year)), size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Water yield anomaly (mm)")+
  theme_cowplot()+
  theme(legend.position="none")

#Ref gages: Precip anomaly ~ Defol
Precipmod_ref <- lmer(precip_norm ~ defol_mean + (1 | year),
            data = dischargePrecip_ref, control = lmerControl(optimizer ="Nelder_Mead"))
dischargePrecip_ref$Precipmod_ref<- predict(Precipmod_ref) 

FIGS2B <- ggplot(dischargePrecip_ref) +
  geom_point(aes(x = defol_mean, y = precip_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = Precipmod_ref, group = as.factor(year), 
                color = as.factor(year)), size = 1, linetype = "dashed") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Precipitation anomaly (mm)")+
  theme_cowplot()+
  theme(legend.position="none")

# Ref Gages: Yield:Precip Ratio ~ Defol
YPmod_ref <- lmer(yieldratio_norm ~ defol_mean + (1 | year),
            data = dischargePrecip_ref, control = lmerControl(optimizer ="Nelder_Mead"))
dischargePrecip_ref$YPmod_ref<- predict(YPmod_ref) 

FIGS2C <- ggplot(dischargePrecip_ref) +
  geom_point(aes(x = defol_mean, y = yieldratio_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = YPmod_ref, group = as.factor(year), 
                color = as.factor(year)), size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Yield:Precip anomaly (mm)", col = 'Year')+
  theme_cowplot()+
  theme(legend.position="none")

FIGS2_leg <- get_legend(FIGS2C + theme(legend.position="bottom"))
FIGS2_plots <- plot_grid(FIGS2A, FIGS2B, FIGS2C, labels = c("A", "B", "C"), nrow = 1)

png("output/FIGS2.png", width = 7000, height = 2500, res = 600)
plot_grid(FIGS2_plots, FIGS2_leg, ncol = 1, rel_heights = c(1, 0.2))
dev.off()

# Bind together model stats for output table 
year_effects <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3)),2),
                           response = c(rep("Yield Anom",6),rep("Y:Prcp Anom",6)),
                           year = rep(2015:2017,4),
                           intercept = signif(c(coef(Ymod_all)$year[,1],
                                               coef(Ymod_ref)$year[,1],
                                               coef(YPmod_all)$year[,1],
                                               coef(YPmod_ref)$year[,1]),3),
                           std_err_intercept = signif(c(rep(sqrt(diag(vcov(Ymod_all)))[1],3), 
                                             rep(sqrt(diag(vcov(Ymod_ref)))[1],3), 
                                             rep(sqrt(diag(vcov(YPmod_all)))[1],3), 
                                             rep(sqrt(diag(vcov(YPmod_ref)))[1],3)),3),
                           slope = signif(c(coef(Ymod_all)$year[,2],
                                           coef(Ymod_ref)$year[,2],
                                           coef(YPmod_all)$year[,2],
                                           coef(YPmod_all)$year[,2]),3), 
                           std_err_slope = signif(c(rep(sqrt(diag(vcov(Ymod_all)))[2],3), 
                                                   rep(sqrt(diag(vcov(Ymod_ref)))[2],3), 
                                                   rep(sqrt(diag(vcov(YPmod_all)))[2],3), 
                                                   rep(sqrt(diag(vcov(YPmod_ref)))[2],3)),3))

write.csv(year_effects, file = "output/Table1_YieldDefol_lmer.csv", row.names=FALSE)


# Calculate Flow Duration Curve 2015-2017 baseline departures ------------------------------------------
# Clean up FDC stats by site
FDC_stats_sites <- mutate(FDC_stats_sites, year = datasub)
FDC_stats_sites <- as_tibble(FDC_stats_sites)

# Calculate FDC defoliation year departure from baseline stats
FDC_discharge <- filter(FDC_stats_sites, FDC_type == "discharge",
                        STAID %in% unique(dischargePrecip$STAID))
sites <- unique(FDC_discharge$STAID)

# Find FDC 2015-2017 flow percent change from 10+ year FDC baseline conditions
# for 5, 25, 50, 75, and 95 flow exceedence percentiles
stats_columns <- grep("flow_", colnames(FDC_discharge))
for(g in 1:length(sites)){
  for(s in 1:length(stats_columns)){
    
    # Filter FDC stats by site 
    discharge_site <- filter(FDC_discharge, STAID == sites[g])
    
    # Baseline FDC stats
    baseline <- as.numeric(filter(discharge_site, year == "baseline")[,stats_columns[s]])
    
    # Calculate percent change between 2015-2017 FDC percentile flows and baseline FDC 
    diff_15 <- (as.numeric(filter(discharge_site, year == "2015")[,stats_columns[s]])-
                  baseline)/baseline
    
    diff_16 <- (as.numeric(filter(discharge_site, year == "2016")[,stats_columns[s]])-
                  baseline)/baseline
    
    diff_17 <- (as.numeric(filter(discharge_site, year == "2017")[,stats_columns[s]])-
                  baseline)/baseline
    
    # Combine site, stat, values for year diffs
    tmp <- tibble(STAID = as.character(sites[g]), 
                  stat_type = colnames(FDC_discharge[,stats_columns[s]]),
                  year = as.character(c(2015, 2016, 2017, 2018)),
                  flow_diff = c(diff_15, diff_16, diff_17))
    
    if(g == 1 & s == 1){
      FDC_departures <- tmp
    } else {
      FDC_departures <- rbind(FDC_departures, tmp)
    }
  }
}

# Join FDC stats departure data with defoliation
allgages_defol$year <- as.character(allgages_defol$year)
flow_diff_defol <- left_join(FDC_departures, allgages_defol, 
                             by = c("STAID", "year")) %>%
  filter(year >= 2015 & year < 2018) %>%
  mutate(defol_mean_pos = defol_mean * -1,
         flow_diff = flow_diff * 100)

# FDC Departure ~ Defoliation Models  --------------------------------------------

# 25% FDC exceedence (low probability of flow value) Model 
ALL_FDC25 <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                  data = filter(flow_diff_defol, stat_type == "flow_25"),
                  control = lmerControl(optimizer ="Nelder_Mead"))
flow_diff_FDC25 <- filter(flow_diff_defol, stat_type == "flow_25")
flow_diff_FDC25$pred_FDC25 <- predict(ALL_FDC25)

REF_FDC25 <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                  data = filter(flow_diff_defol, stat_type == "flow_25", ref_gage == "TRUE"), 
                  control = lmerControl(optimizer ="Nelder_Mead"))
flow_diff_FDC25_ref <- filter(flow_diff_defol, stat_type == "flow_25", ref_gage == "TRUE")
flow_diff_FDC25_ref$pred_FDC25 <- predict(REF_FDC25)

# 50% FDC exceedence (medium probability of flow value) Model
ALL_FDC50 <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                      data = filter(flow_diff_defol, stat_type == "flow_50"),
                      control = lmerControl(optimizer ="Nelder_Mead"))
flow_diff_FDC50 <- filter(flow_diff_defol, stat_type == "flow_50")
flow_diff_FDC50$pred_FDC50 <- predict(ALL_FDC50)

REF_FDC50 <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                      data = filter(flow_diff_defol, stat_type == "flow_50", ref_gage == "TRUE"), 
                      control = lmerControl(optimizer ="Nelder_Mead"))
flow_diff_FDC50_ref <- filter(flow_diff_defol, stat_type == "flow_50", ref_gage == "TRUE")
flow_diff_FDC50_ref$pred_FDC50 <- predict(REF_FDC50)

# 75% FDC exceedence (high probability of flow value) Model
ALL_FDC75 <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                  data = filter(flow_diff_defol, stat_type == "flow_75"),
                  control = lmerControl(optimizer ="Nelder_Mead"))
flow_diff_FDC75 <- filter(flow_diff_defol, stat_type == "flow_75")
flow_diff_FDC75$pred_FDC75 <- predict(ALL_FDC75)

REF_FDC75 <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                  data = filter(flow_diff_defol, stat_type == "flow_75", ref_gage == "TRUE"), 
                  control = lmerControl(optimizer ="Nelder_Mead"))
flow_diff_FDC75_ref <- filter(flow_diff_defol, stat_type == "flow_75", ref_gage == "TRUE")
flow_diff_FDC75_ref$pred_FDC75 <- predict(REF_FDC75)

# FIG 4: Plot changes in FDC percentiles compared to baseline FDC, all gages
# Plot changes in 25%tile against defoliation, by year
FIG4A <- ggplot(flow_diff_FDC25) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = pred_FDC25, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, " flow at 25% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 50%tile against defoliation, by year
FIG4B <- ggplot(flow_diff_FDC50) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = pred_FDC50, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", y = expression(paste("Percent ", Delta, " flow at 50% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 75%tile against defoliation, by year
FIG4C <- ggplot(flow_diff_FDC75) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = pred_FDC75, color = year), linetype = "dashed")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", y = expression(paste("Percent ", Delta, " flow at 75% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")

FIG4_leg <- get_legend(FIG4C + theme(legend.position="bottom"))
FIG4_plots <- cowplot::plot_grid(FIG4A, FIG4B, FIG4C, 
                                 labels = c("A", "B", "C"), nrow = 1)

png("output/FIG4.png", width = 7000, height = 2500, res = 600)
plot_grid(FIG4_plots, FIG4_leg, ncol = 1, rel_heights = c(1, 0.2))
dev.off()

# FIGS3: Reference plots FDC percentile changes with defoliation
# Plot changes in 25%tile against defoliation, by year
FIGS3A <- ggplot(flow_diff_FDC25_ref) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = pred_FDC25, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, " flow at 25% prob."))) +
  theme_cowplot() +
  theme(legend.position="none")

# Plot changes in 50%tile against defoliation, by year
FIGS3B <- ggplot(flow_diff_FDC50_ref) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = pred_FDC50, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", y = expression(paste("Percent ", Delta, " flow at 50% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 75%tile against defoliation, by year
FIGS3C <- ggplot(flow_diff_FDC75_ref) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = pred_FDC75, color = year), linetype = "dashed")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", y = expression(paste("Percent ", Delta, " flow at 75% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")

FIGS3_leg <- get_legend(FIGS3C + theme(legend.position="bottom"))
FIGS3_plots <- cowplot::plot_grid(FIGS3A, FIGS3B, FIGS3C, 
                                  labels = c("A", "B", "C"), nrow = 1)

png("output/FIGS3.png", width = 7000, height = 2500, res = 600)
plot_grid(FIGS3_plots, FIGS3_leg, ncol = 1, rel_heights = c(1, 0.2))
dev.off()

# Make table of FDC change ~ defol model output parameters
FDC_25_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = round(c(coef(ALL_FDC25)$year[,1], 
                                               coef(REF_FDC25)$year[,1]),2), 
                           std_err_intercept = round(c(rep(sqrt(diag(vcov(ALL_FDC25)))[1],3), 
                                                       rep(sqrt(diag(vcov(REF_FDC25)))[1],3)),2),
                           slope = round(c(coef(ALL_FDC25)$year[,2],
                                           coef(REF_FDC25)$year[,2]),2), 
                           std_err_slope = round(c(rep(sqrt(diag(vcov(ALL_FDC25)))[2],3),
                                                   rep(sqrt(diag(vcov(REF_FDC25)))[2],3)),2),
                           stat_type = rep("25tile",6))

FDC_50_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = round(c(coef(ALL_FDC50)$year[,1], 
                                               coef(REF_FDC50)$year[,1]),2), 
                           std_err_intercept = round(c(rep(sqrt(diag(vcov(ALL_FDC50)))[1],3), 
                                                       rep(sqrt(diag(vcov(REF_FDC50)))[1],3)),2),
                           slope = round(c(coef(ALL_FDC50)$year[,2],
                                           coef(REF_FDC50)$year[,2]),2), 
                           std_err_slope = round(c(rep(sqrt(diag(vcov(ALL_FDC50)))[2],3),
                                                   rep(sqrt(diag(vcov(REF_FDC50)))[2],3)),2),
                           stat_type = rep("50tile",6))

FDC_75_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = signif(c(coef(ALL_FDC75)$year[,1], 
                                                coef(REF_FDC75)$year[,1]),2), 
                           std_err_intercept = signif(c(rep(sqrt(diag(vcov(ALL_FDC75)))[1],3), 
                                                        rep(sqrt(diag(vcov(REF_FDC75)))[1],3)),2),
                           slope = signif(c(coef(ALL_FDC75)$year[,2],
                                            coef(REF_FDC75)$year[,2]),2), 
                           std_err_slope = signif(c(rep(sqrt(diag(vcov(ALL_FDC75)))[2],3),
                                                    rep(sqrt(diag(vcov(REF_FDC75)))[2],3)),2),
                           stat_type = rep("75tile",6))

FDC_lm_stats <- rbind(FDC_25_stats, FDC_50_stats) #drop 75 stats here because p > 0.05
write.csv(FDC_lm_stats, file = "output/Table2_FDCDefol_lmer.csv")
