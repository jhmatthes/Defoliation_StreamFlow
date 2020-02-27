# This code is the main workflow for reproducing the analysis in: 
#   Smith-Tripp, S., A. Griffith, V. Pasquarella, J. H. Matthes, 'Impacts of a 
#   regional multi-year insect defoliation on seasonal water yield and instantaneous 
#   streamflow characteristics'.
# Before running this code, you 


# Code for flow duration curves and stats
library(dataRetrieval) # for downloading stream gage data
library(tidyverse) 
library(lubridate)
library(daymetr)
library(cowplot)
library(lme4)

# Load the homegrown functions
source("R/Functions/download_15min_USGSgages.R") # Download USGS gage data
source("R/Functions/FDC_calc.R") # Calculate FDC curves
source("R/Functions/FDC_percentiles.R") # Find percentiles on FDC curves

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
  # gages_monthlyDischarge_site <- dat %>% 
  #   filter(month >= start_month & month <= end_month) %>%
  #   mutate(month_year = paste0(month,"-",year)) %>%
  #   group_by(site_no, month_year, month, year) %>%
  #   summarize(discharge_monthly = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE), #ft3/s to m3/s to m3/15min
  #             yield_monthly = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE)/(mean(watershed_area)*1000^2), # yield in m (m3/m2) 
  #             sum_NA = sum(is.na(X_00060_00000)),
  #             sum_good = length(grep("A", X_00060_00000_cd))) %>%
  #   mutate(prop_days_good = sum_good/(4*24),
  #          prop_good_frac = sum_good/(4*24)/days_in_month(month))
  
  gages_seasonalDischarge_site <- dat %>% 
    filter(month >= start_month & month <= end_month) %>%
    group_by(site_no, year) %>%
    summarize(discharge_seasonal = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE), #ft3/s to m3/s to m3/15min
              yield_seasonal = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE)/(mean(watershed_area)*1000^2), # yield in m (m3/m2) 
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

# Seasonal Budgets --------------------------------------------------------

# Gage data filtering: keep months < 90% missing data, years <95% missing data,
# gages with > 10 years baseline 15-min flow data.
#gages_monthlyDischarge$days_good_frac[gages_monthlyDischarge$days_good_frac < 0.90] <- NA

# # Sum monthly values to total seasonal discharge
# gages_seasonalDischarge <- gages_seasonalDischarge %>%
#   mutate(STAID = site_no) %>%
#   group_by(STAID, year) %>%
#   summarize(discharge_seasonal = sum(discharge_monthly, na.rm=TRUE),
#             yield_seasonal = sum(yield_monthly, na.rm=TRUE), 
#             sum_NA = sum(is.na(discharge_monthly)),
#             prop_good = sum(days_good)/(30+31+31+30)) #%>%
#   #filter(prop_good > 0.95) 

# Make years with <95% good data NA
gages_seasonalDischarge <- gages_seasonalDischarge %>%
  mutate(STAID = str_pad(site_no, width = 8, side = "left", pad = "0"))

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

gages_monthlyDischarge <- gages_monthlyDischarge %>%
  mutate(STAID = str_pad(site_no, width = 8, side = "left", pad = "0"))
 
# Load Daymet Precipitation Data ------------------------------------------
# Load set of stream gage locations (lat, lon) 
gage_locations <- read.csv("data/gagelocations.csv", header = T) %>%
  mutate(STAID_string = stringr::str_pad(STAID, width = 8, side = "left", pad = "0")) %>%
  filter(STAID_string %in% unique(gages_seasonalDischarge$STAID))

# Download Daymet data at the USGS stream gage locaitons
precipitationdata <- daymetr::download_daymet_batch(file_location = "data/gagelocations.csv",
                                                start = 1995, end = 2018, internal = T, silent = F)

# Calculate seasonal Daymet precipitation at each stream gage
for(g in 1:nrow(gage_locations)){
  gage_seasonalPrecip <- precipitationdata[[g]]$data %>%
    mutate(date = (as.Date(paste0(precipitationdata[[g]]$data$year,'-01-01')) + 
                     precipitationdata[[g]]$data$yday),
           month = lubridate::month(date)) %>%
    filter(month >= start_month & month <= end_month) %>%
    group_by(year) %>%
    summarize(precip_ann = sum(prcp..mm.day.)/1000, #convert precip mm to m
              num_na = sum(is.na(prcp..mm.day.))) %>%
    mutate(STAID = stringr::str_pad(precipitationdata[[g]]$site, width = 8, side = "left", pad = "0"),
           precip_gage = precip_ann*(allgages_defol$DRAIN_SQKM[g]*1000^2)) #m * m2(drainage area) = total m3 precip/gage
  
  gage_monthlyPrecip <- precipitationdata[[g]]$data %>%
    mutate(date = (as.Date(paste0(precipitationdata[[g]]$data$year,'-01-01')) + 
                     precipitationdata[[g]]$data$yday),
           month = lubridate::month(date)) %>%
    filter(month >= start_month & month <= end_month) %>%
    group_by(month, year) %>%
    summarize(precip_sum = sum(prcp..mm.day.)/1000, #convert mm to m
              num_na = sum(is.na(prcp..mm.day.))) %>%
    mutate(STAID = stringr::str_pad(precipitationdata[[g]]$site, width = 8, side = "left", pad = "0"),
           precip_gage = precip_sum*(allgages_defol$DRAIN_SQKM[g]*1000^2)) #m * m2(drainage area) = m3 precip/gage
  
  # Aggregate seasonal gage precip into one dataframe
  if(g == 1){
    gages_seasonalPrecip <- gage_seasonalPrecip
    gages_monthlyPrecip <- gage_monthlyPrecip
  } else {
    gages_seasonalPrecip <- rbind(gages_seasonalPrecip, gage_seasonalPrecip)
    gages_monthlyPrecip <- rbind(gages_monthlyPrecip, gage_monthlyPrecip)
  }
}

# Combine precipitation and discharge data: calculate yield-to-precip ratio
gages_dischargePrecip <- full_join(gages_seasonalPrecip, gages_seasonalDischarge, 
                                   by = c("STAID", "year")) %>%
  mutate(prop_precip = yield_seasonal/precip_ann)

# Calculate the departures in precipitation & discharge from the 20-year mean (1995-2014)
dischargePrecip_baselineMean <- gages_dischargePrecip %>%
  filter(year < 2015) %>%
  group_by(STAID) %>%
  summarize(prop_precip_mean = mean(prop_precip, na.rm=TRUE),
            precip_mean = mean(precip_ann, na.rm=TRUE),
            precip_gage_mean = mean(precip_gage, na.rm=TRUE),
            discharge_mean = mean(discharge_seasonal, na.rm=TRUE),
            yield_mean = mean(yield_seasonal, na.rm=TRUE),
            yieldratio_mean = mean(yield_seasonal/precip_ann, na.rm=TRUE)) #%>%
#  filter(prop_precip_mean < 2)

dischargePrecip <- left_join(dischargePrecip_baselineMean, gages_dischargePrecip, by = "STAID") %>%
  mutate(prop_precip_norm = prop_precip - prop_precip_mean,
         prop_precip_yield = yield_seasonal/precip_ann,
         precip_norm = precip_ann - precip_mean,
         discharge_norm = discharge_seasonal - discharge_mean,
         yield_norm = yield_seasonal - yield_mean,
         yieldratio_norm = yield_seasonal/precip_ann - yieldratio_mean)

# Combine yield anomalies with defoliation from Landsat data product
gages_defol_vals <- allgages_defol %>% 
  filter(STAID %in% unique(dischargePrecip$STAID)) %>% 
  mutate(defol_mean = defol_mean*-1)

dischargePrecip <- dischargePrecip %>% 
  left_join(gages_defol_vals, by = c("STAID", "year")) %>% 
  filter(STAID != '01105880', year >=2015 & year<2018,
         !is.na(yieldratio_norm))  #remove stream gage on cape cod, that is not a single basin

# Monthly yield & precip anomalies ----------------------------------------

# Departures in monthly discharge from the 20-year mean (1995-2014)
dischargeMonthly_mean <- gages_monthlyDischargePrecip %>%
  filter(year < 2015) %>%
  group_by(STAID, month) %>%
  summarize(prop_precip_mean = mean(prop_precip, na.rm=TRUE),
            precip_mean = mean(precip_gage, na.rm=TRUE),
            precip_mean_mm = mean(precip_sum*1000, na.rm=TRUE),
            discharge_mean = mean(discharge_monthly, na.rm=TRUE),
            yield_mean = mean(yield_monthly, na.rm=TRUE),
            yieldratio_mean = mean(yield_monthly/precip_sum, na.rm=TRUE)) %>%
  filter(abs(yield_mean) < 500)

#Monthly statistics of yield and precip departures from baseline
# Yield ratio calculates the departure of mean water yield that is theoretically independent of total precipitation amount 
dischargeMonthlyPrecip <- left_join(dischargeMonthly_mean, gages_monthlyDischargePrecip, by = "STAID") %>%
  mutate(prop_precip_norm = prop_precip - prop_precip_mean,
         prop_precip_yield = yield_monthly/precip_sum,
         precip_norm = precip_gage - precip_mean,
         precip_norm_mm = precip_sum*1000 - precip_mean_mm,
         discharge_norm = discharge_monthly - discharge_mean,
         yield_norm = yield_monthly - yield_mean,
         yieldratio_norm = yield_monthly/precip_sum - yieldratio_mean,
         month = month.x) 

allgages_defol_vals <- allgages_defol %>% 
  filter(STAID %in% unique(dischargeMonthlyPrecip$STAID)) %>% 
  mutate(defol_mean = defol_mean*-1)
dischargeMonthlyPrecip <- dischargeMonthlyPrecip %>% 
  left_join(allgages_defol_vals, by = c("STAID", "year")) %>% #make defoliation positive (moved out of doing it QGIS)
  filter(STAID != '01105880') #remove stream gage on cape cod, that is not a single basin

monthly_YtoP <- ggplot(filter(dischargeMonthlyPrecip, year >= 2015)) +
  geom_point(aes(x = defol_mean, y = yieldratio_norm, color = as.factor(year))) +
  geom_smooth(aes(x = defol_mean, y = yieldratio_norm, color = as.factor(year)), method = "lm")+
  scale_color_manual(name = "Year" , values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Yield:Precip anomaly") +
  facet_wrap(~month)+
  theme_cowplot()+
  theme(legend.position="none")

# # Look at relationship between the yield:precip ratio
# # and the size of the drainage area
# ggplot(filter(dischargeMonthlyPrecip, DRAIN_SQKM < 100, 
#               year == 2017)) +
#   geom_point(aes(x = DRAIN_SQKM, y = yieldratio_norm, color = as.factor(year))) +
#   scale_x_log10() +
#   theme_cowplot()


# Plots and Linear Models -----------------------------------------------

## Code to add in correct model to plots
#temp dataset to include lmer model predictions (easisest to remove NAs for predicting)
# complete_cases <- function(query, d) { 
#   d <- filter(dischargePrecip, year >= 2015 & year < 2018)
#   cols <- c(grep(query, names(d)))
#   i <- complete.cases(d[,cols],d$defol_mean) #subset to present values for query of interest
#   d <- d[i,]
#   return(d)
# }
#Yearly Defoliation Boxplot
dischargePrecip$ref_gage <- as.factor(dischargePrecip$ref_gage)
levels(dischargePrecip$ref_gage) <- c("Non-Ref", "Ref")
defol_boxplot <- ggplot(dischargePrecip, aes(as.factor(year), defol_mean, 
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

mean_precip <- dischargePrecip %>% 
  filter(year < 2015)# %>% 
  mutate(mean = mean(precip_mean, na.rm = T)) 

#Yearly precipitation boxplot
precip_boxplot <- ggplot(filter(dischargePrecip, year >=2015),aes(as.factor(year),precip_ann*1000, color = as.factor(year))) +
  geom_hline(yintercept = mean(filter(dischargePrecip, year < 2015)$precip_mean, na.rm = T), lty = 2) +
  geom_boxplot(alpha = 0.2, outlier.colour = NA)+
  geom_point(aes(shape = ref_gage), size = 1.5, alpha = 0.8, position = "jitter")+
  scale_color_manual(name = "Year", values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  scale_shape_discrete(name = "Class", breaks = c("FALSE", "TRUE"), labels = c("Non-Reference", "Reference"))+
  labs(x = "Year", y = "Precipitation (mm)")+
  theme_cowplot()

#Yield anomaly
#yield_norm_prediction <- complete_cases('yield_norm', dischargePrecip)
Ymod_all <- lmer(yield_norm ~ defol_mean + (1 | year),
            data = dischargePrecip, 
            control = lmerControl(optimizer ="Nelder_Mead"))
summary(Ymod_all)
dischargePrecip$Ymod_all<- predict(Ymod_all) #cannot calculate predictions with both standard errors and random effects

all_1 <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = yield_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = Ymod_all, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Water yield anomaly (m)")+
  #facet_wrap(~year)+
  theme_cowplot()+
  theme(legend.position="none")

#Precip anomaly
#precip_norm_prediction <- complete_cases('precip_norm', dischargePrecip)
Precipmod_all <- lmer(precip_norm ~ defol_mean + (1 | year),
            data = dischargePrecip)
dischargePrecip$Precipmod_all<- predict(Precipmod_all) #cannot calculate predictions with both standard errors and random effects

all_2 <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = precip_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = Precipmod_all, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "dashed") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Precipitation anomaly (m)")+
  #facet_wrap(~year)+
  theme_cowplot()+
  theme(legend.position="none")

#Yield Ratio
#yieldratio_norm_prediction <- complete_cases('yieldratio_norm', dischargePrecip)
YPmod_all <- lmer(yieldratio_norm ~ defol_mean + (1 | year),
            data = dischargePrecip, control = lmerControl(optimizer ="Nelder_Mead"))
dischargePrecip$YPmod_all<- predict(YPmod_all) #cannot calculate predictions with both standard errors and random effects

all_3 <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = yieldratio_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = YPmod_all, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Yield:Precip anomaly (m)", col = "Year")+
  #facet_wrap(~year)+
  theme_cowplot()+
  theme(legend.position="none")

all_3_leg <- get_legend(all_3 + theme(legend.position="bottom"))

allgages_plots <- plot_grid(all_1, all_2, all_3, labels = c("A", "B", "C"), nrow = 1)
plot_grid(allgages_plots, all_3_leg, ncol = 1, rel_heights = c(1.1, 0.2))


# REF GAGES ONLY: Yield anomaly, precip anomaly, yield:precip anomaly
# Filter just to reference gages
ref_dischargePrecip <- dischargePrecip[dischargePrecip$STAID %in% ref_STAID,]

# ref GAGES: Yield anomaly, precip anomaly, yield:precip anomaly
#yield_norm_prediction_ref <- complete_cases('yield_norm', ref_dischargePrecip)
Ymod_ref <- lmer(yield_norm ~ defol_mean + (1 | year),
            data = ref_dischargePrecip, control = lmerControl(optimizer ="Nelder_Mead"))
ref_dischargePrecip$Ymod_ref<- predict(Ymod_ref) #cannot calculate prediction_refs with both standard errors and random effects

#Yield anomaly
ref_1 <- ggplot(ref_dischargePrecip) +
  geom_point(aes(x = defol_mean, y = yield_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = Ymod_ref, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Water yield anomaly (m)")+
  #facet_wrap(~year)+
  theme_cowplot()+
  theme(legend.position="none")

#Precip anomaly
#precip_norm_prediction_ref <- complete_cases('precip_norm', ref_dischargePrecip)
Precipmod_ref <- lmer(precip_norm ~ defol_mean + (1 | year),
            data = ref_dischargePrecip, control = lmerControl(optimizer ="Nelder_Mead"))
ref_dischargePrecip$Precipmod_ref<- predict(Precipmod_ref) #cannot calculate prediction_refs with both standard errors and random effects

ref_2 <- ggplot(ref_dischargePrecip) +
  geom_point(aes(x = defol_mean, y = precip_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = Precipmod_ref, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "dashed") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Precipitation anomaly (m)")+
  #facet_wrap(~year)+
  theme_cowplot()+
  theme(legend.position="none")

#Yield Ratio
#yieldratio_norm_prediction_ref <- complete_cases('yieldratio_norm', ref_dischargePrecip)
YPmod_ref <- lmer(yieldratio_norm ~ defol_mean + (1 | year),
            data = ref_dischargePrecip, control = lmerControl(optimizer ="Nelder_Mead"))
ref_dischargePrecip$YPmod_ref<- predict(YPmod_ref) #cannot calculate prediction_refs with both standard errors and random effects

ref_3 <- ggplot(ref_dischargePrecip) +
  geom_point(aes(x = defol_mean, y = yieldratio_norm, color = as.factor(year))) +
  geom_line(aes(x = defol_mean, y = YPmod_ref, group = as.factor(year), color = as.factor(year)), size = 1, linetype = "solid") +
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  labs(x = "Defoliation metric", y = "Yield:Precip anomaly (m)", col = 'Year')+
  #facet_wrap(~year)+
  theme_cowplot()+
  theme(legend.position="none")

ref_3_leg <- get_legend(ref_3 + theme(legend.position="bottom"))

refgages_plots <- plot_grid(ref_1, ref_2, ref_3, labels = c("A", "B", "C"), nrow = 1)
plot_grid(refgages_plots, ref_3_leg, ncol = 1, rel_heights = c(1, 0.2))

# Save Plots --------------------------------------------------------------
## save plots figure folder (both main and not main)
#Defoliation Box Plot
#save_plot("Graphs/Supp_Figures/defol_boxplot.jpg", defol_boxplot, base_height = 4, base_width = 6)
#Precipitation Box plot
#save_plot("Graphs/Supp_Figures/precip_boxplot.jpg", precip_boxplot, base_height = 4, base_width = 6)
#All Gages
#save_plot("Graphs/Main_Figures/allgages_plots_lme_int.jpg", plot_grid(allgages_plots, all_3_leg, ncol = 1, rel_heights = c(1.2, 0.1)), base_height = 3.56, base_width = 10)
#Reference Gages
#save_plot("Graphs/Main_Figures/refgages_plots_lme_int.jpg", refgages_plots, base_height = 3.5, base_width = 10)

# Bind together stats for output table 
# And Plot year effects for intercept and slope
year_effects <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3)),2),
                           response = c(rep("Yield Anom",6),rep("Y:Prcp Anom",6)),
                           year = rep(2015:2017,4),
                           intercept = round(c(coef(Ymod_all)$year[,1],
                                               coef(Ymod_ref)$year[,1],
                                               coef(YPmod_all)$year[,1],
                                               coef(YPmod_ref)$year[,1]),2),
                           std_err_intercept = round(c(rep(sqrt(diag(vcov(Ymod_all)))[1],3), 
                                             rep(sqrt(diag(vcov(Ymod_ref)))[1],3), 
                                             rep(sqrt(diag(vcov(YPmod_all)))[1],3), 
                                             rep(sqrt(diag(vcov(YPmod_ref)))[1],3)),2),
                           slope = round(c(coef(Ymod_all)$year[,2],
                                           coef(Ymod_ref)$year[,2],
                                           coef(YPmod_all)$year[,2],
                                           coef(YPmod_all)$year[,2]),2), 
                           std_err_slope = round(c(rep(sqrt(diag(vcov(Ymod_all)))[2],3), 
                                                   rep(sqrt(diag(vcov(Ymod_ref)))[2],3), 
                                                   rep(sqrt(diag(vcov(YPmod_all)))[2],3), 
                                                   rep(sqrt(diag(vcov(YPmod_ref)))[2],3)),2))

#write.csv(year_effects, file = "year_effects_inteceptmodel.csv")
#Cut graphs of yearly intercepts and slopes, located in prior workflow


# FLOW DURATION CURVE STATISTICS ------------------------------------------
# Clean up FDC stats by site
FDC_stats_sites <- mutate(FDC_stats_sites, year = datasub)
FDC_stats_sites <- as_tibble(FDC_stats_sites)

# Calculate FDC defoliation year departure from baseline stats
FDC_discharge <- filter(FDC_stats_sites, FDC_type == "discharge")
sites <- unique(FDC_stats_sites$STAID)

# Get stats column indices for difference calculation
stats_columns <- grep("flow_", colnames(FDC_discharge))

for(g in 1:length(sites)){
  for(s in 1:length(stats_columns)){
    
    # Filter FDC stats by site 
    discharge_site <- filter(FDC_discharge, STAID == sites[g])
    
    # Baseline stat
    baseline <- as.numeric(filter(discharge_site, year == "baseline")[,stats_columns[s]])
    
    # Calculate differences between 2015-2018 and baseline FDC stats
    diff_15 <- (as.numeric(filter(discharge_site, year == "2015")[,stats_columns[s]])-
                  baseline)/baseline
    
    diff_16 <- (as.numeric(filter(discharge_site, year == "2016")[,stats_columns[s]])-
                  baseline)/baseline
    
    diff_17 <- (as.numeric(filter(discharge_site, year == "2017")[,stats_columns[s]])-
                  baseline)/baseline
    
    diff_18 <- (as.numeric(filter(discharge_site, year == "2018")[,stats_columns[s]])-
                  baseline)/baseline      
    
    # Combine site, stat, values for year diffs
    tmp <- tibble(STAID = as.character(sites[g]), 
                  stat_type = colnames(FDC_discharge[,stats_columns[s]]),
                  year = as.character(c(2015, 2016, 2017, 2018)),
                  flow_diff = c(diff_15, diff_16, diff_17, diff_18))
    
    if(g == 1 & s == 1){
      FDC_departures <- tmp
    } else {
      FDC_departures <- rbind(FDC_departures, tmp)
    }
  }
}

# Join defoliation with FDC stats difference data
allgages_defol$year <- as.character(allgages_defol$year)
flow_diff_defol <- left_join(FDC_departures, allgages_defol, 
                             by = c("STAID", "year")) %>%
  filter(year >= 2015 & year < 2018) %>%
  mutate(defol_mean_pos = defol_mean * -1,
         flow_diff = flow_diff * 100)

# Plot changes in 25%tile against defoliation, by year
all_plot <- ggplot(filter(flow_diff_defol, stat_type == "flow_25"), 
       aes(x = defol_mean_pos, y = flow_diff)) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~year) +
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 25% exceedence")),
       title = "All Gages") +
  theme_cowplot()

ref_plot <- ggplot(filter(flow_diff_defol, stat_type == "flow_25", ref_gage == TRUE),
              aes(x = defol_mean_pos, y = flow_diff)) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~year) +
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 25% exceedence")),
       title = "Reference Gages") +
  theme_cowplot()

plot_grid(ref_plot, all_plot, nrow = 1)

# FDC Mixed Statistics Models  --------------------------------------------
#set up a dataframe to store models in 
prediction_models_ref <- data.frame(filter(flow_diff_defol, ref_gage == 'TRUE'  & stat_type == 'flow_05')$defol_mean_pos)
complete_cases_FDC <- function(query, d, flow, ref) {
  d <- filter(flow_diff_defol, year >= 2015 & year < 2018)
  cols <- c(grep(query, names(d)))
  i <-
    complete.cases(d[, cols]) #subset to present values for query of interest
  d <- d[i, ]
  d <- filter(d, stat_type == flow)
  d <- d[, c(1, 3, cols, 10)]
  if (ref == T) {
    d <- filter(d, ref_gage == TRUE)
  }
  return(d)
}
prediction_models_all <- complete_cases_FDC(query = 'defol_mean_pos', d = flow_diff_defol, flow = 'flow_05', ref = F)
prediction_models_ref <- complete_cases_FDC(query = 'defol_mean_pos', d = flow_diff_defol, flow = 'flow_05', ref = T)
# All Gages 50% exceedence (medium probability of flow value)
ALL_FDC50_int <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                      data = filter(flow_diff_defol, stat_type == "flow_50"),
                      control = lmerControl(optimizer ="Nelder_Mead"))

prediction_models_all$flow_50 <- predict(ALL_FDC50_int)
MuMIn::r.squaredGLMM(ALL_FDC50_int)
REF_FDC50_int <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                       data = filter(flow_diff_defol, stat_type == "flow_50", ref_gage == "TRUE"), 
                      control = lmerControl(optimizer ="Nelder_Mead"))
prediction_models_ref$flow_50 <- predict(REF_FDC50_int)

# All Gages 25% exceedence (low probability of flow value)
ALL_FDC25_int <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                      data = filter(flow_diff_defol, stat_type == "flow_25"))
prediction_models_all$flow_25 <- predict(ALL_FDC25_int)
MuMIn::r.squaredGLMM(ALL_FDC25_int)
# Reference gages 25% exceedence (low probability of flow value)
REF_FDC25_int <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                      data = filter(flow_diff_defol, stat_type == "flow_25", ref_gage == "TRUE"))
prediction_models_ref$flow_25 <- predict(REF_FDC25_int)


# All Gages 75% exceedence (high probability of flow value)
ALL_FDC75_int <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                      data = filter(flow_diff_defol, stat_type == "flow_75"))
prediction_models_all$flow_75 <- predict(ALL_FDC75_int)
# Reference Gages 75% exceedence (high probability of flow value)
REF_FDC75_int <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                      data = filter(flow_diff_defol, stat_type == "flow_75", ref_gage == "TRUE"))
prediction_models_ref$flow_75 <- predict(REF_FDC75_int)

#05 Models Highest propability flow 


ALL_FDC05_int <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                      data = filter(flow_diff_defol, stat_type == "flow_05"))
prediction_models_all$flow_05 <- predict(ALL_FDC05_int)
# Reference gages
REF_FDC05_int <- lmer(flow_diff ~ defol_mean_pos + (1 | year), 
                      data = filter(flow_diff_defol, stat_type == "flow_05", ref_gage == "TRUE"))
prediction_models_ref$flow_05 <- predict(REF_FDC05_int)

#Stack Predictions for Graphing 
temp <- prediction_models_ref %>% 
  pivot_longer(cols = `flow_50`:`flow_05`, names_to = "stat_type", values_to = "prediction") %>% 
  mutate(reference_prediction = T) %>% 
  group_by(STAID, year)
prediction_models <- prediction_models_all %>% 
  pivot_longer(cols = `flow_50`:`flow_05`, names_to = "stat_type", values_to = "prediction") %>% 
  mutate(reference_prediction = F) %>% 
  group_by(STAID, year) %>% 
  full_join(temp)

flow_diff_defol <- flow_diff_defol %>% 
  group_by(STAID, year) %>% 
  left_join(prediction_models)
flow_type <- c("25", "50", "75")
all_plot <- list()
ref_plot <- list()
for (i in 1:length(flow_type)) {
  print(i)
  all_plot[[i]] <- ggplot(filter(flow_diff_defol, stat_type == paste0('flow_', flow_type[i]), 
                                                               year > 2015 & year < 2018, reference_prediction == "FALSE")) +
  geom_point(aes(x = defol_mean_pos, y = flow_diff)) +
  geom_line(aes( x = defol_mean_pos, y = prediction)) + 
  facet_wrap(~year) +
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = paste0("% change at ", flow_type[i], " % prob.")) +
         theme_cowplot()
       
  ref_plot[[i]] <- ggplot(filter(flow_diff_defol, stat_type == paste0('flow_', flow_type[i]), 
                                 year > 2015 & year < 2018,ref_gage == TRUE, reference_prediction == TRUE),
                          aes(x = defol_mean_pos, y = flow_diff)) +
         geom_point() +
         geom_line(aes( x = defol_mean_pos, y = prediction)) + 
         facet_wrap(~year) +
         ylim(c(-100,100)) +
         labs(x = "Defoliation metric", 
              y = paste0("% change at ", flow_type[i], " % prob.")) +
         theme_cowplot()
  } 

plot_grid(ref_plot[[1]], all_plot[[1]],
          ref_plot[[2]], all_plot[[2]],
          ref_plot[[3]], all_plot[[3]], nrow = 3)


# Plot changes in 5%tile against defoliation, by year
all_plot_05 <- ggplot(filter(flow_diff_defol, stat_type == "flow_05", 
                             reference_prediction == F)) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = prediction, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  #  facet_wrap(~year) +
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 5% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 25%tile against defoliation, by year
all_plot_25 <- ggplot(filter(flow_diff_defol, stat_type == "flow_25", 
                             reference_prediction == F)) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = prediction, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  #  facet_wrap(~year) +
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 25% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 50%tile against defoliation, by year
all_plot_50 <- ggplot(filter(flow_diff_defol, stat_type == "flow_50", 
                             reference_prediction == F))  + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = prediction, color = year), linetype = "solid")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  #  facet_wrap(~year) +
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 50% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")


# Plot changes in 75%tile against defoliation, by year
all_plot_75 <- ggplot(filter(flow_diff_defol, stat_type == "flow_75", 
                             reference_prediction == F)) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = prediction, color = year), linetype = "dashed")+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  #  facet_wrap(~year) +
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 75% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")

ref_pertiles_leg <- get_legend(all_plot_75 + theme(legend.position="bottom"))
top_row <- cowplot::plot_grid(all_plot_25, all_plot_50, all_plot_75, labels = c("A", "B", "C"), nrow = 1)
FDC_change_plots <- cowplot::plot_grid(top_row, ref_pertiles_leg, rel_heights = c(1,0.1),
                              nrow = 2)
plot_grid(FDC_change_plots, ref_pertiles_leg, ncol = 1, rel_heights = c(1, 0.2))

ref_plot_25 <- ggplot(filter(flow_diff_defol, stat_type == "flow_25", 
                             reference_prediction == T)) + 
  geom_point(aes( x = defol_mean_pos, y = flow_diff, color = year)) +
  geom_line(aes(x = defol_mean_pos, y = prediction, color = year))+
  scale_color_manual(values = c("chartreuse4", "darkgoldenrod1", "lightsalmon4", "mediumpurple4"))+
  #  facet_wrap(~year) +
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 25% prob."))) +
  theme_cowplot()+
  theme(legend.position="none")


plot_grid(ref_plot, all_plot, nrow = 1)

##Save FDC Plots 
# #all plot deviance 
# save_plot("Graphs/Main_Figures/FDC_all_PredictedValues_Int.jpg", FDC_change_plots, base_width = 12, base_height = 4)
# #individual years (facet wrppaed)
# save_plot("Graphs/Main_Figures/2016-2017_facetwrap_predictedValues_Int.jpg",
#           plot_grid(ref_plot[[1]], all_plot[[1]],
#                     ref_plot[[2]], all_plot[[2]],
#                     ref_plot[[3]], all_plot[[3]], nrow = 3),
#           base_width = 8, base_height = 10)

# FDC Curve Stat Graphs ---------------------------------------------------

FDC_75_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = round(c(coef(ALL_FDC75_int)$year[,1], 
                                               coef(REF_FDC75_int)$year[,1]),2), 
                           std_err_intercept = round(c(rep(sqrt(diag(vcov(ALL_FDC75_int)))[1],3), 
                                                       rep(sqrt(diag(vcov(REF_FDC75_int)))[1],3)),2),
                           slope = round(c(coef(ALL_FDC75_int)$year[,2],
                                         coef(REF_FDC75_int)$year[,2]),2), 
                           std_err_slope = round(c(rep(sqrt(diag(vcov(ALL_FDC75_int)))[2],3),
                                   rep(sqrt(diag(vcov(REF_FDC75_int)))[2],3)),2),
                           stat_type = rep("75tile",6))
FDC_50_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = round(c(coef(ALL_FDC50_int)$year[,1], 
                                               coef(REF_FDC50_int)$year[,1]),2), 
                           std_err_intercept = round(c(rep(sqrt(diag(vcov(ALL_FDC50_int)))[1],3), 
                                                       rep(sqrt(diag(vcov(REF_FDC50_int)))[1],3)),2),
                           slope = round(c(coef(ALL_FDC50_int)$year[,2],
                                           coef(REF_FDC50_int)$year[,2]),2), 
                           std_err_slope = round(c(rep(sqrt(diag(vcov(ALL_FDC50_int)))[2],3),
                                                   rep(sqrt(diag(vcov(REF_FDC50_int)))[2],3)),2),
                           stat_type = rep("50tile",6))
FDC_25_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = round(c(coef(ALL_FDC25_int)$year[,1], 
                                               coef(REF_FDC25_int)$year[,1]),2), 
                           std_err_intercept = round(c(rep(sqrt(diag(vcov(ALL_FDC25_int)))[1],3), 
                                                       rep(sqrt(diag(vcov(REF_FDC25_int)))[1],3)),2),
                           slope = round(c(coef(ALL_FDC25_int)$year[,2],
                                           coef(REF_FDC25_int)$year[,2]),2), 
                           std_err_slope = round(c(rep(sqrt(diag(vcov(ALL_FDC25_int)))[2],3),
                                                   rep(sqrt(diag(vcov(REF_FDC25_int)))[2],3)),2),
                           stat_type = rep("25tile",6))
FDC_lm_stats <- rbind(FDC_25_stats, FDC_50_stats) #drop 75 stats here because the model is not well supported
##Statistics Graphs moved and can be found in prior version of code 
#write.csv(FDC_lm_stats, file = "FDC_lm_stats_int.csv")
