# 02_calculateBudgetsFDCs.R
# This code requires that you run 01_loadData.R first. 
# 1. Load and process 15-minute 1995-2017 stream gage data.
# 1a. Calculate monthly discharge sums 
# 1b. Calculate baseline FDC (1995-2014) and departures in 2015-2017
# 2. Aggregate monthly discharge to growing season discharge
# 3. Calculate monthly precipitation inputs 1995-2017
# 4. Combine growing season metrics for stream and precipitation data
# 5. Calculate growing season 2015-2017 departures from baseline
# 6. Connect water budget data to defoliation data

# Source code for FDC calculations
source("R/Functions/FDC_calc.R") # Calculate FDC curves
source("R/Functions/FDC_percentiles.R") # Find percentiles on FDC curves

# Load and process 15-minute instantaneous stream gage data (previously downloaded locally):
# Aggregate data for monthly discharge & calculate statistics for flow duration curves
plot_FDC <- FALSE # Plot & save each gage's FDC curve in a PDF file? (set file path below)
start_month <- 6
end_month <- 9 
gage15min_path <- "data/Gage_Data/" 

for(g in 1:length(gages_STAID$STAID)){
  
  print(paste0("Working on ",g))
  
  # Read 15-min gage discharge data
  dat <- read.csv(paste0(gage15min_path,gages_STAID$STAID[g],"_15min.csv"), header=T) 
  dat <- tidyr::separate(dat, dateTime, into = c("date","time"), sep = "T", remove = F) 
  dat$time <- substr(dat$time,1,8)
  dat <- dplyr::mutate(dat, date = as.Date(date, format = "%Y-%m-%d"),
                       time = lubridate::hms(time),
                       month = lubridate::month(date), 
                       year = lubridate::year(date))
  
  # Calculate seasonal budgets of discharge for each stream gage x year
  watershed_area <- gages_STAID$DRAIN_SQKM[g]
  gages_monthlyDischarge_hires_site <- dat %>% 
    dplyr::filter(month >= start_month & month <= end_month) %>%
    dplyr::mutate(month_year = paste0(month,"-",year)) %>%
    dplyr::group_by(site_no, month_year, month, year) %>%
    dplyr::summarize(discharge_monthly = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE), #ft3/s to m3/s to m3/15min per subwatershed
                     yield_monthly = sum(X_00060_00000*0.0283168*60*15, na.rm=TRUE)/(mean(watershed_area)*1000^2), # yield in m (m3/m2) 
                     sum_NA = sum(is.na(X_00060_00000)),
                     sum_good = length(grep("A", X_00060_00000_cd))) %>%
    dplyr::mutate(days_good = sum_good/(4*24),
                  days_good_frac = sum_good/(4*24)/lubridate::days_in_month(month))
  
  # Calculate baseline (1995-2014) FDC and pull percentile flow stats
  FDC_colnames <- c("STAID", "datasub", "FDC_type", "flow_05", "flow_25", 
                    "flow_50", "flow_75", "flow_95")
  dat_baseline <- dplyr::filter(dat, date < "2015-01-01")
  baseline <- FDC_calc(dat_baseline$X_00060_00000*0.0283168) #ft3/s to m3/s
  baseline_stats <- FDC_pertiles(baseline, gages_STAID$STAID[g], "baseline")
  
  # Add 2015, 2016, 2017, 2018 flow duration curves & pull percentile flow stats
  dat_2015 <- dplyr::filter(dat, date >= paste0("2015-0",start_month,"-01") & 
                              date <= paste0("2015-0",end_month,"-30"))
  if(nrow(dat_2015) != 0){
    FDC_2015 <- FDC_calc(dat_2015$X_00060_00000*0.0283168) #ft3/s to m3/s
    F2015_stats <- FDC_pertiles(FDC_2015, gages_STAID$STAID[g], "2015")
  } else {
    F2015_stats_vals <- c(gages_STAID$STAID[g], 2015, rep(NA, 6))
    F2015_stats <- rbind(F2015_stats_vals, F2015_stats_vals)
    colnames(F2015_stats) <- FDC_colnames
  }
  
  dat_2016 <- dplyr::filter(dat, date >= paste0("2016-0",start_month,"-01") & 
                              date <= paste0("2016-0",end_month,"-30"))
  if(nrow(dat_2016) != 0){
    FDC_2016 <- FDC_calc(dat_2016$X_00060_00000*0.0283168) #ft3/s to m3/s
    F2016_stats <- FDC_pertiles(FDC_2016, gages_STAID$STAID[g], "2016")
  } else {
    F2016_stats_vals <- c(gages_STAID$STAID[g], 2016, rep(NA, 6))
    F2016_stats <- rbind(F2016_stats_vals, F2016_stats_vals)
    colnames(F2016_stats) <- FDC_colnames
  }
  
  dat_2017 <- dplyr::filter(dat, date >= paste0("2017-0",start_month,"-01") & 
                              date <= paste0("2017-0",end_month,"-30"))
  if(nrow(dat_2017) != 0){
    FDC_2017 <- FDC_calc(dat_2017$X_00060_00000*0.0283168) #ft3/s to m3/s
    F2017_stats <- FDC_pertiles(FDC_2017, gages_STAID$STAID[g], "2017")
  } else {
    F2017_stats_vals <- c(gages_STAID$STAID[g], 2017, rep(NA, 6))
    F2017_stats <- rbind(F2017_stats_vals, F2017_stats_vals)
    colnames(F2017_stats) <- FDC_colnames
  }
  
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
    pdf(paste0("FDC_plots/FDC_plot_",dat$site_no[g],".pdf"))
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

# Make months with more than 90% missing monthly discharge data NA
gages_monthlyDischarge_hires$days_good_frac[gages_monthlyDischarge_hires$days_good_frac < 0.90] <- NA

# Summarize monthly to seasonal yearly discharge, only keep years with > 95% non-NA data
gages_seasonalDischarge <- gages_monthlyDischarge_hires %>%
  dplyr::mutate(STAID = stringr::str_pad(as.factor(site_no),side="left",pad="0",width=8)) %>%
  dplyr::group_by(STAID, year)  %>% 
  dplyr::summarize(discharge_seasonal = sum(discharge_monthly), #m3/subwatershed
                   yield_seasonal = sum(yield_monthly), #m at gage
                   sum_NA = sum(is.na(discharge_monthly)),
                   prop_good = sum(days_good)/(30+31+31+30)) %>%
  dplyr::filter(prop_good > 0.95) 

# Make table of coverage stats for min/max year and number of years per gage
gages_coverage <- gages_seasonalDischarge %>%
  dplyr::group_by(STAID) %>%
  dplyr::summarize(yr_start = min(year),
                   yr_end = max(year),
                   duration = max(year)-min(year) +1,
                   sum_NA = sum(is.na(discharge_seasonal)),
                   prop_yrs_missing = sum_NA/duration) 

# Remove gages with fewer than 10 years of baseline data
gages_seasonalDischarge <- dplyr::left_join(gages_seasonalDischarge, gages_coverage)  %>%
  dplyr::ungroup() %>% 
  dplyr::filter(duration > 14) %>%
  dplyr::mutate(STAID = stringr::str_pad(STAID, width = 8, side = "left", pad = "0"))

# Calculate seasonal Daymet precipitation at each stream gage
for(g in 1:nrow(gages_STAID)){
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
  
  # Aggregate seasonal gage precip into one dataframe
  if(g == 1){
    gages_seasonalPrecip <- gage_seasonalPrecip
  } else {
    gages_seasonalPrecip <- rbind(gages_seasonalPrecip, gage_seasonalPrecip)
  }
}

# Combine precipitation and discharge data: calculate growing season yield-to-precip ratio
gages_dischargePrecip <- dplyr::full_join(gages_seasonalPrecip, gages_seasonalDischarge, 
                                          by = c("STAID", "year")) %>%
  dplyr::mutate(runoff_ratio = yield_seasonal/precip_ann, na.rm=TRUE) 

# Remove gages that had bad FDC curve data (controlled flow)
bad_gages <- c("01184100", "01186000", "01202501", "01205500", "01206900")
gages_dischargePrecip <- gages_dischargePrecip[!(gages_dischargePrecip$STAID %in% bad_gages),]

# Calculate the departures in precipitation & discharge from the 20-year mean (1995-2014)
# Calculate the standard deviation in the baseline period for precip & discharge
dischargePrecip_mean <- gages_dischargePrecip %>%
  dplyr::filter(year < 2015) %>%
  dplyr::group_by(STAID) %>%
  dplyr::summarize(precip_mean = mean(precip_ann, na.rm=TRUE),
                   precip_sd = sd(precip_ann, na.rm=TRUE),
                   precip_gage_mean = mean(precip_gage, na.rm=TRUE),
                   precip_gage_sd = sd(precip_gage, na.rm=TRUE),
                   discharge_mean = mean(discharge_seasonal, na.rm=TRUE),
                   discharge_sd = sd(discharge_seasonal, na.rm=TRUE),
                   yield_mean = mean(yield_seasonal, na.rm=TRUE),
                   yield_sd = sd(yield_seasonal, na.rm=TRUE),
                   runoff_ratio_mean = mean(runoff_ratio, na.rm=TRUE),
                   runoff_ratio_sd = sd(runoff_ratio, na.rm=TRUE)) %>%
  dplyr::filter(runoff_ratio_mean < 2)

# Calculate baseline period departures in runoff ratio, yield, precipitation
dischargePrecip_baseline <- gages_dischargePrecip %>%
  dplyr::filter(year < 2015) %>%
  dplyr::left_join(dischargePrecip_mean) %>%
  dplyr::mutate(runoff_ratio_anom = runoff_ratio - runoff_ratio_mean,
                yield_anom = yield_seasonal-yield_mean,
                precip_anom = precip_ann - precip_mean)

# Calculate the 2015-2017 growing season departures from baseline mean
dischargePrecip <- dplyr::left_join(dischargePrecip_mean, gages_dischargePrecip, by = "STAID") %>%
  dplyr::mutate(precip_anom = precip_ann - precip_mean,
                discharge_anom = discharge_seasonal - discharge_mean,
                yield_anom = yield_seasonal - yield_mean,
                runoff_ratio_anom = runoff_ratio - runoff_ratio_mean)

# Connect stream data to defoliation data                   
gages_defol_vals <- gages_defol %>% 
  dplyr::filter(STAID %in% unique(dischargePrecip$STAID)) %>% 
  dplyr::mutate(defol_mean = defol_mean*-1)

dischargePrecip <- dischargePrecip %>% 
  dplyr::left_join(gages_defol_vals, by = c("STAID", "year")) %>% #make defoliation positive 
  dplyr::filter(STAID != '01105880', year >=2015 & year<2018,
                !is.na(runoff_ratio_anom))  #remove stream gage on cape cod, that is not a single basin

# Format factors
dischargePrecip$ref_gage <- as.factor(dischargePrecip$ref_gage)
levels(dischargePrecip$ref_gage) <- c("Non-Ref", "Ref")
dischargePrecip$year_fac <- as.factor(dischargePrecip$year)

# next up: 03_budgetModelsFigures.R to estimate relationships between defoliation and 
# growing season runoff ratio, precipitation, yield budgets and make figures.
# OR: 04_fdcModelsFigures.R to estimate FDC departure relationships iwth defoliation. 
