# Assess the percent departure from the baseline FDC curve for each 
# of the 20 years in the baseline period.
# Need FDC_calc() and FDC_pertiles() functions loaded from R/ folder.

FDC_baseline_variability <- function(start_month, end_month, 
                                     gages_STAID, gage15min_path){
  
  # Remove gages determined to be bad (dammed, flashy) from baseline FDC analysis
  bad_gages <- c("01184100", "01186000", "01202501", "01205500", "01206900")
  gages_STAID <- gages_STAID[!(gages_STAID$STAID %in% bad_gages),]
  
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
    
    # Calculate baseline (pre-2015) flow duration curve & pull percentile flow stats 
    dat_baseline <- dplyr::filter(dat, date < "2015-01-01")
    baseline <- FDC_calc(dat_baseline$X_00060_00000*0.0283168) #ft3/s to m3/s
    baseline_stats <- FDC_pertiles(baseline, gages_STAID$STAID[g], "baseline")
    baseline_discharge <- dplyr::filter(baseline_stats, FDC_type == "discharge")
    
    # # Print baseline FDCs for visual record
    # if(print_FDC == TRUE){
    #   pdf(paste0("figures/Baseline_FDCs/",gages_STAID$STAID[g],".pdf"))
    #   p1 <- ggplot2::ggplot(baseline) +
    #     ggplot2::geom_line(ggplot2::aes(x = exceed, y = q_dis)) +
    #     ggplot2::geom_hline(yintercept = median(baseline$q_dis), lty=2) +
    #     ggplot2::geom_vline(xintercept = 0.5, lty=2) +
    #     ggplot2::scale_y_log10() +
    #     ggplot2::labs(y = expression(paste("Discharge (",m^{3},s^{-1},")")),
    #                   x = "Probability of exceeding discharge rate",
    #                   title = paste0("Baseline at gage: ",gages_STAID$STAID[g]))
    #   print(p1)
    #   dev.off()
    # }
  
    # Calculate FDC percentiles for each year of baseline data
    baseline_years <- unique(dat_baseline$year)
    
    for(y in 1:length(baseline_years)){
      dat_year <- dplyr::filter(dat, 
                                date >= paste0(baseline_years[y],"-0",start_month,"-01") &
                                  date <= paste0(baseline_years[y],"-0",end_month,"-30"))
      
      if(nrow(dat_year) != 0){
        FDC_year <- FDC_calc(dat_year$X_00060_00000*0.0283168) #ft3/s to m3/s
        Fyear_stats <- FDC_pertiles(FDC_year, gages_STAID$STAID[g], paste(baseline_years[y]))
        
        # Calculate FDC defoliation year departure from baseline stats
        FDC_discharge <- dplyr::filter(Fyear_stats, FDC_type == "discharge")
        
        # Get stats column indices for difference-from-baseline calculation
        stats_columns <- grep("flow_", colnames(FDC_discharge))
        
        for(s in 1:length(stats_columns)){
          
          # Difference between year and baseline FDC discharge
          if(nrow(FDC_discharge) != 0) {
            year_diff <- as.numeric(FDC_discharge[,stats_columns[s]]-
                                      baseline_discharge[,stats_columns[s]])/baseline_discharge[,stats_columns[s]]
          } else {
            year_diff <- NA
          }
          
          # Build data frame with sites, stats, years
          tmp <- data.frame(STAID = as.character(gages_STAID$STAID[g]), 
                            stat_type = colnames(FDC_discharge)[stats_columns[s]],
                            year = baseline_years[y],
                            flow_diff = year_diff)
          
          if(g == 1 & y == 1 & s == 1){
            FDC_baselinedepartures <- tmp
          } else {
            FDC_baselinedepartures <- rbind(FDC_baselinedepartures, tmp)
          }
        } # end stat cols
      } else { # if year & station are missing data
        tmp <- data.frame(STAID = rep(as.character(gages_STAID$STAID[g]), 5), 
                          stat_type = colnames(baseline_discharge)[4:8],
                          year = rep(baseline_years[y], 5),
                          flow_diff = rep(NA, 5))
        if(g == 1 & y == 1){
          FDC_baselinedepartures <- tmp
        } else {
          FDC_baselinedepartures <- rbind(FDC_baselinedepartures, tmp)
        }
      } 
    } # end years
  }# end gages
}

# Summarize & plot regional FDC baseline departures to 
# as median +/- SD at percentiles by year (averaged across STAID) 
FDC_baselinesummary_year <- FDC_baselinedepartures %>%
  dplyr::filter(stat_type %in% c("flow_25","flow_50","flow_75")) %>%
  dplyr::group_by(year, stat_type) %>%
  dplyr::summarize(dep_median = median(flow_diff, na.rm=TRUE),
                   dep_mean = mean(flow_diff, na.rm=TRUE),
                   dep_sd = sd(flow_diff, na.rm=TRUE),
                   dep_na = sum(is.na(flow_diff)),
                   dep_count = n())

FDC_defolsummary_year <- flow_diff_defol %>%
  dplyr::filter(stat_type %in% c("flow_25","flow_50","flow_75")) %>%
  dplyr::group_by(year, stat_type) %>%
  dplyr::summarize(dep_median = median(flow_diff, na.rm=TRUE),
                   dep_mean = mean(flow_diff, na.rm=TRUE),
                   dep_sd = sd(flow_diff, na.rm=TRUE),
                   dep_na = sum(is.na(flow_diff)),
                   dep_count = n())
  
f25_baselinesd <- ggplot(filter(FDC_baselinesummary_year, stat_type == "flow_25")) +
  geom_hline(aes(yintercept=0), lty = 2) +
  geom_point(aes(x = year, y = dep_mean*100)) +
  geom_errorbar(aes(x = year, ymax = (dep_mean + dep_sd)*100, 
                    ymin = (dep_mean - dep_sd)*100)) +
  geom_point(data = filter(FDC_defolsummary_year, stat_type == "flow_25"), 
                           aes(x = as.numeric(year), 
                               y = dep_mean),col="red") +
  geom_errorbar(data = filter(FDC_defolsummary_year, stat_type == "flow_25"),
                aes(x = as.numeric(year), ymax = (dep_mean + dep_sd), 
                    ymin = (dep_mean - dep_sd)), col = "red") +
  labs(x = "Year", y = "", title = "25% Exceedance") +
  cowplot::theme_cowplot() 

f50_baselinesd <- ggplot(filter(FDC_baselinesummary_year, stat_type == "flow_50")) +
  geom_hline(aes(yintercept=0), lty = 2) +
  geom_point(aes(x = year, y = dep_median*100)) +
  geom_errorbar(aes(x = year, ymax = (dep_median + dep_sd)*100, 
                    ymin = (dep_median - dep_sd)*100)) +
  geom_point(data = filter(FDC_defolsummary_year, stat_type == "flow_50"), 
             aes(x = as.numeric(year), 
                 y = dep_median),col="red") +
  geom_errorbar(data = filter(FDC_defolsummary_year, stat_type == "flow_50"),
                aes(x = as.numeric(year), ymax = (dep_median + dep_sd), 
                    ymin = (dep_median - dep_sd)), col = "red") +
  labs(x = "Year", y = "", title = "50% Exceedance") +
  cowplot::theme_cowplot() 

f75_baselinesd <- ggplot(filter(FDC_baselinesummary_year, stat_type == "flow_75")) +
  geom_hline(aes(yintercept=0), lty = 2) +
  geom_point(aes(x = year, y = dep_median*100)) +
  geom_errorbar(aes(x = year, ymax = (dep_median + dep_sd)*100, 
                    ymin = (dep_median - dep_sd)*100)) +
  geom_point(data = filter(FDC_defolsummary_year, stat_type == "flow_75"), 
             aes(x = as.numeric(year), 
                 y = dep_median),col="red") +
  geom_errorbar(data = filter(FDC_defolsummary_year, stat_type == "flow_75"),
                aes(x = as.numeric(year), ymax = (dep_median + dep_sd), 
                    ymin = (dep_median - dep_sd)), col = "red") +
  labs(x = "Year", y = "", title = "75% Exceedance") +
  cowplot::theme_cowplot() 

plot_column <- cowplot::plot_grid(f25_baselinesd, f50_baselinesd, 
                              f75_baselinesd, labels = c("A", "B", "C"), nrow = 1)

y.grob <- grid::textGrob("Regional % Baseline Flow Anomaly", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)

#png("figures/FIGS_RegionalFDCBaselineDep.png", width = 6000, height = 2500, res = 600)
gridExtra::grid.arrange(arrangeGrob(plot_column, left = y.grob))
#dev.off()

# Calculate regional SD in flow anomalies in baseline vs 2015-2017
FDC_baselinesummary_allyears <- FDC_baselinedepartures %>%
  dplyr::filter(stat_type %in% c("flow_25","flow_50","flow_75")) %>%
  dplyr::group_by(stat_type) %>%
  dplyr::summarize(dep_median = median(flow_diff, na.rm=TRUE),
                   dep_mean = mean(flow_diff, na.rm=TRUE),
                   dep_sd = sd(flow_diff, na.rm=TRUE),
                   dep_na = sum(is.na(flow_diff)),
                   dep_count = n())

FDC_defolsummary_year <- flow_diff_defol %>%
  dplyr::filter(stat_type %in% c("flow_25","flow_50","flow_75")) %>%
  dplyr::group_by(year, stat_type) %>%
  dplyr::summarize(dep_median = median(flow_diff/100, na.rm=TRUE),
                   dep_mean = mean(flow_diff/100, na.rm=TRUE),
                   dep_sd = sd(flow_diff/100, na.rm=TRUE),
                   dep_na = sum(is.na(flow_diff/100)),
                   dep_count = n()) %>%
  dplyr::arrange(stat_type)

#write.csv(FDC_defolsummary_year, "figures/FDC_defolsummary_year.csv")  
#write.csv(FDC_baselinesummary_year, "figures/FDC_baselinesummary_year.csv")  
#write.csv(FDC_baselinedefol_staid, "figures/FDC_baselinedefol_staid.csv")

# Calculate interannual variation in FDC baseline departures 
# by SD at percentiles x STAID (across year)
FDC_baselinesummary_STAID <- FDC_baselinedepartures %>%
  dplyr::filter(stat_type %in% c("flow_25","flow_50","flow_75")) %>%
  dplyr::group_by(STAID, stat_type) %>%
  dplyr::summarize(dep_median = round(median(flow_diff*100, na.rm=TRUE),2),
                   dep_mean = round(mean(flow_diff*100, na.rm=TRUE),2),
                   dep_sd = round(sd(flow_diff*100, na.rm=TRUE),2),
                   dep_na = sum(is.na(flow_diff)),
                   dep_count = n()) %>%
  dplyr::select(STAID, stat_type, dep_sd) %>%
  tidyr::pivot_wider(names_from = stat_type, values_from = dep_sd) %>%
  dplyr::mutate(data = "1995-2014")

FDC_defolsummary_STAID <- flow_diff_defol %>%
  dplyr::filter(stat_type %in% c("flow_25","flow_50","flow_75")) %>%
  dplyr::select(STAID, year, stat_type, flow_diff) %>%
  dplyr::group_by(STAID, stat_type) %>%
  dplyr::summarize(dep_median = median(flow_diff, na.rm=TRUE),
                   dep_mean = mean(flow_diff, na.rm=TRUE),
                   dep_sd = sd(flow_diff, na.rm=TRUE),
                   dep_na = sum(is.na(flow_diff)),
                   dep_count = n()) %>%
  dplyr::select(STAID, stat_type, dep_sd) %>%
  tidyr::pivot_wider(names_from = stat_type, values_from = dep_sd) %>%
  dplyr::mutate(data = "2015-2017")

FDC_baselinedefol_staid <- bind_rows(FDC_baselinesummary_STAID, 
                                     FDC_defolsummary_STAID) %>%
  tidyr::pivot_longer(flow_25:flow_75, names_to="percentile",
                      values_to="SD_departure") %>%
  dplyr::arrange(STAID) 

# New facet label names
FDC_labs <- c("25% Exceedance", "50% Exceedance","75% Exceedance")
names(FDC_labs) <- c("flow_25", "flow_50", "flow_75")

#png("figures/FIGS_GageFDCBaselineDep.png", width = 5000, height = 3500, res = 600)
ggplot(FDC_baselinedefol_staid) +
  geom_point(aes(x = STAID, y = SD_departure, color = data)) +
  #labs(y = "σ FDC 25th percentile anomaly") +
  ylim(c(0,200)) +
  scale_color_manual(name = "Data", values = c("darkgrey","black")) +
  labs(x = "Stream gage (STAID)",y = "Std Dev of flow anomaly %") +
  facet_grid(~percentile, labeller = labeller(percentile = FDC_labs)) +
  cowplot::theme_cowplot() +
  theme(axis.text.x=element_blank()) 
#dev.off()


#### Code for checking gages with suspicious FDC curves
# For wonky FDC curves, investigate if there was a weird year or weird all the time: 
# STAIDs 0117550 01184100 01184490 01186000 01202501 01205500 01206900
# For gages with consistent anomalous flows, test: 01102000, 01101500, 01101000, 01095220
STAID <- "01101500"
g <- which(gages_STAID$STAID == STAID)

# Read 15-min gage discharge data
dat <- read.csv(paste0(gage15min_path,gages_STAID$STAID[g],"_15min.csv"), header=T) 
dat <- tidyr::separate(dat, dateTime, into = c("date","time"), sep = "T", remove = F) 
dat$time <- substr(dat$time,1,8)
dat <- dplyr::mutate(dat, date = as.Date(date, format = "%Y-%m-%d"),
                     time = lubridate::hms(time),
                     month = lubridate::month(date), 
                     year = lubridate::year(date))

# Calculate baseline (pre-2015) flow duration curve & pull percentile flow stats 
dat_baseline <- dplyr::filter(dat, date < "2015-01-01")
baseline <- FDC_calc(dat_baseline$X_00060_00000*0.0283168, baseline = T) #ft3/s to m3/s
baseline_stats <- FDC_pertiles(baseline, gages_STAID$STAID[g], "baseline")
baseline_discharge <- dplyr::filter(baseline_stats, FDC_type == "discharge")

# Calculate FDC percentiles for each year of baseline data
baseline_years <- unique(dat_baseline$year) # Diagnose number of years

for(y in 1:length(baseline_years)){
  dat_year <- dplyr::filter(dat, 
                            date >= paste0(baseline_years[y],"-0",start_month,"-01") &
                              date <= paste0(baseline_years[y],"-0",end_month,"-30"))
  
  if(nrow(dat_year) != 0){
    FDC_year <- FDC_calc(dat_year$X_00060_00000*0.0283168) #ft3/s to m3/s
    
    # Print baseline FDCs for visual record
    pdf(paste0("figures/FDCs_BehavingBadly/",gages_STAID$STAID[g],"_",
               baseline_years[y],".pdf"))
    p1 <- ggplot2::ggplot(baseline) +
      ggplot2::geom_line(ggplot2::aes(x = exceed, y = q_dis)) +
      ggplot2::geom_hline(yintercept = median(baseline$q_dis), lty=2) +
      ggplot2::geom_vline(xintercept = 0.5, lty=2) + 
      ggplot2::geom_line(data = FDC_year, aes(x = exceed, y = q_dis), color = "red") +
      ggplot2::scale_y_log10() +
      ggplot2::labs(y = expression(paste("Discharge (",m^{3},s^{-1},")")), 
                    x = "Probability of exceeding discharge rate",
                    title = paste0("Gage ",gages_STAID$STAID[g]," in year ", baseline_years[y]))
    print(p1)
    dev.off()
    
  }
}


# # Did not include in paper
# # Summarize FDC baseline departures to SD at percentiles x STAID (across year)
# FDC_baselinesummary_STAID <- FDC_baselinedepartures %>%
#   dplyr::filter(stat_type %in% c("flow_25","flow_50","flow_75")) %>%
#   dplyr::group_by(STAID, stat_type) %>%
#   dplyr::summarize(dep_median = median(flow_diff, na.rm=TRUE),
#                    dep_mean = mean(flow_diff, na.rm=TRUE),
#                    dep_sd = sd(flow_diff, na.rm=TRUE),
#                    dep_na = sum(is.na(flow_diff)),
#                    dep_count = n())
# 
# FDC_defolsummary_STAID <- flow_diff_defol %>%
#   dplyr::filter(stat_type %in% c("flow_25","flow_50","flow_75")) %>%
#   dplyr::group_by(STAID, stat_type) %>%
#   dplyr::summarize(dep_median = median(flow_diff, na.rm=TRUE),
#                    dep_mean = mean(flow_diff, na.rm=TRUE),
#                    dep_sd = sd(flow_diff, na.rm=TRUE),
#                    dep_na = sum(is.na(flow_diff)),
#                    dep_count = n())
# 
# f25_baselinesd_staid_base <- ggplot(filter(FDC_baselinesummary_STAID, 
#                                            stat_type == "flow_25")) +
#   geom_hline(aes(yintercept=0), lty = 2) +
#   geom_point(aes(x = STAID, y = dep_mean*100)) +
#   geom_errorbar(aes(x = STAID, ymax = (dep_mean + dep_sd)*100, 
#                     ymin = (dep_mean - dep_sd)*100)) +
#   labs(x = "Stream Gage (STAID)", y = "High flow", 
#        title = "Baseline: 1995-2014") +
#   ylim(c(-100,150)) +
#   cowplot::theme_cowplot() +
#   theme(axis.text.x=element_blank())
# 
# f25_baselinesd_staid_defol <- ggplot(filter(FDC_defolsummary_STAID, 
#                                             stat_type == "flow_25")) +
#   geom_hline(aes(yintercept=0), lty = 2) +
#   geom_point(aes(x = STAID, y = dep_mean)) +
#   geom_errorbar(aes(x = STAID, ymax = (dep_mean + dep_sd), 
#                     ymin = (dep_mean - dep_sd))) +
#   ylim(c(-100,150)) +
#   labs(x = "Stream Gage (STAID)", y = "High flow", 
#        title = "Defoliation: 2015-2017") +
#   cowplot::theme_cowplot() +
#   theme(axis.text.x=element_blank())
# 
