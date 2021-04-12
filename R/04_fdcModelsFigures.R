# 04_fdcModelsFigures.R
# This code requires that you run 02_calculateBudgetsFDC.R first. 
# 1. Calculate the FDC percentile departures in 2015-2017 from baseline
# 2. Estimate FDC percentile departures vs defoliation models
# 3. Figure 4: Defoliation vs FDC departures for reference and non-reference gages
# 4. Model output table with FDC stats

# FLOW DURATION CURVE STATISTICS ------------------------------------------
# Clean up FDC stats by site
FDC_stats_sites <- dplyr::mutate(FDC_stats_sites, year = datasub) %>%
  dplyr::filter(STAID %in% unique(dischargePrecip$STAID))

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
    if(nrow(dplyr::filter(discharge_site, year == "2015")) != 0) {
      diff_15 <- (as.numeric(dplyr::filter(discharge_site, year == "2015")[,stats_columns[s]])-
                    baseline)/baseline
    } else {diff_15 <- NA}
    
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
FDC_data <- tidyr::pivot_wider(flow_diff_defol, names_from = stat_type, values_from = flow_diff) %>%
  dplyr::filter(!is.na(flow_50))

# All Gages 50% exceedence (medium probability of flow value) Model
ALL_FDC50_int <- nlme::lme(flow_50 ~ defol_mean, random = ~1 | year, 
                           data = FDC_data, method = "REML")

ALL_FDC50_sl <- nlme::lme(flow_50 ~ defol_mean, random = ~defol_mean | year, 
                          data = FDC_data, method = "REML")
anova(ALL_FDC50_int, ALL_FDC50_sl)

FDC_data$flow_50_model <- as.vector(predict(ALL_FDC50_sl)) # slope model is better

# All Gages 50% exceedence (medium probability of flow value) Model
FDC_REFdata <- dplyr::filter(FDC_data, ref_gage == "Ref")
REF_FDC50_int <- nlme::lme(flow_50 ~ defol_mean, random = ~1 | year, 
                           data = FDC_REFdata, method = "REML")
REF_FDC50_sl <- nlme::lme(flow_50 ~ defol_mean, random = ~defol_mean | year, 
                          data = FDC_REFdata, method = "REML")
anova(REF_FDC50_int, REF_FDC50_sl)
nlme::ranef(REF_FDC50_sl)
FDC_REFdata$flow_50_model <- predict(REF_FDC50_sl)

# All Gages 25% exceedence (low probability of flow value) Model
ALL_FDC25_int <- nlme::lme(flow_25 ~ defol_mean, random = ~1 | year, 
                           data = FDC_data, method = "ML")
ALL_FDC25_sl <- nlme::lme(flow_25 ~ defol_mean, random = ~defol_mean | year, 
                          data = FDC_data, method = "ML")
anova(ALL_FDC25_int, ALL_FDC25_sl)

FDC_data$flow_25_model <- predict(ALL_FDC25_sl)

# Reference gages 25% exceedence (low probability of flow value) Model
REF_FDC25_int <- nlme::lme(flow_25 ~ defol_mean, random = ~1 | year, 
                           data = FDC_REFdata, method = "REML")
REF_FDC25_sl <- nlme::lme(flow_25 ~ defol_mean, random = ~defol_mean | year, 
                          data = FDC_REFdata, method = "REML")
anova(REF_FDC25_int, REF_FDC25_sl)
FDC_REFdata$flow_25_model <- predict(REF_FDC25_sl)

# All Gages 75% exceedence (high probability of flow value) Model
ALL_FDC75_int <- nlme::lme(flow_75 ~ defol_mean, random = ~1 | year, 
                           data = FDC_data, method = "REML")

ALL_FDC75_sl <- nlme::lme(flow_75 ~ defol_mean, random = ~defol_mean | year, 
                          data = FDC_data, method = "REML")
anova(ALL_FDC75_int, ALL_FDC75_sl)

FDC_data$flow_75_model <- predict(ALL_FDC75_sl)

# Reference Gages 75% exceedence (high probability of flow value) Model
REF_FDC75_int <- nlme::lme(flow_75 ~ defol_mean, random = ~1 | year, 
                           data = FDC_REFdata, method = "ML")

REF_FDC75_sl <- nlme::lme(flow_75 ~ defol_mean, random = ~defol_mean | year, 
                          data = FDC_REFdata, method = "REML")

FDC_REFdata$flow_75_model <- predict(REF_FDC75_sl)

#### Plot defoliation ~ FDC percentile changes
#FDC_data <- mutate(FDC_data, Year = year)
# Plot changes in 25%tile against defoliation, by year
all_plot_25 <- ggplot(FDC_data) + 
  geom_point(aes( x = defol_mean, y = flow_25, color = year)) +
  geom_line(aes(x = defol_mean, y = flow_25_model, color = year), linetype = "solid")+
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
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
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 50% prob."))) +
  cowplot::theme_cowplot()+
  theme(legend.position="none")

# Plot changes in 75%tile against defoliation, by year
all_plot_75 <- ggplot(FDC_data) + 
  geom_point(aes( x = defol_mean, y = flow_75, color = year)) +
  geom_line(aes(x = defol_mean, y = flow_75_model, color = year), linetype = "solid")+
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 75% prob."))) +
  cowplot::theme_cowplot()+
  theme(legend.position="none")

all_FDC_leg <- cowplot::get_legend(all_plot_75 + theme(legend.position="bottom"))
FDC_all <- cowplot::plot_grid(all_plot_25, all_plot_50, all_plot_75, labels = c("A", "B", "C"), nrow = 1)

#png("figures/FIG4.png", width = 6000, height = 2500, res = 600)
cowplot::plot_grid(FDC_all, all_FDC_leg, ncol = 1, rel_heights = c(1.1, 0.2))
#dev.off()

# REFERENCE PLOTS ONLY
ref_plot_25 <- ggplot(FDC_REFdata) + 
  geom_point(aes( x = defol_mean, y = flow_25, color = year)) +
  geom_line(aes(x = defol_mean, y = flow_25_model, color = year))+
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
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
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
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
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
  ylim(c(-100,100)) +
  labs(x = "Defoliation metric", 
       y = expression(paste("Percent ", Delta, 
                            " flow at 75% prob."))) +
  cowplot::theme_cowplot()+
  theme(legend.position="none")

ref_FDC_leg <- cowplot::get_legend(ref_plot_75 + theme(legend.position="bottom"))
FDC_ref <- cowplot::plot_grid(ref_plot_25, ref_plot_50, ref_plot_75, labels = c("D", "E", "F"), nrow = 1)

#png("figures/Fig4.png", width = 6000, height = 4500, res = 600)
#cowplot::plot_grid(FDC_all, FDC_ref, top_row, ref_FDC_leg, ncol = 1, rel_heights = c(1.1, 1.1, 0.2))
cowplot::plot_grid(FDC_all, FDC_ref, ref_FDC_leg, ncol = 1, rel_heights = c(1.1, 1.1, 0.2))
#dev.off()

# FDC Statistics Tables ---------------------------------------------------
#Bind LMM outputs in data.frame 
FDC_75_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = signif(c(coef(ALL_FDC75_int)[,1],
                                                coef(REF_FDC75_int)[,1]),2),
                           std_err_intercept = signif(c(rep(sqrt(diag(vcov(ALL_FDC75_int)))[1],3),
                                                        rep(sqrt(diag(vcov(ALL_FDC75_int)))[1],3)),2),
                           slope = signif(c(coef(ALL_FDC75_int)[,2],
                                            coef(REF_FDC75_int)[,2]),2), 
                           std_err_slope = signif(c(rep(sqrt(diag(vcov(ALL_FDC75_int)))[2],3), 
                                                    rep(sqrt(diag(vcov(REF_FDC75_int)))[2],3)),2),
                           stat_type = rep("75tile",6),
                           model_type = rep("int",6))

FDC_50_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = signif(c(coef(ALL_FDC50_sl)[,1], 
                                                coef(REF_FDC50_sl)[,1]),2), 
                           std_err_intercept = signif(c(rep(sqrt(diag(vcov(ALL_FDC50_sl)))[1],3),
                                                        rep(sqrt(diag(vcov(REF_FDC50_sl)))[1],3)),2),
                           slope = signif(c(coef(ALL_FDC50_sl)[,2],
                                            coef(REF_FDC50_sl)[,2]),2), 
                           std_err_slope = signif(c(rep(sqrt(diag(vcov(ALL_FDC50_sl)))[2],3),
                                                    rep(sqrt(diag(vcov(REF_FDC50_sl)))[2],3)),2),
                           stat_type = rep("50tile",6),
                           model_type = rep("slope & int",6))
FDC_25_stats <- data.frame(data = rep(c(rep("All",3),rep("Ref Only",3))),
                           year = rep(2015:2017,2),
                           intercept = signif(c(coef(ALL_FDC25_sl)[,1], 
                                                coef(REF_FDC25_sl)[,1]),2), 
                           std_err_intercept = signif(c(rep(sqrt(diag(vcov(ALL_FDC25_sl)))[1],3),
                                                        rep(sqrt(diag(vcov(REF_FDC25_sl)))[1],3)),2),
                           slope = signif(c(coef(ALL_FDC25_sl)[,2],
                                            coef(REF_FDC25_sl)[,2]),2), 
                           std_err_slope = signif(c(rep(sqrt(diag(vcov(ALL_FDC25_sl)))[2],3),
                                                    rep(sqrt(diag(vcov(REF_FDC25_sl)))[2],3)),2),
                           stat_type = rep("25tile",6),
                           model_type = rep("slope & int",6))
FDC_lm_stats <- dplyr::bind_rows(FDC_25_stats, FDC_50_stats, FDC_75_stats) 
#write.csv(FDC_lm_stats, file = "figures/FDC_lm_stats_int_v2.csv", row.names=F)
