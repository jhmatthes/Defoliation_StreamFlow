# 03_budgetModelsFigures.R
# This code requires that you run 02_calculateBudgetsFDC.R first. 
# This code esimates the stats models for defoliation ~ budget anomaly relationships
# and plots 
# 1. Fig 2: Plot defoliation by year in each watershed
# 2. Supp Fig 2: Plot growing season precipitation
# 3. Fig 3: Plot defoliation ~ runoff ratio / yield / precip anomalies for ref/non-ref gages
# 4. Supp Fig 1: Plot runoff ratio anomaly baseline vs 2015-2017 data
# 5. Model output table with estimated coefficients

library(ggplot2)

# FIG 2: Defoliation Boxplot x Year with Ref/Non-Ref gages marked
#png("figures/FIG2_v4.png", width = 4000, height = 2500, res = 600)
defol_boxplot <- ggplot(dischargePrecip, aes(as.factor(year), defol_mean, 
                                             color = ref_gage)) + 
  geom_boxplot(aes(color = ref_gage),
               outlier.colour = NA, show.legend=FALSE) +
  geom_hline(yintercept = 0, lty = 2) + 
  geom_point(aes(color = ref_gage, group = ref_gage), 
             alpha = 0.8, position = position_jitterdodge())+
  scale_color_manual(name = "Gage class", 
                     values = c("darkgrey","black"),
                     labels = c("Non-Reference", "Reference"))+
  labs(x = "Year", y = "Defoliation Metric")+
  cowplot::theme_cowplot()+
  theme(legend.key = element_rect(fill = NA, color = NA))
defol_boxplot
#dev.off()

# SUPP FIG 2: Yearly Precipitation Boxplot x Year with Ref/Non-Ref gages marked
#png("figures/SuppPrecip.png", width = 4000, height = 2500, res = 600)
precip_boxplot <- ggplot(dplyr::filter(dischargePrecip, year >=2015),
                         aes(year_fac, precip_ann*1000, color = ref_gage))+
  geom_hline(yintercept = mean(dischargePrecip_mean$precip_mean)*1000, lty = 2) + 
  geom_boxplot(aes(color = ref_gage),alpha = 0.2, 
               outlier.colour = NA, show.legend=FALSE)+
  geom_point(aes(color = ref_gage), size = 1.5, alpha = 0.8, 
             position = position_jitterdodge())+
  scale_color_manual(name = "Gage class", values = c("darkgrey", "black"),
                     labels = c("Non-Reference", "Reference"))+
  labs(x = "Year", y = "Precipitation (mm)")+
  cowplot::theme_cowplot()
precip_boxplot
#dev.off()

# FIG 3: A) Runoff Ratio ~ Defoliation, B) Yield and C) Precip
# Year random intercept
YPmod_all <- nlme::lme(runoff_ratio_anom ~ defol_mean, random= ~ 1 | year,
                       data = dischargePrecip, method="REML")

# Year random slope & intercept
YPmod_all2 <- nlme::lme(runoff_ratio_anom ~ defol_mean, random= ~ defol_mean | year,
                        data = dischargePrecip, method="REML")

# Test random intercept vs slope & intercept models
anova(YPmod_all, YPmod_all2)
dischargePrecip$YPmod_all<- predict(YPmod_all) # Random int is better 

RunoffRatio_All <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = runoff_ratio_anom, color = year_fac)) +
  geom_line(aes(x = defol_mean, y = YPmod_all, group = year_fac, color = year_fac), 
            size = 1, linetype = "solid") +
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
  labs(x = "Defoliation metric", y = "Runoff ratio anomaly", col = "Year")+
  ylim(c(-0.25, 0.5))+
  xlim(c(-1, 2))+
  cowplot::theme_cowplot()+
  theme(legend.position="none")

## FIG3b: Yield Anomaly ~ Defoliation and Precip Anomaly ~ Defoliation
Ymod_all <- nlme::lme(yield_anom ~ defol_mean, random= ~ 1 | year,
                data = dischargePrecip, method="REML")
dischargePrecip$Ymod_all<- predict(Ymod_all)

# Plot water yield anomalies ~ defoliation, w/ year random intercept
all_yieldanom <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = yield_anom, color = year_fac)) +
  geom_line(aes(x = defol_mean, y = Ymod_all, group = year_fac, color = year_fac), 
            size = 1, linetype = "solid") +
  scale_color_manual(name = "Year", values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
  labs(x = "Defoliation metric", y = "Water yield anomaly (m)")+
  cowplot::theme_cowplot()+
  theme(legend.position="none")

#Precip Anomaly Model
Precipmod_all <- nlme::lme(precip_anom ~ defol_mean, random= ~ 1 | year,
                     data = dischargePrecip, method="REML")
dischargePrecip$Precipmod_all<- predict(Precipmod_all) 

all_precipanom <- ggplot(dischargePrecip) +
  geom_point(aes(x = defol_mean, y = precip_anom, color = year_fac)) +
  geom_line(aes(x = defol_mean, y = Precipmod_all, group = year_fac, 
                color = year_fac), size = 1, linetype = "dashed") +
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
  labs(x = "Defoliation metric", y = "Precipitation anomaly (m)")+
  cowplot::theme_cowplot()+
  theme(legend.position="none")

RRYieldPrecip_leg <- cowplot::get_legend(all_yieldanom + theme(legend.position="bottom"))
RRYieldPrecip_plots <- cowplot::plot_grid(RunoffRatio_All, all_yieldanom, all_precipanom, labels = c("A", "B", "C"), nrow = 1)
#png("figures/Fig3_v3.png", width = 6000, height = 2500, res = 600)
cowplot::plot_grid(RRYieldPrecip_plots, RRYieldPrecip_leg, 
                   ncol = 1, rel_heights = c(1.1, 0.2))
#dev.off()

# REF GAGES ONLY: Runoff Ratio ~ defoliation
ref_dischargePrecip <- dplyr::filter(dischargePrecip, ref_gage == "Ref")
YPmod_ref <- nlme::lme(runoff_ratio_anom ~ defol_mean, random = ~1 | year,
                 data = ref_dischargePrecip, method = "REML")
ref_dischargePrecip$YPmod_ref <- predict(YPmod_ref) 

RunoffRatio_Ref <- ggplot(ref_dischargePrecip) +
  geom_point(aes(x = defol_mean, y = runoff_ratio_anom, color = year_fac)) +
  geom_line(aes(x = defol_mean, y = YPmod_ref, group = year_fac, color = year_fac), 
            size = 1, linetype = "solid") +
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
  labs(x = "Defoliation metric", y = "Runoff ratio anomaly", col = 'Year')+
  cowplot::theme_cowplot()+
  ylim(c(-0.25, 0.5)) +
  xlim(c(-1, 2)) +
  theme(legend.position="none")

## FIG3b: Yield Anomaly ~ Defoliation and Precip Anomaly ~ Defoliation
Ymod_ref <- nlme::lme(yield_anom ~ defol_mean, random = ~1 | year,
                data = ref_dischargePrecip, method = "REML")
ref_dischargePrecip$Ymod_ref<- predict(Ymod_ref)

# Plot water yield anomalies ~ defoliation, w/ year random intercept
ref_yieldanom <- ggplot(ref_dischargePrecip) +
  geom_point(aes(x = defol_mean, y = yield_anom, color = year_fac)) +
  geom_line(aes(x = defol_mean, y = Ymod_ref, group = year_fac, 
                color = year_fac), size = 1, linetype = "solid") +
  scale_color_manual(name = "Year", values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
  labs(x = "Defoliation metric", y = "Water yield anomaly (m)")+
  cowplot::theme_cowplot()+
  theme(legend.position="none")

#Precip Anomaly Model
Precipmod_ref <- nlme::lme(precip_anom ~ defol_mean, random = ~1 | year,
                     data = ref_dischargePrecip, method = "REML")
ref_dischargePrecip$Precipmod_ref <- predict(Precipmod_ref) 

ref_precipanom <- ggplot(ref_dischargePrecip) +
  geom_point(aes(x = defol_mean, y = precip_anom, color = year_fac)) +
  geom_line(aes(x = defol_mean, y = Precipmod_ref, group = year_fac, 
                color = year_fac), size = 1, linetype = "dashed") +
  scale_color_manual(values = c("#B8DE29FF","#20A387FF", "#440154FF"))+
  labs(x = "Defoliation metric", y = "Precipitation anomaly (m)")+
  cowplot::theme_cowplot()+
  theme(legend.position="none")

RRYieldPrecip_leg_ref <- cowplot::get_legend(ref_yieldanom + theme(legend.position="bottom"))
RRYieldPrecip_plots_ref <- cowplot::plot_grid(RunoffRatio_Ref, ref_yieldanom, ref_precipanom, labels = c("D", "E", "F"), nrow = 1)
#png("figures/Fig4_REFRRanoms.png", width = 6000, height = 2500, res = 600)
#cowplot::plot_grid(RRYieldPrecip_plots_ref, RRYieldPrecip_leg_ref, 
#                   ncol = 1, rel_heights = c(1.1, 0.2))
#dev.off()

#png("figures/Fig3_ALLanoms_300dpi.png", width = 3000, height = 2250, res = 300)
cowplot::plot_grid(RRYieldPrecip_plots,
                   RRYieldPrecip_plots_ref, RRYieldPrecip_leg_ref, 
                   ncol = 1, rel_heights = c(1.1, 1.1, 0.2))
#dev.off()

# SUPP FIG 1: Boxplot of RR anomalies and each year of defoliation
RunoffRatio_boxplot <- ggplot(data = dischargePrecip, 
                              aes(x = year_fac, y = runoff_ratio_anom, 
                                  color = defol_mean)) +
  geom_boxplot(data = dischargePrecip_baseline, aes(x = "1995-2014", y = runoff_ratio_anom),
               outlier.colour = NA, show.legend=FALSE, color = "black") +
  geom_boxplot(data = dischargePrecip, aes(x = year_fac, y = runoff_ratio_anom), 
               outlier.colour = NA, show.legend=FALSE, alpha = 0.8) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.8, position = "jitter")+
  ylim(-0.25, 0.4) +
  scale_color_gradient2(name = "Defoliation\nMetric", high = scales::muted("red"), 
                        mid = "lightgrey", low = scales::muted("blue")) +
  labs(x = "Year", y = "Runoff Ratio Anomaly")+
  cowplot::theme_cowplot()+
  theme(legend.key = element_rect(fill = NA, color = NA))

#png("figures/FigS_RRAnomalyDefol.png", width = 4000, height = 2500, res = 600)
RunoffRatio_boxplot
#dev.off()

# Bind together LMM models outputs in table format
RRAnoms_YReffects <- data.frame(data = c(rep("All", 3),rep("Ref Only", 3)),
                                response = rep("RR Anom",6),
                                year = rep(2015:2017,2),
                                intercept = signif(c(coef(YPmod_all)[,1],
                                                     coef(YPmod_ref)[,1]),2),
                                std_err_intercept = signif(c(rep(intervals(YPmod_all)$sigma[2],3),
                                                             rep(intervals(YPmod_ref)$sigma[2],3)),2),
                                slope = signif(c(coef(YPmod_all)[,2],
                                                 coef(YPmod_ref)[,2]),2), 
                                std_err_slope = signif(c(rep(sqrt(diag(vcov(YPmod_all)))[2],3), 
                                                         rep(sqrt(diag(vcov(YPmod_ref)))[2],3)),2))

YAnoms_YReffects <- data.frame(data = c(rep("All", 3),rep("Ref Only", 3)),
                               response = rep("Yield Anom",6),
                               year = rep(2015:2017,2),
                               intercept = signif(c(coef(Ymod_all)[,1],
                                                    coef(Ymod_ref)[,1]),2),
                               std_err_intercept = signif(c(rep(nlme::intervals(Ymod_all)$sigma[2],3),
                                                            rep(nlme::intervals(Ymod_ref)$sigma[2],3)),2),
                               slope = signif(c(coef(Ymod_all)[,2],
                                                coef(Ymod_ref)[,2]),2), 
                               std_err_slope = signif(c(rep(sqrt(diag(vcov(Ymod_all)))[2],3), 
                                                        rep(sqrt(diag(vcov(Ymod_ref)))[2],3)),2))

PAnoms_YReffects <- data.frame(data = c(rep("All", 3),rep("Ref Only", 3)),
                               response = rep("Precip Anom",6),
                               year = rep(2015:2017,2),
                               intercept = signif(c(coef(Precipmod_all)[,1],
                                                    coef(Precipmod_ref)[,1]),2),
                               std_err_intercept = signif(c(rep(nlme::intervals(Precipmod_all)$sigma[2],3),
                                                            rep(nlme::intervals(Precipmod_ref)$sigma[2],3)),2),
                               slope = signif(c(coef(Precipmod_all)[,2],
                                                coef(Precipmod_ref)[,2]),2), 
                               std_err_slope = signif(c(rep(sqrt(diag(vcov(Precipmod_all)))[2],3), 
                                                        rep(sqrt(diag(vcov(Precipmod_ref)))[2],3)),2))

ALLanoms_YReffects <- dplyr::bind_rows(RRAnoms_YReffects, YAnoms_YReffects, PAnoms_YReffects)

#write.csv(ALLanoms_YReffects, file = "figures/RRanoms_year_effects.csv", row.names = F)
