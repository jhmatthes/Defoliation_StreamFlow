# Pull flow at percentiles from FDC curve

FDC_pertiles <- function(FDC_df, site, dataname) {
  
  flow_05 <- FDC_df$q_dis[which.min(abs(FDC_df$exceed - 0.05))] 
  flow_25 <- FDC_df$q_dis[which.min(abs(FDC_df$exceed - 0.25))] 
  flow_50 <- FDC_df$q_dis[which.min(abs(FDC_df$exceed - 0.5))] 
  flow_75 <- FDC_df$q_dis[which.min(abs(FDC_df$exceed - 0.75))] 
  flow_95 <- FDC_df$q_dis[which.min(abs(FDC_df$exceed - 0.95))] 
  
  ptiles_df <- data.frame(STAID = as.character(site), datasub = dataname, FDC_type = "discharge",
                          flow_05, flow_25, flow_50, flow_75, flow_95)
  
  flow_05 <- FDC_df$q_per[which.min(abs(FDC_df$exceed - 0.05))] 
  flow_25 <- FDC_df$q_per[which.min(abs(FDC_df$exceed - 0.25))] 
  flow_50 <- FDC_df$q_per[which.min(abs(FDC_df$exceed - 0.5))] 
  flow_75 <- FDC_df$q_per[which.min(abs(FDC_df$exceed - 0.75))] 
  flow_95 <- FDC_df$q_per[which.min(abs(FDC_df$exceed - 0.95))] 
  
  ptiles_df_2  <- data.frame(STAID = as.character(site), datasub = dataname, FDC_type = "percent_mean",
                          flow_05, flow_25, flow_50, flow_75, flow_95)
  
  ptiles_df_3 <- rbind(ptiles_df, ptiles_df_2)
  
  return(ptiles_df_3)
}
  