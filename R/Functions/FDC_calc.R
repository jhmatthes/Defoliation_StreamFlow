#' Flow Duration Curve
#'
#' Produces a flow duration curve plot 
#' @param flow daily streamflow time series
#' @param baseline whether or not to do baseline calculations
#' @references Gustard, A., Bullock, A., and Dixon, J.M. (1992). Report No. 108:
#'   Low flow estimation in the United Kingdom. 
#'   Oxfordshire, United Kingdom: Institute of Hydrology.
#' @author Jackie Matthes

FDC_calc <- function(flow, baseline, mean_flow) {
  
  if(baseline == TRUE){
    # Sort the flow values and calculate the FDC statistics
    rank <- rank(flow, ties.method="max")
    rank <- max(rank) - rank
    exceedtime <- 1*(rank / (length(flow) + 1))
    q_per <- sort(100 * flow / mean(flow), decreasing=FALSE)
    q_dis <- sort(flow, decreasing = F)
    exceed <- sort(exceedtime, decreasing=TRUE)
    fdc_df <- data.frame(exceed, q_per, q_dis, mean_flow = mean(flow))
  } else {
    hist_mean <- mean_flow
    rank <- rank(flow, ties.method="max")
    rank <- max(rank) - rank
    exceedtime <- 1*(rank / (length(flow) + 1))
    q_per <- sort(100 * flow / hist_mean, decreasing=FALSE)
    q_dis <- sort(flow, decreasing = F)
    exceed <- sort(exceedtime, decreasing=TRUE)
    fdc_df <- data.frame(exceed, q_per, q_dis, mean_flow = hist_mean)
  }
  
  return(fdc_df)
}
  