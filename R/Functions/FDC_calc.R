#' Flow Duration Curve
#'
#' Produces a flow duration curve plot 
#' @param flow daily streamflow time series
#' @param baseline whether or not to do baseline calculations
#' @references Gustard, A., Bullock, A., and Dixon, J.M. (1992). Report No. 108:
#'   Low flow estimation in the United Kingdom. 
#'   Oxfordshire, United Kingdom: Institute of Hydrology.

FDC_calc <- function(flow, mean_flow) {
  
  # Sort the flow values and calculate the FDC statistics
  rank <- rank(flow, ties.method="max")
  rank <- max(rank) - rank
  exceedtime <- 1*(rank / (length(flow) + 1))
  q_per <- sort(100 * flow / mean(flow), decreasing=FALSE)
  q_dis <- sort(flow, decreasing = F)
  exceed <- sort(exceedtime, decreasing=TRUE)
  fdc_df <- data.frame(exceed, q_per, q_dis, mean_flow = mean(flow))
  
  return(fdc_df)
}
  