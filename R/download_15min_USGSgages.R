# Download 15-minute resolution data from the USGS stream gages
# in southern New England

# Download 15-minute instantaneous discharge USGS stream gage data (just run once!)
start_date <- "05-31"
end_date <- "11-01"

## Load table with USGS stream gage codes in southern New England
#allgages_defol <- read_csv("allgages_defol.csv")

# Load table with info for the stream gages that had 15-min data 
gages_dischargeDuration <- read_csv("data/streamGagesCoverage.csv") 

# Check if data/Gage_Data/ folder exists, if not, create it
if(dir.exists("data/Gage_Data/")){
  print("Will download files to data/Gage_Data/ folder.") 
} else{
  dir.create("data/Gage_Data/")
  print("Created a data/Gage_Data/ folder to hold downloaded 15-min stream gage data.") 
}

# Download data 
for(g in 1:nrow(gages_dischargeDuration)){
  
  # Keep track of progress
  print(paste("Downloading gage",g,"of",nrow(gages_dischargeDuration),":",gages_dischargeDuration$STAID[g]))
  
  years <- gages_dischargeDuration$yr_start[g]:gages_dischargeDuration$yr_end[g]

   for(y in 1:length(years)){

     tmp_yr <- dataRetrieval::readNWISuv(gages_dischargeDuration$STAID[g],
                          parameterCd =  "00060",
                          startDate = paste0(years[y],"-",start_date),
                          endDate = paste0(years[y],"-",end_date))

     if(nrow(tmp_yr)!=0){
       if(y == 1){
         site_yrs <- tmp_yr
       } else {
         site_yrs <- rbind(site_yrs, tmp_yr)
       }
     }
   }
   write.csv(site_yrs, paste0("data/Gage_Data/",
                              gages_dischargeDuration$STAID[g],"_15min.csv"),
             row.names = FALSE)
}
