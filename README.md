# Impacts of defoliation by gypsy moth caterpillars on streamflow in Southern New England

This code is associated with the publication: Smith-Tripp, S., A. Griffith, V. Pasquarella, J.H. Matthes. Impacts of a regional multi-year insect defoliation event on growing-season runoff ratios and instantaneous streamflow characteristics. In Review.

This analysis assesses the connection between defoliation in southern New England during a 2015-2017 outbreak of gypsy moth (*Lymantria dispar*) caterpillars and growing season (June-September) anomalies in watershed runoff ratio, water yield, and instaneous streamflow characterics derived from gage-specific flow duration curves.

The code is organized into R files with analysis pieces: 
1. 01_loadData.R: Load/download defoliation, stream gage, and precipitaiton data.
2. 02_calculateBudgetsFDC.R: Calculate the baseline (1995-2014) and 2015-2017 growing season statistics for yield, precipitation, and runoff ratios, and FDC curves. 
3. 03_budgetModelsFigures.R: Estimate statistical models and make figures for departures in growing-season runoff ratios/yield/precipitation with defoliation intensity.
4. 04_fdcModelsFigures.R: Estimate statistical models and make figures for departures in growing-season FDC statistics with defoliation intensity.
5. 05_studyLocationMap.R: Make a map of the study area within the U.S. and with watershed shapfiles showing defoliation intensity in 2015-2017 growing seasons.

