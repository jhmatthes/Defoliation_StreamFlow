# Impacts of defoliation by gypsy moth caterpillars on streamflow in Southern New England

This code is associated with the publication: Smith-Tripp, S., A. Griffith, V. Pasquarella, J.H. Matthes. TBD. Impacts of a regional multi-year insect defoliation on seasonal water yield and instantaneous streamflow characteristics. In Review.

This analysis assesses the connection between defoliation in southern New England during a 2015-2017 outbreak of gypsy moth (*Lymantria dispar*) caterpillars and growing season (June-October) anomalies in water yield, the yield-to-precipitation ratio, and instaneous streamflow characterics from flow duration curves.

The code is organized by: 
1. Load defoliation values extracted for 2015-2017 by overlaying sub-watershed shapefiles on Landsat data product. 
2. Download 15-minute streamflow data for 1995-2017 from USGS stream gages in the Southern New England region. 
3. Download daymet precipitation data for 1995-2017 co-located with the USGS stream gages.
4. Process the streamflow and precipitation data into seasonal sums.
5. Calculate runoff ratio, yield, and precipitation anomalies at each stream gage for the 2015-2017 time period compared to a 1995-2014 mean baseline period.
6. Assess runoff anomaly ~ defoliation relationships.
7. Calulate stream flow duration curves at each stream gage for the 2015-2017 period and compare to a 1995-2014 baseline flow duration curve.
8. Assess FDC anomaly ~ defoliation relationships.

