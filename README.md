# Mapping potential multi-pollutant wildfire emissions and population exposure in R

This repository contains the R scripts for mapping the emission potential of eight wildfire pollutants (PM2.5, PM10, CO, CO2, CH4, NOx, SO2, and NH3), assessing how these emissions are distributed across the population, and identifying hotspots of population-weighted emission potential.

The input data are available on Zenodo: [DOI to be added].

## Scripts

**[Emission Potential Distribution](Emission%20Potential%20Distribution)** contains the script that maps the spatial distribution of wildfire emission potential for a selected pollutant across Central Chile.

**[Emission Hotspots](Emission%20Hotspots)** contains the script that identifies hotspots of population-weighted emission potential for each pollutant. It combines the 1-ha emission points with census-block population, estimates a kernel density surface weighted by emission × population, and defines hotspots as the areas above the 90th percentile of density.

**[Cumulative Exposure](Cumulative%20Exposure)** contains the script that assesses how the emission potential of each pollutant is distributed across the population of each study area. Using emissions and population aggregated by census block, it builds cumulative exposure (Lorenz) curves by ranking blocks from highest to lowest emission and plotting the cumulative share of population against the cumulative share of emission for each site.
