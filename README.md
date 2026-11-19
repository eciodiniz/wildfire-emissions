# Mapping potential multi-pollutant wildfire emissions and population exposure in R

This repository contains the R scripts used in the study *Potential multi-pollutant wildfire emissions and population exposure across Central Chile*. The scripts map the emission potential of eight wildfire pollutants (PM2.5, PM10, CO, CO2, CH4, NOx, SO2 and NH3), assess how these emissions are distributed across the population, and identify hotspots of population-weighted emission potential.

The input data are available on Zenodo: [DOI to be added].

## Scripts

**Wildfire_emission_potential_v2.R** maps where wildfire emission potential is highest across Central Chile. It reads the 1-ha emission points of each pollutant (`All_pollutants_list_corrected.rds`), sums them within 0.01° (~1 km) cells, and overlays the boundaries of the eight study areas to produce the regional emission maps (JPEG, 300 dpi).

**Cumulative_Curves_Exposure.R** assesses how the emission potential of each pollutant is distributed across the population of each study area. Using emissions and population aggregated by census block (`Emissions_by_manzana.rds`), it derives the exposed and relative population per site and builds cumulative exposure curves by ranking blocks from highest to lowest emission.

**Emission_Exposure_Hotspots.R** identifies where emission potential and population coincide most strongly. It distributes the population of each census block (`Chile_population_2024.shp`) equally among its 1-ha emission points, estimates a kernel density surface weighted by emission × population (400 m bandwidth, 100 m resolution), and delimits hotspots at the 90th percentile of density (GeoTIFF rasters and maps per pollutant).
