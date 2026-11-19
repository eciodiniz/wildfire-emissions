# =========================================================
# POPULATION-EXPOSURE-WEIGHTED KDE (Emi_kton x population)
# Uses the ORIGINAL 1-ha points from the RDS (actual emission
# resolution), not the census-block (manzana) centroid. The
# population of each block is distributed among the emission
# points that fall inside it (point-in-polygon join).
# =========================================================

library(sf)
library(dplyr)
library(terra)
library(spatstat.geom)
library(spatstat.explore)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(patchwork)
library(data.table)
library(future.apply)

plan(multisession, workers = min(3, parallel::detectCores() - 1))

setwd("D:/AXA-Chile/Emissions_tests/Output")

out_dir <- "KDE_exposure_maps"
if (!dir.exists(out_dir)) dir.create(out_dir)

pts_crs <- 32719
manzana_id_col <- "ID"
pop_col <- "TOT_C_2024"

# =========================================================
# 1) LOAD DATA (1-ha points, actual resolution)
# =========================================================
rds <- readRDS("Pollutant_emissions_RDS/All_pollutants_list_corrected.rds")
pollutant_names <- names(rds)

# wide table: 1 row per point, coordinates + Emi_kton of each pollutant
wide <- rds[[pollutant_names[1]]][, .(FID, lon, lat, X, Y, Site)]
for (p in pollutant_names) wide[[p]] <- rds[[p]]$Emi_kton

# =========================================================
# 2) POINT-IN-POLYGON JOIN
# =========================================================
pop_sf <- st_read("Chile_population_2024.shp", quiet = TRUE)
pop_sf <- st_make_valid(pop_sf)
pop_vect <- vect(pop_sf)

pts_vect <- vect(wide, geom = c("lon", "lat"), crs = "EPSG:4326")
ext_res <- terra::extract(pop_vect, pts_vect)
wide[, manzana_id := ext_res[[manzana_id_col]]] 

# =========================================================
# 3) DISTRIBUTE BLOCK POPULATION AMONG ITS EMISSION POINTS
#    (each point gets an equal share of the population of the block
#    it falls in; preserves the exact population total and grid resolution)
# =========================================================
manzana_pop <- as.data.table(st_drop_geometry(pop_sf))[, .(manzana_id = get(manzana_id_col),
                                                           pop_manzana = get(pop_col))]
wide <- merge(wide, manzana_pop, by = "manzana_id", all.x = TRUE)
wide <- wide[!is.na(manzana_id)]

wide[, n_pts_manzana := .N, by = manzana_id]
wide[, pop_share := pop_manzana / n_pts_manzana]

# =========================================================
# 4) CHILE BASE MAP AND SITE POLYGONS 
# =========================================================
chile <- ne_countries(country = "Chile", scale = "medium", returnclass = "sf")
world <- ne_countries(scale = "medium", returnclass = "sf")
chile_neig <- world[st_touches(world, chile, sparse = FALSE), ]

site_names <- c("Talca", "Chillan", "concepcion", "curanilahue",
                "los_angeles", "san_antonio", "santiago", "valparaiso")

site_shp_dir <- "D:/Burn_probability_results/KDE map"

sites_list <- lapply(site_names, function(s) {
  shp <- st_read(file.path(site_shp_dir, paste0(s, ".shp")), quiet = TRUE)
  if (st_crs(shp)$epsg != 4326) shp <- st_transform(shp, 4326)
  shp
})
sites_polygons <- bind_rows(sites_list)
sites_polygons[8, ]$Name <- "Valparaiso"

labs_df <- st_centroid(sites_polygons) |>
  mutate(lon = st_coordinates(geometry)[, 1], lat = st_coordinates(geometry)[, 2])

# =========================================================
# 5) FUNCTION: KDE WEIGHTED BY Emi_kton (point) x pop_share, PER POLLUTANT
# =========================================================
run_kde_pollutant <- function(p, wide, chile, chile_neig, sites_polygons, labs_df, out_dir) {
  
  library(sf); library(terra); library(spatstat.geom); library(spatstat.explore)
  library(ggplot2); library(ggrepel); library(patchwork); library(data.table)
  
  w <- wide[[p]] * wide$pop_share
  keep <- w > 0
  if (sum(keep) < 5) { return(list(msg = paste("Too few points with weight > 0 for", p, "- skipping."), plot = NULL)) }
  
  coords <- as.matrix(wide[keep, .(X, Y)])
  weights <- w[keep]
  
  win <- owin(xrange = range(coords[, 1]), yrange = range(coords[, 2]))
  pp <- ppp(x = coords[, 1], y = coords[, 2], window = win)
  
  kde <- density.ppp(pp, sigma = 400, eps = 100, weights = weights)
  
  kde_raster <- rast(kde)
  crs(kde_raster) <- "EPSG:32719"
  
  writeRaster(kde_raster, file.path(out_dir, paste0("kde_exposure_", p, ".tif")), overwrite = TRUE)
  
  # ---- HOTSPOTS: top 10% (90th percentile) of density, as polygons ----
  p90 <- global(kde_raster, quantile, probs = 0.90, na.rm = TRUE)[1, 1]
  hotspots <- ifel(kde_raster >= p90, 1, NA)
  hotspots_poly <- as.polygons(hotspots, dissolve = TRUE, values = FALSE)
  hotspots_poly <- st_as_sf(hotspots_poly)
  hotspots_poly <- st_transform(hotspots_poly, 4326)
  hotspots_poly <- suppressWarnings(st_intersection(hotspots_poly, st_union(st_make_valid(sites_polygons))))
  hotspots_poly <- st_make_valid(hotspots_poly)
  
  study_bb <- st_bbox(sites_polygons)
  xpad <- (study_bb["xmax"] - study_bb["xmin"]) * 0.10
  ypad <- (study_bb["ymax"] - study_bb["ymin"]) * 0.08
  xlim_p <- c(study_bb["xmin"] - xpad, study_bb["xmax"] + xpad)
  ylim_p <- c(study_bb["ymin"] - ypad, study_bb["ymax"] + ypad)
  
  # ---- MAIN MAP (lightweight: polygons only, no country-wide tiles) ----
  main_map <- ggplot() +
    geom_sf(data = chile, fill = "grey96", color = "grey85", linewidth = 0.25) +
    geom_sf(data = sites_polygons, fill = "#E7DED3", color = "grey35", linewidth = 0.65) +
    geom_sf(data = hotspots_poly, fill = "firebrick2", color = "firebrick4",
            linewidth = 0.18, alpha = 0.90) +
    ggrepel::geom_text_repel(
      data = labs_df, aes(x = lon, y = lat, label = Name),
      color = "steelblue4", size = 3.5, fontface = "bold",
      segment.color = "grey50", segment.size = 0.5, segment.alpha = 0.6,
      box.padding = unit(0.6, "lines"), point.padding = unit(0.5, "lines"),
      force = 2, max.overlaps = Inf, bg.color = "ivory", bg.r = 0.2
    ) +
    labs(title = paste(p, " ")) +
    coord_sf(xlim = xlim_p, ylim = ylim_p, expand = FALSE) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "grey82", linetype = 2, linewidth = 0.25),
      panel.background = element_rect(fill = "lightcyan"),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 9)
    )
  
  p_final <- main_map
  
  #ggsave(file.path(out_dir, paste0("Chile_hotspots_exposure_", p, ".png")),
  #       plot = p_final, width = 12, height = 8, dpi = 300)
  
  list(msg = paste("Map created:", p), plot = p_final)
}

# =========================================================
# 6) PARALLEL EXECUTION (one pollutant per core) + TIMER
# =========================================================
t0 <- Sys.time()

results <- future_lapply(pollutant_names, function(p) {
  run_kde_pollutant(p, wide, chile, chile_neig, sites_polygons, labs_df, out_dir)
}, future.seed = TRUE)

t1 <- Sys.time()
elapsed <- difftime(t1, t0, units = "mins")

invisible(lapply(results, cat, "\n"))
cat("Total processing time:", round(as.numeric(elapsed), 2), "minutes\n")
cat("Done. Output in:", out_dir, "\n")

p_final <- results[[which(pollutant_names == "PM10")]]$plot   # choose the pollutant
print(p_final)
