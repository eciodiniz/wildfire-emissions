##############################################################
# POPULATION EXPOSURE TO POLLUTANTS - CUMULATIVE CURVE
# Based on the per-manzana table 
##############################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(snakecase)
library(data.table)
library(scales)

# =========================
# PATHS
# =========================
base_path <- "D:/AXA-Chile/Emissions_tests/Output"
out_dir   <- file.path(base_path, "Pollutant_Population_RDS")
fig_dir   <- file.path(base_path, "Exposure_figures")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

path_in <- file.path(out_dir, "Emissions_by_manzana.rds")

df_all <- as.data.table(readRDS(path_in))

pollutant_names <- c("CO2", "PM10", "PM2.5", "NH3", "NOx", "SO2", "CH4", "CO")

# Site -> formatted label map, computed once (avoids recomputing for each pollutant)
site_map <- setNames(to_any_case(unique(df_all$Site), case = "title"), unique(df_all$Site))

cum_plots <- list()   # cumulative curves per pollutant (view: cum_plots[["PM10"]])

# =========================
# LOOP BY POLLUTANT
# =========================
for (p in pollutant_names) {
  
  # "exposed" = manzana population with emission of this pollutant > 0
  df_p <- df_all[, .(Site, pop, exposed_flag = get(p) > 0)]
  
  expo_tab <- df_p |>
    group_by(Site) |>
    summarise(
      total_pop = sum(pop, na.rm = TRUE),
      exposed_pop = sum(pop[exposed_flag], na.rm = TRUE),
      relative_exposure = sum(pop[exposed_flag], na.rm = TRUE) / sum(pop, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(exposed_pop)) |>
    mutate(
      exposed_million = exposed_pop / 1e6,
      relative_percent = relative_exposure * 100,
      site_label = site_map[Site]
    )
  
   #write.csv(expo_tab, file.path(fig_dir, paste0("Exposure_summary_", p, ".csv")), row.names = FALSE)
  
  # ---- PLOT 3: LORENZ CURVE (cumulative % population vs. cumulative % emission) ----
  # Manzanas sorted from HIGHEST to LOWEST emission (version that produced the PM10 figure)
  df_cum <- df_all[, .(Site, pop, emi = get(p))]
  df_cum <- df_cum[order(Site, -emi)]
  df_cum[, cum_pop := cumsum(pop) / sum(pop), by = Site]
  df_cum[, cum_emi := cumsum(emi) / sum(emi), by = Site]
  df_cum[, site_label := site_map[Site]]
  
  p3 <- ggplot(df_cum, aes(x = cum_pop, y = cum_emi, color = site_label)) +
    geom_line(linewidth = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = paste("Cumulative exposure curve -", p),
         x = "Cumulative % of population",
         y = "Cumulative % of emission",
         color = "Site") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  
  #ggsave(file.path(fig_dir, paste0("Cumulative_curve_", p, ".jpeg")),
  #       plot = p3, width = 8, height = 6, dpi = 300)
  
  cum_plots[[p]] <- p3 
  
  cat("Pollutant processed:", p, "\n")
}

cat("Figures and tables saved to:", fig_dir, "\n")

# View curve (change the pollutant)
print(cum_plots[["PM2.5"]])
