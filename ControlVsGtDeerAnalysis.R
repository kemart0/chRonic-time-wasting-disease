# ============================================================
# WT control vs GtDeer prion-infected regional PrP burden map
# Blob size and color intensity both scale with median % PrP
# using ONE shared global scale across both conditions
# ============================================================

library(magick)
library(tidyverse)
library(scales)
library(patchwork)

# ------------------------------------------------------------
# 1. Load Allen atlas and make it grey
# ------------------------------------------------------------

atlas <- image_read("allen_coronal_atlas_section.jpg")

atlas_gray <- image_modulate(
  atlas,
  saturation = 0,
  brightness = 112
)

atlas_info <- image_info(atlas_gray)

atlas_width <- atlas_info$width
atlas_height <- atlas_info$height

atlas_raster <- as.raster(atlas_gray)

# ------------------------------------------------------------
# 2. Load WT control and GtDeer infected results
# ------------------------------------------------------------

control_results <- read.csv(
  "/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_single_channel_results.csv"
)

infected_results <- read.csv(
  "/Users/xbeatty/Desktop/Summer 2026/GtDeer_PrP_single_channel_results.csv"
)

# ------------------------------------------------------------
# 3. Summarize each condition by brain region
# Using median so it matches the boxplot center line
# ------------------------------------------------------------

control_summary <- control_results %>%
  group_by(brain_region) %>%
  summarise(
    n_images = n(),
    median_percent_prp = median(percent_prp_positive, na.rm = TRUE),
    mean_percent_prp = mean(percent_prp_positive, na.rm = TRUE),
    sd_percent_prp = sd(percent_prp_positive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(condition = "WT control")

infected_summary <- infected_results %>%
  group_by(brain_region) %>%
  summarise(
    n_images = n(),
    median_percent_prp = median(percent_prp_positive, na.rm = TRUE),
    mean_percent_prp = mean(percent_prp_positive, na.rm = TRUE),
    sd_percent_prp = sd(percent_prp_positive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(condition = "GtDeer prion infected")

control_summary
infected_summary

# ------------------------------------------------------------
# 4. Combine summaries so scaling is global across both groups
# ------------------------------------------------------------

all_summary <- bind_rows(control_summary, infected_summary)

# Shared value range for both color and size
global_min <- min(all_summary$median_percent_prp, na.rm = TRUE)
global_max <- max(all_summary$median_percent_prp, na.rm = TRUE)

global_min
global_max

# Add one shared size scaling variable
all_summary <- all_summary %>%
  mutate(
    size_scale = rescale(median_percent_prp, to = c(0.6, 1.5))
  )

all_summary

# Split back into control and infected
control_summary_scaled <- all_summary %>%
  filter(condition == "WT control")

infected_summary_scaled <- all_summary %>%
  filter(condition == "GtDeer prion infected")

# ------------------------------------------------------------
# 5. Function to create oval region blobs
# ------------------------------------------------------------

make_blob <- function(region, x_center, y_center, width, height, angle = 0, n = 100) {
  
  theta <- seq(0, 2 * pi, length.out = n)
  
  x <- (width / 2) * cos(theta)
  y <- (height / 2) * sin(theta)
  
  angle_rad <- angle * pi / 180
  
  x_rot <- x * cos(angle_rad) - y * sin(angle_rad)
  y_rot <- x * sin(angle_rad) + y * cos(angle_rad)
  
  tibble(
    brain_region = region,
    x = x_center + x_rot,
    y = y_center + y_rot
  )
}

# ------------------------------------------------------------
# 6. Function to build blobs with globally scaled size
# ------------------------------------------------------------

make_condition_blobs <- function(summary_df) {
  
  septum_scale <- summary_df %>%
    filter(brain_region == "septum") %>%
    pull(size_scale)
  
  hippocampus_scale <- summary_df %>%
    filter(brain_region == "hippocampus") %>%
    pull(size_scale)
  
  midbrain_scale <- summary_df %>%
    filter(brain_region == "midbrain") %>%
    pull(size_scale)
  
  cerebellum_scale <- summary_df %>%
    filter(brain_region == "cerebellum") %>%
    pull(size_scale)
  
  bind_rows(
    
    # Septum
    make_blob(
      region = "septum",
      x_center = 0.36 * atlas_width,
      y_center = 0.44 * atlas_height,
      width = 0.12 * atlas_width * septum_scale,
      height = 0.18 * atlas_height * septum_scale,
      angle = -10
    ),
    
    # Hippocampus
    make_blob(
      region = "hippocampus",
      x_center = 0.50 * atlas_width,
      y_center = 0.33 * atlas_height,
      width = 0.13 * atlas_width * hippocampus_scale,
      height = 0.09 * atlas_height * hippocampus_scale,
      angle = -20
    ),
    
    # Midbrain
    make_blob(
      region = "midbrain",
      x_center = 0.63 * atlas_width,
      y_center = 0.43 * atlas_height,
      width = 0.15 * atlas_width * midbrain_scale,
      height = 0.20 * atlas_height * midbrain_scale,
      angle = 10
    ),
    
    # Cerebellum
    make_blob(
      region = "cerebellum",
      x_center = 0.86 * atlas_width,
      y_center = 0.42 * atlas_height,
      width = 0.16 * atlas_width * cerebellum_scale,
      height = 0.28 * atlas_height * cerebellum_scale,
      angle = 5
    )
  ) %>%
    left_join(summary_df, by = "brain_region")
}

# ------------------------------------------------------------
# 7. Label positions
# ------------------------------------------------------------

base_labels <- tibble(
  brain_region = c("septum", "hippocampus", "midbrain", "cerebellum"),
  x = c(
    0.30 * atlas_width,
    0.48 * atlas_width,
    0.61 * atlas_width,
    0.84 * atlas_width
  ),
  y = c(
    0.63 * atlas_height,
    0.20 * atlas_height,
    0.66 * atlas_height,
    0.70 * atlas_height
  )
)

control_labels <- base_labels %>%
  left_join(control_summary_scaled, by = "brain_region")

infected_labels <- base_labels %>%
  left_join(infected_summary_scaled, by = "brain_region")

# ------------------------------------------------------------
# 8. Build blob data
# ------------------------------------------------------------

control_blobs <- make_condition_blobs(control_summary_scaled)
infected_blobs <- make_condition_blobs(infected_summary_scaled)

# ------------------------------------------------------------
# 9. Function to make one map
# ------------------------------------------------------------

make_brain_map <- function(blob_data, label_data, low_col, high_col, plot_title, legend_title) {
  
  ggplot() +
    annotation_raster(
      atlas_raster,
      xmin = 0,
      xmax = atlas_width,
      ymin = atlas_height,
      ymax = 0
    ) +
    geom_polygon(
      data = blob_data,
      aes(
        x = x,
        y = y,
        group = brain_region,
        fill = median_percent_prp
      ),
      color = NA,
      alpha = 0.65
    ) +
    geom_polygon(
      data = blob_data,
      aes(
        x = x,
        y = y,
        group = brain_region
      ),
      fill = NA,
      color = "black",
      linewidth = 0.4,
      alpha = 0.7
    ) +
    geom_label(
      data = label_data,
      aes(
        x = x,
        y = y,
        label = paste0(
          str_to_title(brain_region),
          "\n",
          round(median_percent_prp, 2),
          "%"
        )
      ),
      size = 3.8,
      fill = "white",
      label.size = 0.2,
      alpha = 0.9
    ) +
    scale_fill_gradient(
      name = legend_title,
      low = low_col,
      high = high_col,
      limits = c(global_min, global_max)
    ) +
    coord_fixed(
      xlim = c(0, atlas_width),
      ylim = c(atlas_height, 0),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    labs(
      title = plot_title,
      subtitle = "Grey = atlas background; blob size and color intensity = median % PrP-positive area"
    )
}

# ------------------------------------------------------------
# 10. Make maps
# ------------------------------------------------------------

p_control_map <- make_brain_map(
  blob_data = control_blobs,
  label_data = control_labels,
  low_col = "#bdfcff",
  high_col = "#00bcd4",
  plot_title = "WT control",
  legend_title = "Median % PrP-positive area"
)

p_infected_map <- make_brain_map(
  blob_data = infected_blobs,
  label_data = infected_labels,
  low_col = "#ffb3b3",
  high_col = "#d7191c",
  plot_title = "GtDeer prion infected",
  legend_title = "Median % PrP-positive area"
)

print(p_control_map)
print(p_infected_map)

# ------------------------------------------------------------
# 11. Combine maps
# ------------------------------------------------------------

combined_map <- p_control_map + p_infected_map +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Regional PrP burden map",
    subtitle = "Aqua = WT control; Red = GtDeer prion infected"
  )

print(combined_map)

# ------------------------------------------------------------
# 12. Save figure
# ------------------------------------------------------------

ggsave(
  filename = "/Users/xbeatty/Desktop/Summer 2026/WT_vs_GtDeer_brain_maps_GLOBAL_SIZE_SCALE.png",
  plot = combined_map,
  width = 16,
  height = 6,
  dpi = 300,
  bg = "white"
)

file.exists(
  "/Users/xbeatty/Desktop/Summer 2026/WT_vs_GtDeer_brain_maps_GLOBAL_SIZE_SCALE.png"
)