---
title: "README: Generating Regional PrP Burden Maps"
author: "Xuan Beatty"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 3
    number_sections: true
---

# Overview

This README describes the R workflow used to generate a regional PrP burden map comparing **WT control** and **GtDeer prion-infected** brain samples.

The final figure uses a grey Allen mouse brain atlas as an anatomical reference. Regional overlays are placed on the atlas for the **septum**, **hippocampus**, **midbrain**, and **cerebellum**. The overlays are colored by condition:

- **Aqua** = WT control
- **Red** = GtDeer prion infected

The **blob size** and **color intensity** represent the **median percentage of PrP-positive area** for each brain region. Larger and darker blobs indicate greater regional PrP burden.

# Input data

## Image organization

The image analysis begins with PrP IHC images organized by brain region. The control and infected samples are stored in separate folders with matching subfolder structures.

```text
/Users/xbeatty/Desktop/Summer 2026/control_images/
├── cerebellum_4x
├── hippocampus_4x
├── midbrain_4x
└── septum_4x

/Users/xbeatty/Desktop/Summer 2026/GtDeer_images/
├── cerebellum_4x
├── hippocampus_4x
├── midbrain_4x
└── septum_4x
```

Image files are named in a way that includes sample information and brain region. Example:

```text
WT_control_midbrain_4x_01.tif
```

This name contains:

| File name part | Meaning |
|---|---|
| `WT` | Genotype/sample type |
| `control` | Experimental group |
| `midbrain` | Brain region |
| `4x` | Magnification |
| `01` | Image number |

# Image analysis workflow

Each image is analyzed to estimate the amount of PrP-positive staining. Because the images are treated as single-channel images, the analysis is based on staining intensity rather than color deconvolution.

For each image, the code performs the following steps:

1. Loads the image.
2. Converts the image to a numeric grayscale matrix.
3. Normalizes pixel intensity from 0 to 1.
4. Inverts the image so darker staining becomes higher signal.
5. Creates a tissue mask to remove bright slide background.
6. Thresholds dark signal within tissue to identify PrP-positive pixels.
7. Removes small noise using morphological opening.
8. Labels connected PrP-positive objects.
9. Calculates PrP-positive area as a percentage of total tissue area.

The main measurement is:

```text
% PrP-positive area = PrP-positive pixels / total tissue pixels × 100
```

# R packages

```{r packages, eval=FALSE}
library(magick)
library(EBImage)
library(tidyverse)
library(stringr)
library(scales)
library(patchwork)
```

# Step 1: Analyze WT control images

## Locate control images

```{r control-files, eval=FALSE}
control_folder <- "/Users/xbeatty/Desktop/Summer 2026/control_images"

image_files <- list.files(
  path = control_folder,
  pattern = "\\.(tif|tiff|png|jpg|jpeg)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

metadata <- tibble(file_path = image_files) %>%
  mutate(
    file_name = basename(file_path),
    file_name_no_ext = str_remove(file_name, "\\.[^.]+$")
  ) %>%
  separate(
    file_name_no_ext,
    into = c("genotype", "group", "brain_region", "magnification", "image_number"),
    sep = "_",
    remove = FALSE
  )
```

## Function to quantify PrP staining

```{r analysis-function, eval=FALSE}
analyze_prp_single_channel <- function(file_path,
                                       tissue_threshold = 0.95,
                                       prp_threshold = 0.45,
                                       min_object_size = 20,
                                       max_object_size = 5000) {
  
  img <- readImage(file_path)
  img_array <- as.array(img)
  
  if (length(dim(img_array)) == 2) {
    img_gray <- img_array
  } else if (length(dim(img_array)) == 3) {
    if (dim(img_array)[3] >= 3) {
      img_gray <- apply(img_array[,,1:3], c(1, 2), mean)
    } else {
      img_gray <- img_array[,,1]
    }
  } else if (length(dim(img_array)) == 4) {
    img_gray <- img_array[,,1,1]
  } else {
    stop("Unsupported image dimensions for file: ", file_path)
  }
  
  img_gray <- normalize(img_gray)
  img_inverted <- 1 - img_gray
  
  tissue_mask <- img_gray < tissue_threshold
  prp_mask_raw <- tissue_mask & img_inverted > prp_threshold
  prp_mask_clean <- opening(prp_mask_raw, makeBrush(3, shape = "disc"))
  
  prp_objects <- bwlabel(prp_mask_clean)
  object_features <- computeFeatures.shape(prp_objects)
  
  if (nrow(object_features) == 0) {
    return(tibble(
      tissue_area_pixels = sum(tissue_mask),
      prp_positive_pixels = 0,
      percent_prp_positive = 0,
      number_prp_deposits = 0,
      mean_deposit_size_pixels = NA,
      total_deposit_area_pixels = 0,
      mean_prp_intensity = NA
    ))
  }
  
  object_sizes <- object_features[, "s.area"]
  
  keep_objects <- which(
    object_sizes >= min_object_size &
      object_sizes <= max_object_size
  )
  
  if (length(keep_objects) == 0) {
    return(tibble(
      tissue_area_pixels = sum(tissue_mask),
      prp_positive_pixels = 0,
      percent_prp_positive = 0,
      number_prp_deposits = 0,
      mean_deposit_size_pixels = NA,
      total_deposit_area_pixels = 0,
      mean_prp_intensity = NA
    ))
  }
  
  prp_mask_final <- prp_objects %in% keep_objects
  
  tissue_area_pixels <- sum(tissue_mask)
  prp_positive_pixels <- sum(prp_mask_final)
  percent_prp_positive <- prp_positive_pixels / tissue_area_pixels * 100
  number_prp_deposits <- length(keep_objects)
  mean_deposit_size_pixels <- mean(object_sizes[keep_objects], na.rm = TRUE)
  total_deposit_area_pixels <- sum(object_sizes[keep_objects], na.rm = TRUE)
  mean_prp_intensity <- mean(img_inverted[prp_mask_final], na.rm = TRUE)
  
  tibble(
    tissue_area_pixels = tissue_area_pixels,
    prp_positive_pixels = prp_positive_pixels,
    percent_prp_positive = percent_prp_positive,
    number_prp_deposits = number_prp_deposits,
    mean_deposit_size_pixels = mean_deposit_size_pixels,
    total_deposit_area_pixels = total_deposit_area_pixels,
    mean_prp_intensity = mean_prp_intensity
  )
}
```

## Run analysis on WT control images

```{r control-analysis, eval=FALSE}
results_list <- list()

for (i in seq_len(nrow(metadata))) {
  
  cat("Analyzing control image", i, "of", nrow(metadata), ":", metadata$file_name[i], "\n")
  
  results_list[[i]] <- analyze_prp_single_channel(
    metadata$file_path[i],
    tissue_threshold = 0.95,
    prp_threshold = 0.45,
    min_object_size = 20,
    max_object_size = 5000
  )
}

analysis_results <- bind_rows(results_list)
control_results <- bind_cols(metadata, analysis_results)

write.csv(
  control_results,
  "/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_single_channel_results.csv",
  row.names = FALSE
)
```

# Step 2: Analyze GtDeer prion-infected images

The same analysis function is applied to the GtDeer prion-infected images. The same threshold settings are used so the WT control and infected samples are comparable.

```{r infected-analysis, eval=FALSE}
infected_folder <- "/Users/xbeatty/Desktop/Summer 2026/GtDeer_images"

infected_image_files <- list.files(
  path = infected_folder,
  pattern = "\\.(tif|tiff|png|jpg|jpeg)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

infected_metadata <- tibble(file_path = infected_image_files) %>%
  mutate(
    file_name = basename(file_path),
    file_name_no_ext = str_remove(file_name, "\\.[^.]+$")
  ) %>%
  separate(
    file_name_no_ext,
    into = c("genotype", "group", "brain_region", "magnification", "image_number"),
    sep = "_",
    remove = FALSE,
    fill = "right"
  )

infected_results_list <- list()

for (i in seq_len(nrow(infected_metadata))) {
  
  cat("Analyzing infected image", i, "of", nrow(infected_metadata), ":", infected_metadata$file_name[i], "\n")
  
  infected_results_list[[i]] <- analyze_prp_single_channel(
    infected_metadata$file_path[i],
    tissue_threshold = 0.95,
    prp_threshold = 0.45,
    min_object_size = 20,
    max_object_size = 5000
  )
}

infected_analysis_results <- bind_rows(infected_results_list)
infected_results <- bind_cols(infected_metadata, infected_analysis_results)

write.csv(
  infected_results,
  "/Users/xbeatty/Desktop/Summer 2026/GtDeer_PrP_single_channel_results.csv",
  row.names = FALSE
)
```

# Step 3: Summarize PrP burden by brain region

After image analysis, the WT control and GtDeer infected results are summarized by brain region.

The **median percent PrP-positive area** is used because it matches the center line of the boxplot and is less affected by unusually high artifact-heavy or highly stained images.

```{r summarize-results, eval=FALSE}
control_results <- read.csv(
  "/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_single_channel_results.csv"
)

infected_results <- read.csv(
  "/Users/xbeatty/Desktop/Summer 2026/GtDeer_PrP_single_channel_results.csv"
)

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
```

# Step 4: Load the Allen atlas background

The Allen atlas image is used as a grey anatomical reference. The original atlas is converted to grayscale so the PrP burden overlays are easier to see.

```{r load-atlas, eval=FALSE}
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
```

# Step 5: Create regional atlas overlays

Approximate oval overlays are manually positioned over the four analyzed brain regions:

- Septum
- Hippocampus
- Midbrain
- Cerebellum

The overlays do not represent exact atlas segmentation masks. They are visual markers for the regions analyzed from the image folders.

```{r blob-function, eval=FALSE}
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
```

# Step 6: Scale blob size and color intensity

The WT control and GtDeer infected summaries are combined before scaling. This creates one shared scale across both groups.

This is important because it prevents the control and infected maps from being scaled separately.

```{r scaling, eval=FALSE}
all_summary <- bind_rows(control_summary, infected_summary)

global_min <- min(all_summary$median_percent_prp, na.rm = TRUE)
global_max <- max(all_summary$median_percent_prp, na.rm = TRUE)

all_summary <- all_summary %>%
  mutate(
    size_scale = rescale(median_percent_prp, to = c(0.6, 1.5))
  )

control_summary_scaled <- all_summary %>%
  filter(condition == "WT control")

infected_summary_scaled <- all_summary %>%
  filter(condition == "GtDeer prion infected")
```

In the final map:

| Visual feature | Meaning |
|---|---|
| Aqua overlay | WT control |
| Red overlay | GtDeer prion infected |
| Larger blob | Higher median % PrP-positive area |
| Darker color | Higher median % PrP-positive area |
| Grey atlas | Anatomical background |

# Step 7: Build condition-specific blobs

```{r build-blobs, eval=FALSE}
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
    make_blob(
      region = "septum",
      x_center = 0.36 * atlas_width,
      y_center = 0.44 * atlas_height,
      width = 0.12 * atlas_width * septum_scale,
      height = 0.18 * atlas_height * septum_scale,
      angle = -10
    ),
    make_blob(
      region = "hippocampus",
      x_center = 0.50 * atlas_width,
      y_center = 0.33 * atlas_height,
      width = 0.13 * atlas_width * hippocampus_scale,
      height = 0.09 * atlas_height * hippocampus_scale,
      angle = -20
    ),
    make_blob(
      region = "midbrain",
      x_center = 0.63 * atlas_width,
      y_center = 0.43 * atlas_height,
      width = 0.15 * atlas_width * midbrain_scale,
      height = 0.20 * atlas_height * midbrain_scale,
      angle = 10
    ),
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

control_blobs <- make_condition_blobs(control_summary_scaled)
infected_blobs <- make_condition_blobs(infected_summary_scaled)
```

# Step 8: Add labels

```{r labels, eval=FALSE}
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
```

# Step 9: Generate brain maps

```{r map-function, eval=FALSE}
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
```

# Step 10: Combine and save the final figure

```{r save-map, eval=FALSE}
combined_map <- p_control_map + p_infected_map +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Regional PrP burden map",
    subtitle = "Aqua = WT control; Red = GtDeer prion infected"
  )

print(combined_map)

ggsave(
  filename = "/Users/xbeatty/Desktop/Summer 2026/WT_vs_GtDeer_brain_maps_GLOBAL_SIZE_SCALE.png",
  plot = combined_map,
  width = 16,
  height = 6,
  dpi = 300,
  bg = "white"
)
```

# Interpretation of the final image

The final map compares regional PrP burden between WT control and GtDeer prion-infected brain samples.

- WT control regions are shown in aqua.
- GtDeer prion-infected regions are shown in red.
- The grey atlas provides anatomical context.
- Larger and darker overlays indicate higher median percent PrP-positive area.

In the example final figure, the WT control regions show low median PrP-positive area, while the GtDeer prion-infected regions show higher median PrP-positive area across all four analyzed brain regions. This indicates greater PrP deposition in the prion-infected samples compared with controls.

# Limitations

This map is a visual regional summary, not a precise anatomical registration.

Important limitations:

- The region overlays are approximate manually placed blobs.
- The overlays are not exact Allen atlas segmentation masks.
- The original microscopy images are not directly warped or registered to the atlas.
- The map summarizes region-level PrP burden from image folders rather than pixel-perfect anatomical coordinates.

For a more anatomically precise map, exact atlas region masks or image registration software would be needed.

# Suggested figure caption

**Regional PrP burden map comparing WT control and GtDeer prion-infected brain regions.** PrP IHC images were analyzed by brain region to calculate the percentage of PrP-positive tissue area. Median percent PrP-positive area was summarized for the cerebellum, hippocampus, midbrain, and septum. Values were visualized on a grey Allen mouse brain atlas background using approximate region overlays. Aqua indicates WT control samples and red indicates GtDeer prion-infected samples. Blob size and color intensity represent median percent PrP-positive area, with larger and darker overlays indicating greater regional PrP burden.
