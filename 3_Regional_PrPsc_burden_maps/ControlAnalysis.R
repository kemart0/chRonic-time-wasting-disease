install.packages("magick")
install.packages("EBImage")
install.packages("BiocManager")
install.packages("imager")
install.packages("magrittr")
install.packages("tidyverse")

library(magick)
library(EBImage)
library(imager)
library(tidyverse)
library(stringr)

# ------------------------------------------------------------
# 1. Set the main folder path
# ------------------------------------------------------------

control_folder <- "/Users/xbeatty/Desktop/Summer 2026/control_images"

# Check that R can see your region folders
list.files(control_folder)

# ------------------------------------------------------------
# 2. Find all image files inside subfolders
# ------------------------------------------------------------

image_files <- list.files(
  path = control_folder,
  pattern = "\\.(tif|tiff|png|jpg|jpeg)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

image_files

# ------------------------------------------------------------
# 3. Extract metadata from file names
# Example: WT_control_midbrain_4x_01.tif
# ------------------------------------------------------------

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

metadata

# ------------------------------------------------------------
# 4. Function to analyze single-channel PrP staining
# ------------------------------------------------------------

analyze_prp_single_channel <- function(file_path,
                                       tissue_threshold = 0.95,
                                       prp_threshold = 0.45,
                                       min_object_size = 20,
                                       max_object_size = 5000) {
  
  # Read image
  img <- readImage(file_path)
  
  # Convert image to numeric array
  img_array <- as.array(img)
  
  # Handle different image dimensions
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
  
  # Normalize intensity from 0 to 1
  img_gray <- normalize(img_gray)
  
  # Invert image so dark staining becomes high signal
  img_inverted <- 1 - img_gray
  
  # Tissue mask: removes mostly white background
  tissue_mask <- img_gray < tissue_threshold
  
  # PrP-positive mask: detects darker staining inside tissue
  prp_mask_raw <- tissue_mask & img_inverted > prp_threshold
  
  # Clean small noise
  prp_mask_clean <- opening(prp_mask_raw, makeBrush(3, shape = "disc"))
  
  # Label connected objects
  prp_objects <- bwlabel(prp_mask_clean)
  
  # Measure object features
  object_features <- computeFeatures.shape(prp_objects)
  
  # If no objects are detected
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
  
  # Object sizes
  object_sizes <- object_features[, "s.area"]
  
  # Keep objects within size range
  keep_objects <- which(
    object_sizes >= min_object_size &
      object_sizes <= max_object_size
  )
  
  # If no objects pass the size filter
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
  
  # Final PrP-positive mask
  prp_mask_final <- prp_objects %in% keep_objects
  
  # Measurements
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

results_list <- list()

for (i in seq_len(nrow(metadata))) {
  
  cat("Analyzing image", i, "of", nrow(metadata), ":", metadata$file_name[i], "\n")
  
  results_list[[i]] <- analyze_prp_single_channel(
    metadata$file_path[i],
    tissue_threshold = 0.95,
    prp_threshold = 0.45,
    min_object_size = 20,
    max_object_size = 5000
  )
}

analysis_results <- bind_rows(results_list)

results <- bind_cols(metadata, analysis_results)

results

write.csv(
  results,
  "/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_single_channel_results.csv",
  row.names = FALSE
)

p1 <- ggplot(results, aes(x = brain_region, y = percent_prp_positive)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "PrP-positive area in WT control brain regions",
    x = "Brain region",
    y = "% PrP-positive area"
  )

print(p1)

ggsave(
  filename = "/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_percent_positive_by_region.png",
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

file.exists("/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_percent_positive_by_region.png")

control_summary <- results %>%
  group_by(brain_region) %>%
  summarise(
    n_images = n(),
    mean_percent_prp = mean(percent_prp_positive, na.rm = TRUE),
    sd_percent_prp = sd(percent_prp_positive, na.rm = TRUE),
    median_percent_prp = median(percent_prp_positive, na.rm = TRUE),
    min_percent_prp = min(percent_prp_positive, na.rm = TRUE),
    max_percent_prp = max(percent_prp_positive, na.rm = TRUE),
    mean_deposits = mean(number_prp_deposits, na.rm = TRUE),
    mean_intensity = mean(mean_prp_intensity, na.rm = TRUE)
  )

control_summary

write.csv(
  control_summary,
  "/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_summary_by_region.csv",
  row.names = FALSE
)

p1 <- ggplot(results, aes(x = brain_region, y = percent_prp_positive)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  theme_classic() +
  labs(
    title = "WT control baseline: PrP-positive area by brain region",
    x = "Brain region",
    y = "% PrP-positive area"
  )

print(p1)

ggsave(
  filename = "/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_percent_positive_by_region.png",
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

top_outliers <- results %>%
  arrange(desc(percent_prp_positive)) %>%
  select(
    file_name,
    brain_region,
    percent_prp_positive,
    number_prp_deposits,
    mean_prp_intensity,
    tissue_area_pixels
  ) %>%
  head(20)

top_outliers

write.csv(
  top_outliers,
  "/Users/xbeatty/Desktop/Summer 2026/WT_control_top_possible_artifacts.csv",
  row.names = FALSE
)

p2 <- ggplot(results, aes(x = percent_prp_positive)) +
  geom_histogram(bins = 30) +
  theme_classic() +
  labs(
    title = "Distribution of detected PrP-positive area in WT controls",
    x = "% PrP-positive area",
    y = "Number of images"
  )

print(p2)

ggsave(
  filename = "/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_percent_positive_histogram.png",
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)
