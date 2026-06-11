# ------------------------------------------------------------
# Code for prion-infected GtDeer images
# ------------------------------------------------------------

# ------------------------------------------------------------
# 1. Set the main folder path for prion-infected samples
# ------------------------------------------------------------

infected_folder <- "/Users/xbeatty/Desktop/Summer 2026/GtDeer_images"

# Check that R can see your region folders
list.files(infected_folder)

# ------------------------------------------------------------
# 2. Find all image files inside subfolders
# ------------------------------------------------------------

infected_image_files <- list.files(
  path = infected_folder,
  pattern = "\\.(tif|tiff|png|jpg|jpeg)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

length(infected_image_files)
infected_image_files

# ------------------------------------------------------------
# 3. Extract metadata from infected image file names
# ------------------------------------------------------------

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

infected_metadata

# ------------------------------------------------------------
# 4. Analyze every infected image using same threshold as controls
# ------------------------------------------------------------

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

infected_results

# ------------------------------------------------------------
# 5. Save infected results
# ------------------------------------------------------------

write.csv(
  infected_results,
  "/Users/xbeatty/Desktop/Summer 2026/GtDeer_PrP_single_channel_results.csv",
  row.names = FALSE
)

file.exists("/Users/xbeatty/Desktop/Summer 2026/GtDeer_PrP_single_channel_results.csv")

# ------------------------------------------------------------
# 6. Summary table by brain region
# ------------------------------------------------------------

infected_summary <- infected_results %>%
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

infected_summary

write.csv(
  infected_summary,
  "/Users/xbeatty/Desktop/Summer 2026/GtDeer_PrP_summary_by_region.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Plot infected PrP-positive area by brain region
# ------------------------------------------------------------

p_infected <- ggplot(infected_results, aes(x = brain_region, y = percent_prp_positive)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "PrP-positive area in prion-infected brain regions",
    x = "Brain region",
    y = "% PrP-positive area"
  )

print(p_infected)

ggsave(
  filename = "/Users/xbeatty/Desktop/Summer 2026/GtDeer_PrP_percent_positive_by_region.png",
  plot = p_infected,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Combine WT control and infected data
# ------------------------------------------------------------

control_results <- read.csv("/Users/xbeatty/Desktop/Summer 2026/WT_control_PrP_single_channel_results.csv")
infected_results <- read.csv("/Users/xbeatty/Desktop/Summer 2026/GtDeer_PrP_single_channel_results.csv")

control_results$condition <- "WT control"
infected_results$condition <- "Prion infected"

combined_results <- bind_rows(control_results, infected_results)

combined_results

# ------------------------------------------------------------
# 9. Plot WT control vs prion-infected samples
# ------------------------------------------------------------

p_compare <- ggplot(
  combined_results,
  aes(x = brain_region, y = percent_prp_positive, fill = condition)
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(color = condition),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.8,
    alpha = 0.6
  ) +
  theme_classic() +
  labs(
    title = "PrP-positive area: WT control vs GtDeer prion-infected brain regions",
    x = "Brain region",
    y = "% PrP-positive area",
    fill = "Condition",
    color = "Condition"
  )

print(p_compare)

ggsave(
  filename = "/Users/xbeatty/Desktop/Summer 2026/WT_vs_GtDeer_PrP_percent_positive_by_region.png",
  plot = p_compare,
  width = 8,
  height = 5,
  dpi = 300
)

file.exists("/Users/xbeatty/Desktop/Summer 2026/WT_vs_GtDeer_PrP_percent_positive_by_region.png")