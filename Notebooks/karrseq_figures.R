#!/usr/bin/Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))


# 1. Prepare the background predicted Alus from master-table
exons.df <- fread("master_table_updated.tsv")

# Assign each row and unique ID
exons.df$id <- paste0("id", rownames(exons.df))

# Filter Alus
alu.df <- exons.df %>%
  dplyr::filter(repeat_pair == TRUE) %>%
  dplyr::filter(!str_detect(repeat_1, "No overlap")) %>%
  dplyr::filter(!str_detect(repeat_2, "No overlap")) %>%
  dplyr::filter(str_detect(repeat_1, "Alu") & str_detect(repeat_2, "Alu")) %>%
  mutate(id = paste(id, pair_type, sep = "_"))


background_counts.df <- alu.df %>%
  group_by(pair_type) %>%
  summarise(count = n()) %>%
  mutate(group = "background") %>%
  ungroup()


# 2. Load the repeats overlapped KARR-seq 
karrseq.dt <- fread("karrseq_hek293_g1_crssant.hits.bedpe", header = F)
karrseq.dt <- karrseq.dt %>%
  dplyr::rename(id = V7)

karrseq_counts.df <- karrseq.dt %>%
  mutate(pair_type = case_when(str_detect(id, "direct") ~ "direct",
                               TRUE ~ "inverted")) %>%
  group_by(pair_type) %>%
  summarise(count = n_distinct(id)) %>%
  mutate(group = "KARR-seq") %>%
  ungroup()


# 3. Compare groups
data_sub <- bind_rows(background_counts.df, karrseq_counts.df) %>%
  dplyr::rename(orientation = pair_type)

data_sub <- data.table(data_sub)

# (Re)calculate totals and percentages by group
data_sub[, total := sum(count), by = group]
data_sub[, perc := count / total]

# Create a contingency table directly from your data
#    - Rows: groups (e.g., "background", "KARR-seq")
#    - Columns: pair types (e.g., "inverted", "direct")

contingency_table <- dcast(
  data_sub, 
  group ~ orientation, 
  value.var = "count"
)

# Reorder rows to match desired order
contingency_table <- contingency_table[
  match(c("background", "KARR-seq"), contingency_table$group)
]

# Convert to matrix, excluding the 'group' column
contingency_mat <- as.matrix(contingency_table[, -1])

# Run Fisher's exact test
fisher_result <- fisher.test(contingency_mat)
p_value <- fisher_result$p.value
sig_label <- ifelse(p_value < 0.001, "***",
                    ifelse(p_value < 0.01, "**",
                           ifelse(p_value < 0.05, "*", "ns")))

print(fisher_result)


# Plot
karrseq_Fig1E <- ggplot(data_sub, aes(x = group, y = perc, fill = orientation)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  
  scale_fill_manual(
    name = "Alu-pair",
    values = c("inverted" = "red", "direct"  = "blue")
  ) +
  
  
  # new: show raw counts inside each bar segment
  geom_text(aes(label = count),
            position = position_fill(vjust = 0.5),
            color = "white", size = 3) +
  
  # (totals annotation removed!)
  
  scale_x_discrete(labels = c(
    "KARR-seq"  = "KARR-seq captured Alu",
    "background"       = "Background Alu"
  )) +
  theme_classic() +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.direction= "vertical"
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  
  geom_signif(
    comparisons = list(c("background", "KARR-seq")),
    annotations = sig_label,
    y_position  = 1.0,
    tip_length  = 0.02,
    textsize    = 3
  ) +
  
  scale_y_continuous(
    breaks  = seq(0, 1, by = 0.2),
    labels  = scales::percent_format(accuracy = 1),
    expand  = c(0, 0),
    limits  = c(0, 1.15)
  )

print(karrseq_Fig1E)










