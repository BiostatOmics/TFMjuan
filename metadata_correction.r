library(tidyverse)

full_cell_metadata <- read.csv("B:/TFM/AMP_PD_cluster/04-meta_data/full_cell_metadata.tsv", row.names=1)


braak_score <- read.csv("B:/TFM/AMP_PD_cluster/04-meta_data/releases_2023_v4release_1027_clinical_LBD_Cohort_Path_Data.csv")


participants = str_to_upper(unique(full_cell_metadata$participant_id))

braak_score <- braak_score %>%
  mutate(participant_id = str_to_upper(participant_id)) %>% 
  filter(participant_id %in% participants) %>%
  select(c(participant_id, path_braak_nft, path_braak_lb))

full_cell_metadata <- left_join(full_cell_metadata, braak_score, by = "participant_id")

rownames(full_cell_metadata) <- full_cell_metadata$barcodekey

full_cell_metadata <- full_cell_metadata %>%
  select(-c(nCount_RNA, nFeature_RNA, leiden_1, ROIGroupFine, dissection, supercluster_term, cell_type)) %>% 
  rename(case_control = case_control_other_latest)

write.csv(
  full_cell_metadata,
  "04-meta_data/corrected_full_cell_metadata.tsv",
  sep = "\t",
  row.names = T,
  col.names = T,
)
