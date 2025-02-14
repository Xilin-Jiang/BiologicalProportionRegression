library(data.table)
library(tidyr)
library(dplyr)
## code to prepare `create_nonidentify_data` dataset goes here
all_phenotypes <- fread("../../PROJECTS/protein_components/protein_components/Mediation_analysis/output/Weissbrod_2022_meta_file.txt")
traits_list <- names(all_phenotypes)[7:ncol(all_phenotypes)]
# load protein level (note in this analysis, it is important to distinguish first and second visit )
long_olink_data_2923 <- read.table("../../PROJECTS/protein_components/protein_components/raw_data/olink_data_20231204.txt" , header = T, sep = "\t", quote = "")
# only use ins_index == 0 as the later visit are COVID-19 recall, will cause some issue
long_olink_data_2923 <- long_olink_data_2923 %>%
  mutate(protein_id = paste0("pid", protein_id)) %>%
  filter(ins_index == 0)
outcome_protein_list <- long_olink_data_2923 %>%
  select(protein_id) %>%
  distinct()

# filter to british ancestry
british_ancestry <- all_phenotypes %>%
  filter(genetic_white_british) %>%
  pull(eid) %>%
  unique()

protein_data_all <- long_olink_data_2923 %>%
  filter(eid %in% british_ancestry) %>%
  filter(ins_index == 0)

# pivot it to be wider -- we could do this because only instance == 0 are kept.
protein_ins0_british <- protein_data_all %>%
  pivot_wider(names_from = protein_id, values_from = result)

X_id <- "pid1980" # we are expecting strong eta_j for pid767 (this two proteins have really high correlation)
X <- protein_ins0_british %>%
  select(!!X_id) %>%
  as.matrix() %>%
  scale()

Y_id <- "pid1968"
Y <- protein_ins0_british %>%
  select(!!Y_id) %>%
  as.matrix() %>%
  scale()

Z <- protein_ins0_british %>%
  select(-!!X_id, -!!Y_id, -eid, - ins_index) %>%
  select(1:200) %>%
  as.matrix()
Z <- apply(Z, 2, scale)

# summary level covariance -- could be shared
EXAMPLE_COV <- padding_X_Y_cov(X, Y, Z, X_padding_num = 1, Y_padding_num = 1)
usethis::use_data(EXAMPLE_COV, overwrite = TRUE)

# add noise to X, Y, and Z
indiviual_ids <- sample(1:length(X), size = 5000, replace = F)
X_EXAMPLE <- X[indiviual_ids] + rnorm(length(indiviual_ids)) * 0.03 * sd(X, na.rm = T)
Y_EXAMPLE <- Y[indiviual_ids] + rnorm(length(indiviual_ids)) * 0.3 * sd(Y, na.rm = T)
Z_EXAMPLE <- Z[indiviual_ids, ] + 0.5 * matrix(rnorm(length(indiviual_ids) * dim(Z)[2]), nrow = length(indiviual_ids) , ncol = dim(Z)[2])
usethis::use_data(X_EXAMPLE, overwrite = TRUE)
usethis::use_data(Y_EXAMPLE, overwrite = TRUE)
usethis::use_data(Z_EXAMPLE, overwrite = TRUE)



