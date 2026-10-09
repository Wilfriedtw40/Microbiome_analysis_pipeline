#####setup work environment
knitr::opts_knit$set(root.dir = 'C:/Users/riri4/OneDrive/Bureau/Stage_Bernard/16S_NCCR')

######install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("dada2", quietly = TRUE))
  BiocManager::install("dada2")
if (!requireNamespace("VennDiagram", quietly = TRUE))
  BiocManager::install("VennDiagram")

######load required libraries
rm(list=ls()) 
library(dada2)
library(Biostrings)
library(VennDiagram)
library(tidyverse)

######define paths, core variables and create directories
my_root <- "."
my_plots <- file.path(my_root, "plots", "02_denoising")
common_files <- "C:/Users/riri4/OneDrive/Bureau/Stage_Bernard/16S_NCCR/data"

dir.create(my_plots, showWarnings = F, recursive = T)

######define helper functions 
getN <- function(x) {
  sum(getUniques(x))
}

count_variants <- function(x) {
  length(getUniques(x))
}

extract_inferred_sequences <- function(
    dada_list,
    direction,
    pooling_method
) {
  
  names(dada_list) <- sample.names
  
  purrr::imap_dfr(
    dada_list,
    function(dada_object, sample_name) {
      
      uniques <- getUniques(dada_object)
      
      tibble(
        sample = sample_name,
        direction = direction,
        pooling = pooling_method,
        sequence = names(uniques),
        abundance = as.numeric(uniques)
      )
    }
  )
}

safe_jaccard <- function(x, y) {
  
  union_size <- sum(x | y)
  
  if (union_size == 0) {
    return(NA_real_)
  }
  
  sum(x & y) / union_size
}

getMergedN <- function(x) {
  sum(x$abundance)
}

getMergedVariants <- function(x) {
  length(unique(x$sequence))
}

summarise_merging <- function(
    dadaF,
    dadaR,
    mergers,
    method
) {
  
  if (
    length(dadaF) != length(sample.names) ||
    length(dadaR) != length(sample.names) ||
    length(mergers) != length(sample.names)
  ) {
    stop(
      "Unexpected number of samples for method: ",
      method
    )
  }
  
  tibble(
    sample = sample.names,
    denoising_method = method,
    
    denoised_F = unname(
      vapply(dadaF, getN, numeric(1))
    ),
    
    denoised_R = unname(
      vapply(dadaR, getN, numeric(1))
    ),
    
    merged_reads = unname(
      vapply(mergers, getMergedN, numeric(1))
    ),
    
    merged_variants = unname(
      vapply(
        mergers,
        getMergedVariants,
        integer(1)
      )
    )
  ) %>%
    mutate(
      # Conservative estimate of reads available in both directions
      pairable_denoised_reads = pmin(
        denoised_F,
        denoised_R
      ),
      
      perc_merged = if_else(
        pairable_denoised_reads > 0,
        100 * merged_reads /
          pairable_denoised_reads,
        NA_real_
      )
    )
}

summarise_dada_results <- function(
    dada_list,
    direction,
    pooling_method
) {
  
  if (length(dada_list) != length(sample.names)) {
    stop(
      "Unexpected number of samples for ",
      direction,
      " / ",
      pooling_method
    )
  }
  names(dada_list) <- sample.names
  
  tibble(
    sample = sample.names,
    direction = direction,
    pooling = pooling_method,
    
    denoised_reads = unname(
      vapply(
        dada_list,
        getN,
        numeric(1)
      )
    ),
    
    inferred_variants = unname(
      vapply(
        dada_list,
        count_variants,
        integer(1)
      )
    )
  )
}

#Chimera QC analysis function
flag_chimera_results <- function(
    original_table,
    nochim_table,
    denoising_method,
    sample_chimera_threshold = 20,
    candidate_relative_threshold = 5,
    candidate_chimera_share_threshold = 50
) {
  #input validation
  if (is.null(rownames(original_table)) ||
      is.null(rownames(nochim_table))) {
    stop("Both ASV tables must contain sample names as row names.")
  }
  
  if (!identical(
    rownames(original_table),
    rownames(nochim_table)
  )) {
    stop(
      "Original and chimera-free tables have ",
      "different sample names or sample order."
    )
  }
  
  if (!all(
    colnames(nochim_table) %in%
    colnames(original_table)
  )) {
    stop(
      "The chimera-free table contains sequences ",
      "not present in the original table."
    )
  }
  #ASVs removed by chimera filtering
  removed_sequences <- setdiff(
    colnames(original_table),
    colnames(nochim_table)
  )
  if (length(removed_sequences) == 0) {
    return(
      list(
        sample_summary = tibble(),
        flagged_candidates = tibble(),
        all_removed_candidates = tibble()
      )
    )
  }
  
  #Convert removed-ASV counts to long format
  removed_long <- original_table[
    ,
    removed_sequences,
    drop = FALSE
  ] %>%
    as.data.frame(
      check.names = FALSE
    ) %>%
    rownames_to_column("sample") %>%
    as_tibble() %>%
    pivot_longer(
      cols = -sample,
      names_to = "sequence",
      values_to = "abundance"
    ) %>%
    filter(abundance > 0)
  
  #Total reads before chimera removal
  sample_totals <- tibble(
    sample = rownames(original_table),
    total_reads_before_chimera = as.numeric(
      rowSums(original_table)
    )
  )
  
  #Total reads after chimera removal
  nochim_totals <- tibble(
    sample = rownames(original_table),
    nonchim_reads = as.numeric(
      rowSums(nochim_table)
    )
  )
  
  #Candidate level statistics
  candidate_table <- removed_long %>%
    left_join(
      sample_totals,
      by = "sample"
    ) %>%
    group_by(sample) %>%
    mutate(
      total_chimera_reads = sum(abundance),
      
      relative_abundance = if_else(
        total_reads_before_chimera > 0,
        100 * abundance / total_reads_before_chimera,
        NA_real_
      ),
      
      share_of_chimera_reads = if_else(
        total_chimera_reads > 0,
        100 * abundance / total_chimera_reads,
        NA_real_
      )
    ) %>%
    ungroup() %>%
    mutate(
      denoising_method = denoising_method,
      sequence_length = nchar(sequence),
      
      candidate_flag =
        relative_abundance >= candidate_relative_threshold |
        share_of_chimera_reads >= candidate_chimera_share_threshold,
      
      candidate_flag_reason = case_when(
        relative_abundance >= candidate_relative_threshold &
          share_of_chimera_reads >= candidate_chimera_share_threshold ~
          "High sample abundance; dominates chimera-associated reads",
        
        relative_abundance >= candidate_relative_threshold ~
          "High relative abundance in sample",
        
        share_of_chimera_reads >= candidate_chimera_share_threshold ~
          "Dominates chimera-associated reads",
        
        TRUE ~ "Not flagged"
      )
    ) %>%
    arrange(
      sample,
      desc(abundance)
    )
  
  #Sample level QC summary
  sample_summary <- sample_totals %>%
    left_join(
      nochim_totals,
      by = "sample"
    ) %>%
    left_join(
      candidate_table %>%
        group_by(sample) %>%
        summarise(
          n_removed_ASVs_present = n_distinct(sequence),
          chimera_reads = sum(abundance),
          top_removed_sequence = sequence[
            which.max(abundance)
          ],
          top_removed_abundance = max(abundance),
          top_removed_relative_abundance = max(
            relative_abundance
          ),
          top_removed_share_of_chimera_reads =
            max(share_of_chimera_reads),
          .groups = "drop"
        ),
      by = "sample"
    ) %>%
    mutate(
      denoising_method = denoising_method,
      
      n_removed_ASVs_present = replace_na(
        n_removed_ASVs_present,
        0L
      ),
      
      chimera_reads = replace_na(
        chimera_reads,
        0
      ),
      
      top_removed_abundance = replace_na(
        top_removed_abundance,
        0L
      ),
      
      top_removed_relative_abundance = replace_na(
        top_removed_relative_abundance,
        0
      ),
      
      top_removed_share_of_chimera_reads = replace_na(
        top_removed_share_of_chimera_reads,
        0
      ),
      
      percentage_chimeric = if_else(
        total_reads_before_chimera > 0,
        100 * chimera_reads /
          total_reads_before_chimera,
        NA_real_
      ),
      
      sample_flag =
        percentage_chimeric >= sample_chimera_threshold |
        top_removed_relative_abundance >=
        candidate_relative_threshold,
      
      sample_flag_reason = case_when(
        percentage_chimeric >= sample_chimera_threshold &
          top_removed_relative_abundance >=
          candidate_relative_threshold ~
          "High total chimera loss; abundant removed ASV",
        
        percentage_chimeric >= sample_chimera_threshold ~
          "High total chimera loss",
        
        top_removed_relative_abundance >=
          candidate_relative_threshold ~
          "Abundant ASV removed as chimera",
        
        TRUE ~ "Not flagged"
      )
    ) %>%
    arrange(
      desc(sample_flag),
      desc(percentage_chimeric)
    )
  
  flagged_candidates <- candidate_table %>%
    filter(candidate_flag) %>%
    select(
      sample,
      denoising_method,
      sequence_length,
      abundance,
      relative_abundance,
      total_chimera_reads,
      share_of_chimera_reads,
      candidate_flag_reason,
      sequence
    )
  
  list(
    sample_summary = sample_summary,
    flagged_candidates = flagged_candidates,
    all_removed_candidates = candidate_table
  )
}

######locate filtered fastq files 
filtered_dir <- file.path(common_files, "filtered_reads")

if (!dir.exists(filtered_dir)) {
  stop("Filtered-read directory does not exist: ", filtered_dir)
}

filtered_reads_fwd_paths <- list.files(
  filtered_dir,
  pattern = "_R1_filtered\\.fastq\\.gz$",
  full.names = TRUE
)

filtered_reads_rev_paths <- list.files(
  filtered_dir,
  pattern = "_R2_filtered\\.fastq\\.gz$",
  full.names = TRUE
)

if (length(filtered_reads_fwd_paths) == 0) {
  stop("No filtered R1 files found in: ", filtered_dir)
}

if (length(filtered_reads_rev_paths) == 0) {
  stop("No filtered R2 files found in: ", filtered_dir)
}

#####pair filtered R1 and R2 files by sample
filtered_file_info <- tibble(
  file = c(
    filtered_reads_fwd_paths,
    filtered_reads_rev_paths
  ),
  read = c(
    rep("R1", length(filtered_reads_fwd_paths)),
    rep("R2", length(filtered_reads_rev_paths))
  )
) %>%
  mutate(
    filename = basename(file),
    sample = str_remove(
      filename,
      "_R[12]_filtered\\.fastq\\.gz$"
    )
  )
  #filtered pair check 
filtered_pair_check <- filtered_file_info %>%
  dplyr::count(sample, read, name = "n") %>%
  pivot_wider(
    names_from = read,
    values_from = n,
    values_fill = 0
  )

if (!all(c("R1", "R2") %in% names(filtered_pair_check)) ||
    any(filtered_pair_check$R1 != 1) ||
    any(filtered_pair_check$R2 != 1)) {
  stop("Each sample must have exactly one filtered R1 and R2 file.")
}

filtered_paired <- filtered_file_info %>%
  select(sample, read, file) %>%
  pivot_wider(
    names_from = read,
    values_from = file
  ) %>%
  arrange(sample)

sample.names <- filtered_paired$sample
filtered_reads_fwd_paths <- filtered_paired$R1
filtered_reads_rev_paths <- filtered_paired$R2

names(filtered_reads_fwd_paths) <- sample.names
names(filtered_reads_rev_paths) <- sample.names

######build error models for denoising 
  #avoid bulding error models if can load existing ones 
errF_path <- file.path(common_files, "errF.RDS")
errR_path <- file.path(common_files, "errR.RDS")

if (file.exists(errF_path) && file.exists(errR_path)) {
  
  message("Loading existing error models.")
  
  errF <- readRDS(errF_path)
  errR <- readRDS(errR_path)
  
} else {
  
  message("Building new error models.")
  
  #Learn fwr and rev read error models
  set.seed(123)
  
  errF <- learnErrors(
    filtered_reads_fwd_paths,
    nbases = 1e8,
    randomize = TRUE,
    multithread = FALSE,
    verbose = TRUE
  )
  
  errR <- learnErrors(
    filtered_reads_rev_paths,
    nbases = 1e8,
    randomize = TRUE,
    multithread = FALSE,
    verbose = TRUE
  )
  #save error models
  saveRDS(errF, errF_path)
  saveRDS(errR, errR_path)
}

file.info(c(errF_path, errR_path)) %>%
  select(size, mtime)

 #load error models
errF <- readRDS(file.path(common_files, "errF.RDS"))
errR <- readRDS(file.path(common_files, "errR.RDS"))

#plot error models
errF_plot <- plotErrors(
  errF,
  nominalQ = TRUE
) +
  ggtitle("R1 error model")

errR_plot <- plotErrors(
  errR,
  nominalQ = TRUE
) +
  ggtitle("R2 error model")

ggsave(
  file.path(my_plots, "errF_plot.pdf"),
  errF_plot,
  device = "pdf",
  width = 8,
  height = 6
)
ggsave(
  file.path(my_plots, "errR_plot.pdf"),
  errR_plot,
  device = "pdf",
  width = 8,
  height = 6
)

######define directory for denoising outputs 
outputs <- file.path(
  common_files,
  "denoising_outputs"
)

dir.create(
  outputs,
  recursive = TRUE,
  showWarnings = FALSE
)

#######dereplicate reads (dada can accept direct fastaq paths but next paired-end step with mergePairs needs dereplicated objects)
derepFs <- derepFastq(
  filtered_reads_fwd_paths,
  verbose = TRUE
)

derepRs <- derepFastq(
  filtered_reads_rev_paths,
  verbose = TRUE
)

  #Preserve sample identities
names(derepFs) <- sample.names
names(derepRs) <- sample.names

  #Save for later merging
saveRDS(
  derepFs,
  file.path(outputs, "derepFs.RDS")
)

saveRDS(
  derepRs,
  file.path(outputs, "derepRs.RDS")
)

######independent denoising (no pooling, faster)
time_F_nopool <- system.time({
  
  dadaF_nopool <- dada(
    derepFs,
    err = errF,
    pool = FALSE,
    multithread = FALSE,
    verbose = TRUE
  )
})
time_F_nopool

time_R_nopool <- system.time({
  
  dadaR_nopool <- dada(
    derepRs,
    err = errR,
    pool = FALSE,
    multithread = FALSE,
    verbose = TRUE
  )
})
time_R_nopool

saveRDS(
  dadaF_nopool,
  file.path(outputs, "dadaF_nopool.RDS")
)

saveRDS(
  dadaR_nopool,
  file.path(outputs, "dadaR_nopool.RDS")
)

#####pseudo-pooling denoising 
time_F_pseudopool <- system.time({
  
  dadaF_pseudopool <- dada(
    derepFs,
    err = errF,
    pool = "pseudo",
    multithread = FALSE,
    verbose = TRUE
  )
})
time_F_pseudopool


time_R_pseudopool <- system.time({
  
  dadaR_pseudopool <- dada(
    derepRs,
    err = errR,
    pool = "pseudo",
    multithread = FALSE,
    verbose = TRUE
  )
})
time_R_pseudopool

saveRDS(
  dadaF_pseudopool,
  file.path(outputs, "dadaF_pseudopool.RDS")
)

saveRDS(
  dadaR_pseudopool,
  file.path(outputs, "dadaR_pseudopool.RDS")
)

######Full pooling denoising (high cumputational cost)
time_F_pool <- system.time({
  
  dadaF_pool <- dada(
    derepFs,
    err = errF,
    pool = TRUE,
    multithread = FALSE,
    verbose = TRUE
  )
})
time_F_pool

saveRDS(
  dadaF_pool,
  file.path(outputs, "dadaF_pool.RDS")
)

time_R_pool <- system.time({
  
  dadaR_pool <- dada(
    derepRs,
    err = errR,
    pool = TRUE,
    multithread = FALSE,
    verbose = TRUE
  )
})
time_R_pool

saveRDS(
  dadaR_pool,
  file.path(outputs, "dadaR_pool.RDS")
)

######reload filtering table and clean samples name
filter_summary_df <- read_tsv(
  file.path(
    common_files,
    "filtering_retention.tsv"
  ),
  show_col_types = FALSE
)

filter_summary_df <- filter_summary_df %>%
  mutate(
    sample = basename(sample),
    sample = str_remove(
      sample,
      "_R1_primerfree\\.fastq\\.gz$"
    ),
    sample = str_remove(
      sample,
      "_R1_filtered\\.fastq\\.gz$"
    )
  )
#extract denoised counts no pool
denoised_counts <- tibble(
  sample = sample.names,
  
  denoised_F = unname(
    vapply(
      dadaF_nopool,
      getN,
      numeric(1)
    )
  ),
  
  denoised_R = unname(
    vapply(
      dadaR_nopool,
      getN,
      numeric(1)
    )
  )
)

denoising_nopool_summary <- filter_summary_df %>%
  select(
    sample,
    reads_in,
    reads_out,
    retained_percentage
  ) %>%
  rename(
    primerfree_reads = reads_in,
    filtered_reads = reads_out,
    filtering_retention_percent = retained_percentage
  ) %>%
  left_join(
    denoised_counts,
    by = "sample"
  ) %>%
  mutate(
    retained_F_percent =
      100 * denoised_F / filtered_reads,
    
    retained_R_percent =
      100 * denoised_R / filtered_reads
  )

#count inferred variants no pool 
variant_counts <- tibble(
  sample = sample.names,
  
  variants_F = unname(
    vapply(
      dadaF_nopool,
      count_variants,
      integer(1)
    )
  ),
  
  variants_R = unname(
    vapply(
      dadaR_nopool,
      count_variants,
      integer(1)
    )
  )
)

denoising_nopool_summary <- denoising_nopool_summary %>%
  left_join(
    variant_counts,
    by = "sample"
  )

print(
  denoising_nopool_summary,
  n = Inf
)

write_tsv(
  denoising_nopool_summary,
  file.path(
    outputs,
    "denoising_nopool_summary.tsv"
  )
)

#extract denoised counts pseudopool
denoised_counts <- tibble(
  sample = sample.names,
  
  denoised_F = unname(
    vapply(
      dadaF_pseudopool,
      getN,
      numeric(1)
    )
  ),
  
  denoised_R = unname(
    vapply(
      dadaR_pseudopool,
      getN,
      numeric(1)
    )
  )
)

denoising_pseudopool_summary <- filter_summary_df %>%
  select(
    sample,
    reads_in,
    reads_out,
    retained_percentage
  ) %>%
  rename(
    primerfree_reads = reads_in,
    filtered_reads = reads_out,
    filtering_retention_percent = retained_percentage
  ) %>%
  left_join(
    denoised_counts,
    by = "sample"
  ) %>%
  mutate(
    retained_F_percent =
      100 * denoised_F / filtered_reads,
    
    retained_R_percent =
      100 * denoised_R / filtered_reads
  )

#count inferred variants pseudopool
variant_counts <- tibble(
  sample = sample.names,
  
  variants_F = unname(
    vapply(
      dadaF_pseudopool,
      count_variants,
      integer(1)
    )
  ),
  
  variants_R = unname(
    vapply(
      dadaR_pseudopool,
      count_variants,
      integer(1)
    )
  )
)

denoising_pseudopool_summary <- denoising_pseudopool_summary %>%
  left_join(
    variant_counts,
    by = "sample"
  )

print(
  denoising_pseudopool_summary,
  n = Inf
)

write_tsv(
  denoising_pseudopool_summary,
  file.path(
    outputs,
    "denoising_pseudopool_summary.tsv"
  )
)

#extract denoised counts pool
denoised_counts <- tibble(
  sample = sample.names,
  
  denoised_F = unname(
    vapply(
      dadaF_pool,
      getN,
      numeric(1)
    )
  ),
  
  denoised_R = unname(
    vapply(
      dadaR_pool,
      getN,
      numeric(1)
    )
  )
)

denoising_pool_summary <- filter_summary_df %>%
  select(
    sample,
    reads_in,
    reads_out,
    retained_percentage
  ) %>%
  rename(
    primerfree_reads = reads_in,
    filtered_reads = reads_out,
    filtering_retention_percent = retained_percentage
  ) %>%
  left_join(
    denoised_counts,
    by = "sample"
  ) %>%
  mutate(
    retained_F_percent =
      100 * denoised_F / filtered_reads,
    
    retained_R_percent =
      100 * denoised_R / filtered_reads
  )

#count inferred variants pool
variant_counts <- tibble(
  sample = sample.names,
  
  variants_F = unname(
    vapply(
      dadaF_pool,
      count_variants,
      integer(1)
    )
  ),
  
  variants_R = unname(
    vapply(
      dadaR_pool,
      count_variants,
      integer(1)
    )
  )
)

denoising_pool_summary <- denoising_pool_summary %>%
  left_join(
    variant_counts,
    by = "sample"
  )

print(
  denoising_pool_summary,
  n = Inf
)

write_tsv(
  denoising_pool_summary,
  file.path(
    outputs,
    "denoising_pool_summary.tsv"
  )
)

######compared inferred reads and variants across pooling methods 
inference_comparison_long <- bind_rows(
  summarise_dada_results(
    dadaF_nopool,
    "Forward",
    "No pooling"
  ),
  summarise_dada_results(
    dadaF_pseudopool,
    "Forward",
    "Pseudo-pooling"
  ),
  #summarise_dada_results( dadaF_pool,"Forward","Full pooling"),
  
  summarise_dada_results(
    dadaR_nopool,
    "Reverse",
    "No pooling"
  ),
  summarise_dada_results(
    dadaR_pseudopool,
    "Reverse",
    "Pseudo-pooling"
  ),
  #summarise_dada_results(dadaR_pool,"Reverse","Full pooling")
) %>%
  mutate(
    pooling = factor(
      pooling,
      levels = c(
        "No pooling",
        "Pseudo-pooling",
        "Full pooling"
      )
    )
  ) %>%
  arrange(
    sample,
    direction,
    pooling
  )

print(
  inference_comparison_long,
  n = Inf
)

write_tsv(
  inference_comparison_long,
  file.path(
    outputs,
    "inference_pooling_comparison_long.tsv"
  )
)

######create denoising retention table 
denoising_retention_table <- inference_comparison_long %>%
  left_join(
    filter_summary_df %>%
      select(
        sample,
        filtered_reads = reads_out
      ),
    by = "sample"
  ) %>%
  mutate(
    retained_percent =
      100 * denoised_reads / filtered_reads,
    
    reads_removed =
      filtered_reads - denoised_reads,
    
    removed_percent =
      100 * reads_removed / filtered_reads
  ) %>%
  arrange(
    sample,
    direction,
    pooling
  )

print(
  denoising_retention_table,
  n = Inf
)

write_tsv(
  denoising_retention_table,
  file.path(
    outputs,
    "denoising_retention_by_method.tsv"
  )
)

######plot variant counts
variant_count_plot <- inference_comparison_long %>%
  ggplot(
    aes(
      x = pooling,
      y = inferred_variants,
      fill = pooling
    )
  ) +
  geom_col(
    width = 0.7
  ) +
  geom_text(
    aes(label = inferred_variants),
    vjust = -0.4,
    size = 3.3
  ) +
  facet_grid(
    rows = vars(direction),
    cols = vars(sample),
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "No pooling" = "#0072B2",
      "Pseudo-pooling" = "#009E73",
      "Full pooling" = "#D55E00"
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.15)
    )
  ) +
  labs(
    x = NULL,
    y = "Number of inferred variants",
    fill = "Inference method",
    title = "Inferred variants by pooling method"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(
  file.path(
    my_plots,
    "inferred_variants_pooling_comparison.pdf"
  ),
  plot = variant_count_plot,
  width = 11,
  height = 7
)

#####compare inferred variants sequence identity 
inferred_sequences <- bind_rows(
  extract_inferred_sequences(
    dadaF_nopool,
    "Forward",
    "nopool"
  ),
  extract_inferred_sequences(
    dadaF_pseudopool,
    "Forward",
    "pseudopool"
  ),
  #extract_inferred_sequences(dadaF_pool,"Forward","pool"),
  
  extract_inferred_sequences(
    dadaR_nopool,
    "Reverse",
    "nopool"
  ),
  extract_inferred_sequences(
    dadaR_pseudopool,
    "Reverse",
    "pseudopool"
  ),
  #extract_inferred_sequences(dadaR_pool,"Reverse","pool")
)

write_tsv(
  inferred_sequences,
  file.path(
    outputs,
    "inferred_sequences_by_method.tsv"
  )
)
#####calculate sequence overlap
variant_presence <- inferred_sequences %>%
  select(
    sample,
    direction,
    pooling,
    sequence
  ) %>%
  distinct() %>%
  mutate(
    present = TRUE
  ) %>%
  pivot_wider(
    names_from = pooling,
    values_from = present,
    values_fill = FALSE
  )

variant_overlap <- variant_presence %>%
  group_by(
    sample,
    direction
  ) %>%
  summarise(
    variants_nopool = sum(nopool),
    variants_pseudopool = sum(pseudopool),
    #variants_pool = sum(pool),
    
    shared_all = sum(
      nopool &
        pseudopool
    ), #& pool
    
    nopool_only = sum(
      nopool &
        !pseudopool
    ), #& !pool
    
    pseudopool_only = sum(
      !nopool &
        pseudopool
    ), #& !pool
    
    #pool_only = sum(!nopool & !pseudopool & pool),
    
    shared_nopool_pseudopool = sum(
      nopool &
        pseudopool
    ),
    
    #shared_nopool_pool = sum(nopool & pool),
    
    #shared_pseudopool_pool = sum(pseudopool &  pool),
    
    jaccard_nopool_pseudopool =
      safe_jaccard(
        nopool,
        pseudopool
      ),
    
    #jaccard_nopool_pool = safe_jaccard(nopool,pool),
    
    #jaccard_pseudopool_pool = safe_jaccard(pseudopool,pool),
    
    .groups = "drop"
  )

print(
  variant_overlap,
  n = Inf
)

write_tsv(
  variant_overlap,
  file.path(
    outputs,
    "inferred_variant_overlap.tsv"
  )
)

#####plot unique variants 
unique_variant_plot_data <- variant_overlap %>%
  select(
    sample,
    direction,
    nopool_only,
    pseudopool_only
  ) %>%
  pivot_longer(
    cols = c(
      nopool_only,
      pseudopool_only
    ),
    names_to = "method",
    values_to = "unique_variants"
  ) %>%
  mutate(
    method = recode(
      method,
      nopool_only = "No-pooling only",
      pseudopool_only = "Pseudo-pooling only"
    )
  )

unique_variant_plot <- ggplot(
  unique_variant_plot_data,
  aes(
    x = sample,
    y = unique_variants,
    fill = method
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  geom_text(
    aes(label = unique_variants),
    position = position_dodge(width = 0.75),
    vjust = -0.4,
    size = 3.5
  ) +
  facet_wrap(
    ~direction,
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "No-pooling only" = "#0072B2",
      "Pseudo-pooling only" = "#009E73"
    )
  ) +
  scale_y_continuous(
    breaks = scales::breaks_width(2),
    expand = expansion(
      mult = c(0, 0.15)
    )
  ) +
  labs(
    x = "Sample",
    y = "Method-specific variants",
    fill = "Category",
    title = "Variants unique to each inference method",
    subtitle = "Shared variants are excluded to highlight method differences"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(
  file.path(
    my_plots,
    "method_specific_inferred_variants.pdf"
  ),
  plot = unique_variant_plot,
  width = 9,
  height = 5
)

######plot Jaccard similarity (describs how similar the sets of DNA sequences inferred the pooling methods are)
#0: no inferred sequences are shared
#0.5: moderate overlap
#1: exactly the same set of inferred sequences
jaccard_plot <- ggplot(
  variant_overlap,
  aes(
    x = sample,
    y = jaccard_nopool_pseudopool,
    color = direction
  )
) +
  geom_point(
    position = position_dodge(width = 0.35),
    size = 3
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        jaccard_nopool_pseudopool
      )
    ),
    position = position_dodge(width = 0.35),
    vjust = -0.8,
    size = 3.3,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    limits = c(0.90, 1.00),
    breaks = seq(0.90, 1.00, 0.02)
  ) +
  scale_color_manual(
    values = c(
      "Forward" = "#0072B2",
      "Reverse" = "#D55E00"
    )
  ) +
  labs(
    x = "Sample",
    y = "Jaccard similarity",
    color = "Read direction",
    title = "Similarity of inferred variant sets",
    subtitle = "No-pooling versus pseudo-pooling"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(
  file.path(
    my_plots,
    "Jaccard_similarity_sequences_variants.pdf"
  ),
  plot = jaccard_plot,
  width = 9,
  height = 5
)

#######merge reads
#load objects 
derepFs <- readRDS(
  file.path(outputs, "derepFs.RDS")
)

derepRs <- readRDS(
  file.path(outputs, "derepRs.RDS")
)

dadaF_nopool <- readRDS(
  file.path(outputs, "dadaF_nopool.RDS")
)

dadaR_nopool <- readRDS(
  file.path(outputs, "dadaR_nopool.RDS")
)

dadaF_pseudopool <- readRDS(
  file.path(outputs, "dadaF_pseudopool.RDS")
)

dadaR_pseudopool <- readRDS(
  file.path(outputs, "dadaR_pseudopool.RDS")
)

#dadaF_pool <- readRDS(file.path(outputs, "dadaF_pool.RDS"))

#dadaR_pool <- readRDS(file.path(outputs, "dadaR_pool.RDS"))

######merge no-pooling reads
mergers_nopool <- mergePairs(
  dadaF = dadaF_nopool,
  derepF = derepFs,
  dadaR = dadaR_nopool,
  derepR = derepRs,
  minOverlap = 12,
  maxMismatch = 0,
  returnRejects = FALSE,
  verbose = TRUE
)

names(mergers_nopool) <- sample.names

saveRDS(
  mergers_nopool,
  file.path(outputs, "mergers_nopool.RDS")
)

#######merge pseudo-pooled reads
mergers_pseudopool <- mergePairs(
  dadaF = dadaF_pseudopool,
  derepF = derepFs,
  dadaR = dadaR_pseudopool,
  derepR = derepRs,
  minOverlap = 12,
  maxMismatch = 0,
  returnRejects = FALSE,
  verbose = TRUE
)

names(mergers_pseudopool) <- sample.names

saveRDS(
  mergers_pseudopool,
  file.path(outputs, "mergers_pseudopool.RDS")
)

######merge pooled reads 
mergers_pool <- mergePairs(
  dadaF = dadaF_pool,
  derepF = derepFs,
  dadaR = dadaR_pool,
  derepR = derepRs,
  minOverlap = 12,
  maxMismatch = 0,
  returnRejects = FALSE,
  verbose = TRUE
)

names(mergers_pool) <- sample.names

saveRDS(
  mergers_pool,
  file.path(outputs, "mergers_pool.RDS")
)

#####inspect merging 
merging_inspection <- tibble(
  sample = sample.names,
  
  nopool_merged_reads = vapply(
    mergers_nopool,
    function(x) sum(x$abundance),
    numeric(1)
  ),
  
  nopool_merged_sequences = vapply(
    mergers_nopool,
    function(x) length(unique(x$sequence)),
    integer(1)
  ),
  
  pseudopool_merged_reads = vapply(
    mergers_pseudopool,
    function(x) sum(x$abundance),
    numeric(1)
  ),
  
  pseudopool_merged_sequences = vapply(
    mergers_pseudopool,
    function(x) length(unique(x$sequence)),
    integer(1)
  )
)

print(merging_inspection, n = Inf)

#inspect overlap lengths 
overlap_summary <- purrr::imap_dfr(
  mergers_pseudopool,
  function(x, sample_name) {
    
    if (nrow(x) == 0) {
      return(
        tibble(
          sample = sample_name,
          min_overlap = NA_real_,
          median_overlap = NA_real_,
          max_overlap = NA_real_
        )
      )
    }
    
    tibble(
      sample = sample_name,
      min_overlap = min(x$nmatch),
      median_overlap = median(x$nmatch),
      max_overlap = max(x$nmatch)
    )
  }
)

print(overlap_summary)

#######after merging read tracking 
filter_lookup <- filter_summary_df %>%
  dplyr::transmute(
    sample = as.character(sample),
    primerfree_reads = reads_in,
    filtered_reads = reads_out
  ) %>%
  dplyr::distinct(
    sample,
    .keep_all = TRUE
  )

read_counts_merged <- dplyr::bind_rows(
  summarise_merging(
    dadaF = dadaF_nopool,
    dadaR = dadaR_nopool,
    mergers = mergers_nopool,
    method = "No pooling"
  ),
  
  summarise_merging(
    dadaF = dadaF_pseudopool,
    dadaR = dadaR_pseudopool,
    mergers = mergers_pseudopool,
    method = "Pseudo-pooling"
  )
) %>%
  dplyr::mutate(
    sample = as.character(sample),
    denoising_method =
      as.character(denoising_method)
  ) %>%
  dplyr::left_join(
    filter_lookup,
    by = "sample"
  ) %>%
  dplyr::mutate(
    perc_merged_from_filtered = dplyr::if_else(
      filtered_reads > 0,
      100 * merged_reads / filtered_reads,
      NA_real_
    )
  ) %>%
  dplyr::arrange(
    sample,
    denoising_method
  )
#check each sample-method combination occurs once 
tracking_duplicates <- read_counts_merged %>%
  dplyr::count(
    sample,
    denoising_method,
    name = "n"
  ) %>%
  dplyr::filter(n != 1)

if (nrow(tracking_duplicates) > 0) {
  print(tracking_duplicates)
  stop(
    "Duplicated or missing sample-method combinations ",
    "in read_counts_merged."
  )
}

print(
  read_counts_merged,
  n = Inf,
  width = Inf
)

write_tsv(
  read_counts_merged,
  file.path(
    common_files,
    "reads_counts_after_merging.tsv"
  )
)

#######chose type of plot based on number of samples 
density_sample_check <- read_counts_merged %>%
  dplyr::filter(
    !is.na(perc_merged),
    !is.na(denoising_method),
    !is.na(sample)
  ) %>%
  dplyr::group_by(denoising_method) %>%
  dplyr::summarise(
    n_samples = dplyr::n_distinct(sample),
    .groups = "drop"
  )

minimum_density_n <- 5

enough_for_density <- all(
  density_sample_check$n_samples >= minimum_density_n
)

#calculate adaptive limits for narrow percentage cluster 
merging_range <- range(
  read_counts_merged$perc_merged,
  na.rm = TRUE
)

x_lower <- max(
  0,
  floor(min(merging_range) / 5) * 5 - 5
)

x_upper <- min(
  100,
  ceiling(max(merging_range) / 5) * 5 + 5
)

if (enough_for_density) {
  
  merging_density_plot <- read_counts_merged %>%
    filter(
      !is.na(perc_merged),
      !is.na(denoising_method)
    ) %>%
    ggplot(
      aes(
        x = perc_merged,
        fill = denoising_method,
        color = denoising_method
      )
    ) +
    geom_density(
      alpha = 0.25,
      linewidth = 0.8,
      adjust = 1,
      trim = TRUE
    ) +
    geom_rug(
      aes(color = denoising_method),
      alpha = 0.6,
      sides = "b",
      show.legend = FALSE
    ) +
    scale_x_continuous(
      limits = c(x_lower, x_upper),
      labels = function(x) paste0(x, "%")
    ) +
    scale_fill_manual(
      values = c(
        "No pooling" = "#0072B2",
        "Pseudo-pooling" = "#009E73",
        "Full pooling" = "#D55E00"
      ),
      drop = TRUE
    ) +
    scale_color_manual(
      values = c(
        "No pooling" = "#0072B2",
        "Pseudo-pooling" = "#009E73",
        "Full pooling" = "#D55E00"
      ),
      drop = TRUE
    ) +
    labs(
      x = "Denoised reads successfully merged",
      y = "Density",
      fill = "Denoising method",
      color = "Denoising method",
      title = "Distribution of paired-read merging success",
      subtitle = "Minimum overlap = 12 bp; maximum mismatches = 0"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  final_merging_plot <- merging_density_plot
  plot_filename <- "merging_density_perc_reads.pdf"
  
} else {
  
  final_merging_plot <- read_counts_merged %>%
    ggplot(
      aes(
        x = denoising_method,
        y = perc_merged,
        group = sample,
        color = sample
      )
    ) +
    geom_line(
      linewidth = 0.7,
      alpha = 0.7
    ) +
    geom_point(
      size = 3
    ) +
    geom_text(
      aes(
        label = sprintf(
          "%.1f%%",
          perc_merged
        )
      ),
      vjust = -0.8,
      size = 3.2,
      show.legend = FALSE
    ) +
    scale_y_continuous(
      breaks = seq(x_lower, x_upper, by = 5),
      labels = function(x) paste0(x, "%")
    ) +
    coord_cartesian(
      ylim = c(x_lower, x_upper)
    ) +
    labs(
      x = "Denoising method",
      y = "Denoised reads successfully merged (%)",
      color = "Sample",
      title = "Paired-read merging success",
      subtitle = paste0(
        "Individual samples shown because n < ",
        minimum_density_n,
        " per method"
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      ),
      legend.position = "top"
    )
  
  plot_filename <- "merging_points_perc_reads.pdf"
}

ggsave(
  file.path(
    my_plots,
    plot_filename
  ),
  plot = final_merging_plot,
  device = "pdf",
  width = 7,
  height = 5
)

#####make ASV tables from merged reads
ASV_table_nopool <- makeSequenceTable(
  mergers_nopool
)

ASV_table_pseudopool <- makeSequenceTable(
  mergers_pseudopool
)

ASV_table_pool <- makeSequenceTable(
  mergers_pool
)
dim(ASV_table_nopool)#output is number of samples × number of ASVs
dim(ASV_table_pseudopool)
dim(ASV_table_pool)

#######check sequence length before chimera removal 
ASV_length_summary <- bind_rows(
  tibble(
    denoising_method = "No pooling",
    length = nchar(colnames(ASV_table_nopool))
  ),
  tibble(
    denoising_method = "Pseudo-pooling",
    length = nchar(colnames(ASV_table_pseudopool))
  )
  #,tibble(denoising_method = "Full pooling",length = nchar(colnames(ASV_table_pool)))
) %>%
  group_by(denoising_method) %>%
  summarise(
    n_ASVs = n(),
    minimum = min(length),
    median = median(length),
    maximum = max(length),
    .groups = "drop"
  )

print(ASV_length_summary)

saveRDS(
  ASV_table_nopool,
  file.path(common_files, "ASV_table_nopool.RDS")
)

saveRDS(
  ASV_table_pseudopool,
  file.path(common_files, "ASV_table_pseudopool.RDS")
)

saveRDS(
  ASV_table_pool,
  file.path(common_files, "ASV_table_pool.RDS")
)

########remove chimeras 
#method = "consensus" check each sample separatly and make consensus decision across samples with ASV 
ASV_table_nopool.nochim <- removeBimeraDenovo(
  ASV_table_nopool,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

ASV_table_pseudopool.nochim <- removeBimeraDenovo(
  ASV_table_pseudopool,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)
#method = "pooled" pool sample abundance for chimera detection 
ASV_table_pool.nochim <- removeBimeraDenovo(
  ASV_table_pool,
  method = "pooled",
  multithread = FALSE,
  verbose = TRUE
)

saveRDS(
  ASV_table_nopool.nochim,
  file.path(
    common_files,
    "ASV_table_nopool.nochim.RDS"
  )
)

saveRDS(
  ASV_table_pseudopool.nochim,
  file.path(
    common_files,
    "ASV_table_pseudopool.nochim.RDS"
  )
)

saveRDS(
  ASV_table_pool.nochim,
  file.path(
    common_files,
    "ASV_table_pool.nochim.RDS"
  )
)

######create summary tables 
chimera_summary <- tibble(
  denoising_method = c(
    "No pooling",
    "Pseudo-pooling"
  ),
  
  ASVs_before = c(
    ncol(ASV_table_nopool),
    ncol(ASV_table_pseudopool)
  ),
  
  ASVs_after = c(
    ncol(ASV_table_nopool.nochim),
    ncol(ASV_table_pseudopool.nochim)
  ),
  
  reads_before = c(
    sum(ASV_table_nopool),
    sum(ASV_table_pseudopool)
  ),
  
  reads_after = c(
    sum(ASV_table_nopool.nochim),
    sum(ASV_table_pseudopool.nochim)
  )
) %>%
  dplyr::mutate(
    chimeric_ASVs =
      ASVs_before - ASVs_after,
    
    ASVs_retained_percent =
      100 * ASVs_after / ASVs_before,
    
    ASVs_removed_percent =
      100 - ASVs_retained_percent,
    
    chimeric_reads =
      reads_before - reads_after,
    
    reads_retained_percent =
      100 * reads_after / reads_before,
    
    reads_removed_percent =
      100 - reads_retained_percent
  )

chimera_summary_view <- chimera_summary %>%
  dplyr::select(
    denoising_method,
    ASVs_before,
    ASVs_after,
    ASVs_retained_percent,
    reads_before,
    reads_after,
    reads_retained_percent
  )

print(
  chimera_summary_view,
  n = Inf
)

write_tsv(
  chimera_summary,
  file.path(
    common_files,
    "chimera_summary.tsv"
  )
)

######chimera removal QC 
chimera_qc_results <- list(
  
  nopool = flag_chimera_results(
    original_table = ASV_table_nopool,
    nochim_table = ASV_table_nopool.nochim,
    denoising_method = "No pooling",
    sample_chimera_threshold = 20,
    candidate_relative_threshold = 5,
    candidate_chimera_share_threshold = 50
  ),
  
  pseudopool = flag_chimera_results(
    original_table = ASV_table_pseudopool,
    nochim_table = ASV_table_pseudopool.nochim,
    denoising_method = "Pseudo-pooling",
    sample_chimera_threshold = 20,
    candidate_relative_threshold = 5,
    candidate_chimera_share_threshold = 50
  )
)

chimera_sample_method_QC <- purrr::map_dfr(
  chimera_qc_results,
  ~.x$sample_summary
)

flagged_chimera_candidates <- purrr::map_dfr(
  chimera_qc_results,
  ~.x$flagged_candidates
)

all_removed_chimera_candidates <- purrr::map_dfr(
  chimera_qc_results,
  ~.x$all_removed_candidates
)

#create table with flagged samples and inspect outputs 
flagged_sample_methods <- chimera_sample_method_QC %>%
  dplyr::filter(
    !is.na(sample_flag),
    sample_flag
  ) %>%
  dplyr::arrange(
    dplyr::desc(percentage_chimeric),
    dplyr::desc(top_removed_relative_abundance)
  )

print(chimera_sample_method_QC, n = Inf)

print(flagged_sample_methods, n = Inf)

flagged_samples <- chimera_sample_method_QC %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    methods_evaluated = dplyr::n_distinct(
      denoising_method
    ),
    
    methods_flagged = sum(
      sample_flag,
      na.rm = TRUE
    ),
    
    flagged_by_all_methods =
      methods_flagged == methods_evaluated,
    
    denoising_methods_flagged = paste(
      sort(
        unique(
          denoising_method[
            !is.na(sample_flag) &
              sample_flag
          ]
        )
      ),
      collapse = "; "
    ),
    
    minimum_percentage_chimeric = min(
      percentage_chimeric,
      na.rm = TRUE
    ),
    
    maximum_percentage_chimeric = max(
      percentage_chimeric,
      na.rm = TRUE
    ),
    
    maximum_removed_ASV_relative_abundance = max(
      top_removed_relative_abundance,
      na.rm = TRUE
    ),
    
    sample_flag = any(
      sample_flag,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(sample_flag) %>%
  dplyr::arrange(
    dplyr::desc(
      maximum_percentage_chimeric
    )
  )

flagged_chimera_candidates %>%
  select(
    sample,
    denoising_method,
    sequence_length,
    abundance,
    relative_abundance,
    share_of_chimera_reads,
    candidate_flag_reason,
    sequence
  ) %>%
  print(n = Inf)

flagged_candidate_agreement <-
  flagged_chimera_candidates %>%
  dplyr::group_by(
    sample,
    sequence
  ) %>%
  dplyr::summarise(
    sequence_length =
      nchar(sequence[1]),
    
    methods_flagged =
      dplyr::n_distinct(
        denoising_method
      ),
    
    denoising_methods = paste(
      sort(
        unique(denoising_method)
      ),
      collapse = "; "
    ),
    
    minimum_abundance =
      min(abundance),
    
    maximum_abundance =
      max(abundance),
    
    minimum_relative_abundance =
      min(relative_abundance),
    
    maximum_relative_abundance =
      max(relative_abundance),
    
    cross_method_flag =
      methods_flagged > 1,
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      cross_method_flag
    ),
    dplyr::desc(
      maximum_relative_abundance
    )
  )

print(
  flagged_candidate_agreement,
  n = Inf
)

chimera_qc_dir <- file.path(
  common_files,
  "chimera_QC"
)

dir.create(
  chimera_qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_tsv(
  chimera_sample_method_QC,
  file.path(
    chimera_qc_dir,
    "chimera_sample_method_QC.tsv"
  )
)

write_tsv(
  flagged_sample_methods,
  file.path(
    chimera_qc_dir,
    "flagged_sample_methods.tsv"
  )
)

write_tsv(
  flagged_samples,
  file.path(
    chimera_qc_dir,
    "flagged_samples.tsv"
  )
)

write_tsv(
  flagged_chimera_candidates,
  file.path(
    chimera_qc_dir,
    "flagged_candidates.tsv"
  )
)

write_tsv(
  flagged_candidate_agreement,
  file.path(
    chimera_qc_dir,
    "candidate_agreement.tsv"
  )
)

saveRDS(
  all_removed_chimera_candidates,
  file.path(
    chimera_qc_dir,
    "all_removed_chimera_candidates.RDS"
  ),
  compress = "xz"
)

#####read tracking after chimera removal 
chimera_stage_counts <- dplyr::bind_rows(
  tibble(
    sample = rownames(ASV_table_nopool),
    denoising_method = "No pooling",
    
    merged_table_reads = as.numeric(
      rowSums(ASV_table_nopool)
    ),
    
    ASVs_before_chimera = as.integer(
      rowSums(ASV_table_nopool > 0)
    ),
    nonchimeric_reads = as.numeric(
      rowSums(ASV_table_nopool.nochim)
    ),
    
    nonchimeric_ASVs = as.integer(
      rowSums(ASV_table_nopool.nochim > 0)
    )
  ),
  
  tibble(
    sample = rownames(
      ASV_table_pseudopool
    ),
    
    denoising_method = "Pseudo-pooling",
    
    merged_table_reads = as.numeric(
      rowSums(ASV_table_pseudopool)
    ),
    
    ASVs_before_chimera = as.integer(
      rowSums(ASV_table_pseudopool > 0)
    ),
    nonchimeric_reads = as.numeric(
      rowSums(
        ASV_table_pseudopool.nochim
      )
    ),
    
    nonchimeric_ASVs = as.integer(
      rowSums(
        ASV_table_pseudopool.nochim > 0
      )
    )
  )
)

######add after chimera removal read track to final count read tracking table 
read_counts_final <- read_counts_merged %>%
  dplyr::left_join(
    chimera_stage_counts,
    by = c(
      "sample",
      "denoising_method"
    ),
    relationship = "one-to-one"
  ) %>%
  dplyr::mutate(
    chimera_reads =
      merged_table_reads -
      nonchimeric_reads,
    
    perc_nonchimeric_from_merged =
      dplyr::if_else(
        merged_table_reads > 0,
        100 * nonchimeric_reads /
          merged_table_reads,
        NA_real_
      ),
    
    perc_chimeric_from_merged =
      dplyr::if_else(
        merged_table_reads > 0,
        100 * chimera_reads /
          merged_table_reads,
        NA_real_
      ),
    
    perc_final_from_filtered =
      dplyr::if_else(
        filtered_reads > 0,
        100 * nonchimeric_reads /
          filtered_reads,
        NA_real_
      ),
    
    perc_final_from_primerfree =
      dplyr::if_else(
        primerfree_reads > 0,
        100 * nonchimeric_reads /
          primerfree_reads,
        NA_real_
      )
  ) %>%
  dplyr::arrange(
    sample,
    denoising_method
  )

tracking_mismatch <- read_counts_final %>%
  dplyr::filter(
    merged_reads !=
      merged_table_reads
  )

if (nrow(tracking_mismatch) > 0) {
  print(tracking_mismatch)
  stop(
    "Merged-read counts do not match ",
    "the ASV-table row sums."
  )
}

if (anyNA(
  read_counts_final$nonchimeric_reads
)) {
  stop(
    "Missing non-chimeric counts after joining."
  )
}

read_counts_final %>%
  dplyr::select(
    sample,
    denoising_method,
    primerfree_reads,
    filtered_reads,
    denoised_F,
    denoised_R,
    merged_reads,
    perc_merged,
    nonchimeric_reads,
    perc_nonchimeric_from_merged,
    perc_final_from_filtered,
    nonchimeric_ASVs
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

write_tsv(
  read_counts_final,
  file.path(
    common_files,
    "reads_counts_after_denoising_and_chimera_removal.tsv"
  )
)

#####compare read count at each stage 
primer_detection_raw <- read_tsv(
  file.path(
    common_files,
    "primer_detection_long.tsv"
  ),
  show_col_types = FALSE
)

raw_pair_counts <- primer_detection_raw %>%
  dplyr::distinct(
    sample,
    read,
    n_reads
  ) %>%
  tidyr::pivot_wider(
    names_from = read,
    values_from = n_reads
  )

if (anyNA(raw_pair_counts$R1) ||
    anyNA(raw_pair_counts$R2)) {
  stop("Missing raw R1 or R2 count.")
}

if (any(
  raw_pair_counts$R1 !=
  raw_pair_counts$R2
)) {
  stop("Raw R1 and R2 counts differ.")
}

raw_read_lookup <- raw_pair_counts %>%
  dplyr::transmute(
    sample,
    raw_reads = R1
  )
read_counts_final <- read_counts_final %>%
  dplyr::select(
    -dplyr::any_of("raw_reads")
  ) %>%
  dplyr::left_join(
    raw_read_lookup,
    by = "sample",
    relationship = "many-to-one"
  )
required_tracking_columns <- c(
  "sample",
  "denoising_method",
  "raw_reads",
  "primerfree_reads",
  "filtered_reads",
  "pairable_denoised_reads",
  "merged_reads",
  "nonchimeric_reads"
)

missing_tracking_columns <- setdiff(
  required_tracking_columns,
  names(read_counts_final)
)

if (length(missing_tracking_columns) > 0) {
  stop(
    "Missing tracking columns: ",
    paste(
      missing_tracking_columns,
      collapse = ", "
    )
  )
}
tracking_stages <- c(
  "Raw",
  "Primers removed",
  "Filtered",
  "Denoised",
  "Merged",
  "Chimera removal"
)

tracking_stages <- c(
    "Raw",
    "Primers removed",
    "Filtered",
    "Denoised",
    "Merged",
    "Chimera removal"
  )

read_tracking_long <- read_counts_final %>%
  dplyr::select(
    sample,
    denoising_method,
    raw_reads,
    primerfree_reads,
    filtered_reads,
    pairable_denoised_reads,
    merged_reads,
    nonchimeric_reads
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      raw_reads,
      primerfree_reads,
      filtered_reads,
      pairable_denoised_reads,
      merged_reads,
      nonchimeric_reads
    ),
    names_to = "stage",
    values_to = "n_reads"
  ) %>%
  dplyr::mutate(
    stage = dplyr::recode(
      stage,
      raw_reads = "Raw",
      primerfree_reads = "Primers removed",
      filtered_reads = "Filtered",
      pairable_denoised_reads = "Denoised",
      merged_reads = "Merged",
      nonchimeric_reads = "Chimera removal"
    ),
    
    stage = factor(
      stage,
      levels = tracking_stages,
      ordered = TRUE
    ),
    
    denoising_method = factor(
      denoising_method,
      levels = c(
        "No pooling",
        "Pseudo-pooling"
      ),
      ordered = TRUE
    )
  )

#####plot read count for each stage 
final_read_stats_plot <- ggplot(
  read_tracking_long,
  aes(
    x = stage,
    y = n_reads,
    group = sample,
    color = sample
  )
) +
  geom_line(
    alpha = 0.7,
    linewidth = 0.7
  ) +
  geom_point(
    alpha = 0.8,
    size = 2.3
  ) +
  facet_wrap(
    ~denoising_method,
    nrow = 1
  ) +
  scale_y_continuous(
    labels = scales::label_comma()
  ) +
  labs(
    x = NULL,
    y = "Number of reads",
    color = "Sample",
    title = "Read retention through the sequencing pipeline"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    my_plots,
    "final_read_stats_plot.pdf"
  ),
  plot = final_read_stats_plot,
  device = "pdf",
  width = 9,
  height = 5
)

######create relative retention plot (percentage plot)
read_tracking_relative <- read_tracking_long %>%
  dplyr::group_by(
    sample,
    denoising_method
  ) %>%
  dplyr::mutate(
    raw_total = n_reads[
      stage == "Raw"
    ][1],
    
    retention_percentage =
      dplyr::if_else(
        raw_total > 0,
        100 * n_reads / raw_total,
        NA_real_
      )
  ) %>%
  dplyr::ungroup()
read_retention_plot <- ggplot(
  read_tracking_relative,
  aes(
    x = stage,
    y = retention_percentage,
    group = sample,
    color = sample
  )
) +
  geom_line(
    alpha = 0.7,
    linewidth = 0.7
  ) +
  geom_point(size = 2.3) +
  facet_wrap(
    ~denoising_method,
    nrow = 1
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    labels = function(x) {
      paste0(x, "%")
    }
  ) +
  labs(
    x = NULL,
    y = "Reads retained from raw input",
    color = "Sample",
    title = "Relative read retention through the pipeline"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    my_plots,
    "final_read_retention_plot.pdf"
  ),
  plot = read_retention_plot,
  width = 9,
  height = 5
)
write_tsv(
  read_tracking_long,
  file.path(
    common_files,
    "read_tracking_long.tsv"
  )
)

########define comparison objects (list) for comparisons 
comparison_tables <- list(
  nopool =
    ASV_table_nopool.nochim,
  pseudopool =
    ASV_table_pseudopool.nochim
)

method_labels <- c(
  nopool = "No pooling",
  pseudopool = "Pseudo-pooling"
)

comparison_methods <- names(
  comparison_tables
)

######calculate number of overlapping sequences in ASV
ASV_sets <- lapply(
  comparison_tables,
  colnames
)

total_unique_ASVs <- length(
  Reduce(
    union,
    ASV_sets
  )
)

total_unique_ASVs

ASV_overlap_summary <- tibble(
  shared_ASVs = length(
    intersect(
      ASV_sets$nopool,
      ASV_sets$pseudopool
    )
  ),
  
  nopool_only = length(
    setdiff(
      ASV_sets$nopool,
      ASV_sets$pseudopool
    )
  ),
  
  pseudopool_only = length(
    setdiff(
      ASV_sets$pseudopool,
      ASV_sets$nopool
    )
  ),
  
  union_ASVs = total_unique_ASVs
) %>%
  dplyr::mutate(
    jaccard_similarity =
      shared_ASVs / union_ASVs
  )

write_tsv(
  ASV_overlap_summary,
  file.path(
    common_files,
    "ASV_method_overlap.tsv"
  )
)

#######Venn diagram (for no pooling and pseudo-polling methods, need to be adapted to integrate pooling)
VennDiagram::venn.diagram(
  x = list(
    "No pooling" =
      ASV_sets$nopool,
    "Pseudo-pooling" =
      ASV_sets$pseudopool
  ),
  filename = file.path(
    my_plots,
    "ASV_pooling_comparison.png"
  ),
  fill = c(
    "#0072B2",
    "#009E73"
  ),
  alpha = c(
    0.55,
    0.55
  ),
  col = "white",
  cat.col = c(
    "#0072B2",
    "#009E73"
  ),
  category.names = c(
    "No pooling",
    "Pseudo-pooling"
  ),
  output = TRUE
)

#####calculate and plot ASV richness per samples and denoising methods 
#samples are ordered by their mean richness across the 2 denoising methods used
asv_richness <- purrr::imap_dfr(
  comparison_tables,
  function(table, method) {
    
    tibble(
      sample = rownames(table),
      denoising_method =
        method_labels[[method]],
      observed_ASVs = as.integer(
        rowSums(table > 0)
      )
    )
  }
)
richness_order <- asv_richness %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    mean_richness =
      mean(observed_ASVs),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_richness) %>%
  dplyr::mutate(
    sample_rank =
      dplyr::row_number()
  )
asv_richness <- asv_richness %>%
  dplyr::left_join(
    richness_order,
    by = "sample"
  ) %>%
  dplyr::mutate(
    denoising_method = factor(
      denoising_method,
      levels = c(
        "No pooling",
        "Pseudo-pooling"
      )
    )
  )
asv_richness_plot <- ggplot(
  asv_richness,
  aes(
    x = sample_rank,
    y = observed_ASVs,
    color = denoising_method,
    group = denoising_method
  )
) +
  geom_line(
    linewidth = 0.7,
    alpha = 0.7
  ) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "No pooling" = "#0072B2",
      "Pseudo-pooling" = "#009E73"
    )
  ) +
  scale_x_continuous(
    breaks = richness_order$sample_rank,
    labels = richness_order$sample
  ) +
  labs(
    x = "Sample",
    y = "Observed non-chimeric ASVs",
    color = "Denoising method",
    title = "Sample richness by denoising method"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    my_plots,
    "asv_richness_plot.pdf"
  ),
  plot = asv_richness_plot,
  width = 7,
  height = 5
)

write_tsv(
  asv_richness,
  file.path(
    common_files,
    "ASV_richness_by_method.tsv"
  )
)

########create ASV abundance comparison plot per sample and method 
ASV_sums <- purrr::imap_dfr(
  comparison_tables,
  function(table, method) {
    
    tibble(
      ASV = colnames(table),
      denoising_method = method,
      total_abundance =
        as.numeric(
          colSums(table)
        )
    )
  }
)
ASV_sums_wide <- ASV_sums %>%
  tidyr::pivot_wider(
    names_from = denoising_method,
    values_from = total_abundance,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    shared_by = dplyr::case_when(
      nopool > 0 &
        pseudopool > 0 ~
        "Both methods",
      
      nopool > 0 ~
        "No pooling only",
      
      pseudopool > 0 ~
        "Pseudo-pooling only",
      
      TRUE ~
        "Neither"
    )
  )
ASV_sums_long <- ASV_sums_wide %>%
  tidyr::pivot_longer(
    cols = c(
      nopool,
      pseudopool
    ),
    names_to = "method_code",
    values_to = "n_reads"
  ) %>%
  dplyr::filter(
    n_reads > 0
  ) %>%
  dplyr::mutate(
    denoising_method =
      dplyr::recode(
        method_code,
        nopool = "No pooling",
        pseudopool =
          "Pseudo-pooling"
      ),
    
    denoising_method = factor(
      denoising_method,
      levels = c(
        "No pooling",
        "Pseudo-pooling"
      )
    ),
    
    shared_by = factor(
      shared_by,
      levels = c(
        "No pooling only",
        "Both methods",
        "Pseudo-pooling only"
      )
    )
  )
#plot abundance distribution 
ASV_abundance_palette <- c(
  "No pooling only" = "#0072B2",
  "Both methods" = "#CC79A7",
  "Pseudo-pooling only" = "#009E73"
)

denoising_ASV_comparison_plot <- ggplot(
  ASV_sums_long,
  aes(
    x = n_reads,
    fill = shared_by
  )
) +
  geom_histogram(
    bins = 30,
    position = "identity",
    alpha = 0.8
  ) +
  scale_x_log10(
    breaks = 10^(0:6),
    labels = scales::label_number()
  ) +
  scale_fill_manual(
    values =
      ASV_abundance_palette,
    drop = TRUE
  ) +
  facet_wrap(
    ~denoising_method,
    nrow = 1
  ) +
  labs(
    x = "Total ASV abundance",
    y = "Number of ASVs",
    fill = "Detected by",
    title = "ASV abundance by denoising method"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "top",
    panel.grid.minor =
      element_blank()
  )

ggsave(
  filename = file.path(
    my_plots,
    "denoising_ASV_comparison_plot.pdf"
  ),
  plot =
    denoising_ASV_comparison_plot,
  width = 8,
  height = 5
)
write_tsv(
  ASV_sums_wide,
  file.path(
    common_files,
    "ASV_abundance_method_comparison.tsv"
  )
)

######create chimera by sample table 
chimera_abundance <- chimera_sample_method_QC %>%
  dplyr::transmute(
    sample,
    denoising_method,
    chimera_rel_abundance = percentage_chimeric
  ) %>%
  dplyr::arrange(
    sample,
    denoising_method
  )

print(chimera_abundance, n = Inf)

write_tsv(
  chimera_abundance,
  file.path(
    common_files,
    "chimera_abundance_by_sample.tsv"
  )
)

#######check sequences lenght before and after chimera removal and plot 
ASV_length_status <- bind_rows(
  tibble(
    method = "No pooling",
    sequence = colnames(ASV_table_nopool),
    length = nchar(sequence),
    status = if_else(
      sequence %in% colnames(ASV_table_nopool.nochim),
      "Retained",
      "Removed as chimera"
    )
  ),
  
  tibble(
    method = "Pseudo-pooling",
    sequence = colnames(ASV_table_pseudopool),
    length = nchar(sequence),
    status = if_else(
      sequence %in% colnames(ASV_table_pseudopool.nochim),
      "Retained",
      "Removed as chimera"
    )
  )
) %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "No pooling",
        "Pseudo-pooling"
      )
    ),
    status = factor(
      status,
      levels = c(
        "Retained",
        "Removed as chimera"
      )
    )
  )

ASV_length_status %>%
  dplyr::group_by(method, status) %>%
  dplyr::summarise(
    n_ASVs = dplyr::n(),
    minimum = min(length),
    q1 = quantile(length, 0.25),
    median = median(length),
    q3 = quantile(length, 0.75),
    maximum = max(length),
    .groups = "drop"
  )

ASV_length_counts <- ASV_length_status %>%
  dplyr::count(
    method,
    status,
    length,
    name = "n_ASVs"
  )

ASV_length_barplot <- ggplot(
  ASV_length_counts,
  aes(
    x = length,
    y = n_ASVs
  )
) +
  geom_col(
    fill = "#4477AA",
    width = 0.85
  ) +
  facet_grid(
    rows = vars(status),
    cols = vars(method),
    scales = "free_y"
  ) +
  labs(
    x = "Merged ASV length (bp)",
    y = "Number of ASVs",
    title = "Length distribution of retained and chimeric ASVs",
    subtitle = "Each sequence is classified as retained or removed"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = "grey92",
      color = "grey60"
    )
  )

ggsave(
  filename = file.path(
    my_plots,
    "barplot_ASV_length_by_chimera_status.pdf"
  ),
  plot = ASV_length_barplot,
  device = "pdf",
  width = 9,
  height = 6
)

#######create adaptive plot (depends on number of samples)
chimera_sample_check <- chimera_abundance %>%
  filter(
    !is.na(chimera_rel_abundance)
  ) %>%
  group_by(denoising_method) %>%
  summarise(
    n_samples = n_distinct(sample),
    .groups = "drop"
  )

minimum_density_n <- 5

enough_for_density <-
  nrow(chimera_sample_check) > 0 &&
  all(
    chimera_sample_check$n_samples >=
      minimum_density_n
  )
if (enough_for_density) {
  
  chimera_plot <- ggplot(
    chimera_abundance,
    aes(
      x = chimera_rel_abundance,
      fill = denoising_method,
      color = denoising_method
    )
  ) +
    geom_density(
      alpha = 0.25,
      linewidth = 0.8,
      trim = TRUE
    ) +
    geom_rug(
      alpha = 0.6,
      show.legend = FALSE
    ) +
    labs(
      x = "Chimeric reads (%)",
      y = "Density",
      fill = "Denoising method",
      color = "Denoising method",
      title = "Distribution of chimeric-read abundance"
    )
  
  chimera_plot_filename <-
    "chimera_abundance_density.pdf"
  
} else {
  
  chimera_plot <- ggplot(
    chimera_abundance,
    aes(
      x = denoising_method,
      y = chimera_rel_abundance,
      group = sample,
      color = sample
    )
  ) +
    geom_line(
      linewidth = 0.7,
      alpha = 0.7
    ) +
    geom_point(size = 3) +
    geom_text(
      aes(
        label = sprintf(
          "%.2f%%",
          chimera_rel_abundance
        )
      ),
      vjust = -0.8,
      size = 3.1,
      show.legend = FALSE
    ) +
    labs(
      x = "Denoising method",
      y = "Chimeric reads (%)",
      color = "Sample",
      title = "Chimeric-read abundance by sample"
    )
  
  chimera_plot_filename <-
    "chimera_abundance_points.pdf"
}

chimera_plot <- chimera_plot +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "top"
  )

ggsave(
  file.path(
    my_plots,
    chimera_plot_filename
  ),
  plot = chimera_plot,
  device = "pdf",
  width = 7,
  height = 5
)

######save threshold metadata (chimera analysis parameters)
chimera_QC_parameters <- tibble(
  parameter = c(
    "sample_chimera_threshold",
    "candidate_relative_threshold",
    "candidate_chimera_share_threshold"
  ),
  value = c(
    20,
    5,
    50
  ),
  unit = c(
    "percent of sample reads",
    "percent of sample reads",
    "percent of chimera-associated reads"
  )
)

write_tsv(
  chimera_QC_parameters,
  file.path(
    chimera_qc_dir,
    "chimera_QC_parameters.tsv"
  )
)

#######optional investigation of the highest priority candidate 
priority_candidates <- flagged_chimera_candidates %>%
  dplyr::group_by(
    sample,
    sequence
  ) %>%
  dplyr::summarise(
    sequence_length =
      nchar(sequence[1]),
    
    methods_flagged =
      dplyr::n_distinct(denoising_method),
    
    denoising_methods = paste(
      sort(unique(denoising_method)),
      collapse = "; "
    ),
    
    maximum_abundance =
      max(abundance),
    
    maximum_relative_abundance =
      max(relative_abundance),
    
    maximum_share_of_chimera_reads =
      max(share_of_chimera_reads),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(maximum_relative_abundance),
    dplyr::desc(maximum_share_of_chimera_reads)
  )

if (nrow(priority_candidates) > 0) {
  
  priority_candidate <-
    priority_candidates %>%
    dplyr::slice(1)
  
  priority_sample <-
    priority_candidate$sample
  
  priority_sequence <-
    priority_candidate$sequence
  
  print(priority_candidate)
  
} else {
  
  message(
    "No chimera candidates passed the review thresholds."
  )
}

