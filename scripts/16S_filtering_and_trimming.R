#setup work environment
if (!requireNamespace("here", quietly = TRUE)) {
  BiocManager::install("here")
}
library(here)

#install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("dada2", quietly = TRUE)) {
  BiocManager::install("dada2")
}

#load required libraries
rm(list=ls()) 
library(dada2)
library(ShortRead)
library(Biostrings)
library(tidyverse)

#define paths, core variables and create directories
my_root <- "."
my_plots <- here::here("results", "plots", "filtering_and_trimming")
common_files <- here::here("data")

dir.create(my_plots, showWarnings = F, recursive = T)

# define primer sequences
fwd.primer <- "CCTACGGGNGGCWGCAG"
rev.primer <- "GACTACHVGGGTATCTAATCC"
fwd.primer.rev <- as.character(reverseComplement(DNAStringSet(fwd.primer)))
rev.primer.rev <- as.character(reverseComplement(DNAStringSet(rev.primer)))

#define directory with fastq files
fastq_files <- list.files(
  common_files,
  pattern = "\\.(fastq|fq)(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(fastq_files) == 0) {
  stop("No FASTQ files found in: ", common_files)
}
#recognize file names containing _R1_, _R2_, _1. or _2.
file_info <- tibble(
  file = fastq_files,
  filename = basename(fastq_files)
) %>%
  mutate(
    read = case_when(
      str_detect(
        filename,
        regex("_R?1(?:_|\\.)", ignore_case = TRUE)
      ) ~ "R1",
      str_detect(
        filename,
        regex("_R?2(?:_|\\.)", ignore_case = TRUE)
      ) ~ "R2",
      TRUE ~ NA_character_
    ),
  
#Keeps everything before _R1, _R2, _1 or _2 as the sample ID
sample = str_remove(
  filename,
  regex(
    "_R?[12](?:_|\\.).*$",
    ignore_case = TRUE
  )
)
  )

if (anyNA(file_info$read)) {
  print(file_info %>% filter(is.na(read)))
  stop("Some FASTQ filenames could not be classified as R1 or R2. Edit the filename patterns.")
}

# Check that every sample has one R1 and one R2 file
pair_check <- file_info %>%
  count(sample, read) %>%
  tidyr::pivot_wider(
    names_from = read,
    values_from = n,
    values_fill = 0
  )

print(pair_check)

if (!all(c("R1", "R2") %in% names(pair_check)) ||
    any(pair_check$R1 != 1) ||
    any(pair_check$R2 != 1)) {
  stop("Each sample must have exactly one R1 file and one R2 file.")
}

# Define the two biological primers in the orientation expected per read
primer_map <- tribble(
  ~read, ~primer,   ~sequence,
  "R1",  "Forward", fwd.primer,
  "R1",  "Reverse", rev.primer.rev,
  "R2",  "Forward", fwd.primer.rev,
  "R2",  "Reverse", rev.primer
)

# Count primer-containing reads
count_primers_in_fastq <- function(file, read_type, chunk_size = 100000) {
  
  patterns <- primer_map %>%
    filter(read == read_type)
  
  streamer <- FastqStreamer(file, n = chunk_size)
  on.exit(close(streamer))
  
  total_reads <- 0
  #Separate counters for primer at start of read and primer found anywhere in the read (primer position check for trimming optimisation)
  hit_counts_anywhere <- numeric(nrow(patterns))
  hit_counts_start <- numeric(nrow(patterns))
  
  repeat {
    fq_chunk <- yield(streamer)
    
    if (length(fq_chunk) == 0) {
      break
    }
    
    sequences <- sread(fq_chunk)
    total_reads <- total_reads + length(sequences)
    
    for (j in seq_len(nrow(patterns))) {
      
      primer_sequence <- patterns$sequence[j]
      primer_length <- nchar(primer_sequence)
      #look for complete reads
      matches_anywhere <- vcountPattern(
        pattern = DNAString(primer_sequence),
        subject = sequences,
        max.mismatch = 0,
        fixed = FALSE
      )
      # Count reads with at least one occurrence, not total occurrences
      hit_counts_anywhere[j] <- hit_counts_anywhere[j] +
        sum(matches_anywhere > 0)
      #look for only first primer_lenght bases
      long_enough <- width(sequences) >= primer_length
      
      if (any(long_enough)) {
        
        first_bases <- narrow(
          sequences[long_enough],
          start = 1,
          width = primer_length
        )
        
        matches_start <- vcountPattern(
          pattern = DNAString(primer_sequence),
          subject = first_bases,
          max.mismatch = 0,
          fixed = FALSE
        )
        # Count reads with at least one occurrence, not total occurrences
        hit_counts_start[j] <- hit_counts_start[j] +
          sum(matches_start > 0)
      }
    }
  }
  
# Calculate hit percentage per primer (anywhere in the read and at the start)
  if (total_reads > 0) {
    percentage_anywhere <- 100 * hit_counts_anywhere / total_reads
    percentage_start <- 100 * hit_counts_start / total_reads
  } else {
    percentage_anywhere <- rep(NA_real_, nrow(patterns))
    percentage_start <- rep(NA_real_, nrow(patterns))
  }
  
  patterns %>%
    transmute(
      primer,
      sequence,
      n_reads = total_reads,
      n_hits_anywhere = hit_counts_anywhere,
      percentage_anywhere = percentage_anywhere,
      n_hits_start = hit_counts_start,
      percentage_start = percentage_start
    )
}
      
# Run the counter once for each FASTQ file and generate table "primer_detection_long.tsv" to look for both primers in both reads(final table should contain 12 rows: 3 samples x 2 read directions x 2 primers)
primer_detection <- file_info %>%
  mutate(
    detection = map2(
      file,
      read,
      ~count_primers_in_fastq(.x, .y)
    )
  ) %>%
  select(sample, read, filename, detection) %>%
  unnest(detection) %>%
  arrange(sample, read, primer)

print(primer_detection, n = Inf)

write_tsv(
  primer_detection,
  file.path(common_files, "primer_detection_long.tsv")
)

#detect and plot primers
primer_detection <- read_tsv(file.path(common_files, "primer_detection_long.tsv"), show_col_types = F)

plot_primer_detection <- primer_detection %>%
  mutate(
    sample = factor(sample),
    primer = factor(primer, levels = c("Forward", "Reverse"))
  ) %>%
  ggplot(
    aes(
      x = sample,
      y = percentage_anywhere,
      color = primer
    )
  ) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_text(
    aes(label = sprintf("%.3g%%", percentage_anywhere)),
    position = position_dodge(width = 0.4),
    vjust = -0.8,
    size = 3,
    show.legend = FALSE
  ) +
  facet_wrap(~read, nrow = 1) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(
      sigma = 0.01,
      base = 10
    ),
    breaks = c(0, 0.01, 0.1, 1, 10, 100),
    limits = c(0, 110)
  ) +
  scale_color_manual(
    values = c(
      "Forward" = "#0072B2",
      "Reverse" = "#D55E00"
    )
  ) +
  labs(
    x = "Sample",
    y = "Reads containing primer (%)",
    color = "Primer",
    title = "Primer detection in paired-end reads",
    subtitle = "Exact matches found anywhere in each read"
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

plot_primer_detection

ggsave(
  file.path(my_plots, "primer_detection_dotplot.pdf"),
  plot_primer_detection,
  width = 8,
  height = 4.5
)

primer_position_difference <- primer_detection %>%
  mutate(
    expected = case_when(
      read == "R1" & primer == "Forward" ~ TRUE,
      read == "R2" & primer == "Reverse" ~ TRUE,
      TRUE ~ FALSE
    ),
    percentage_not_at_start =
      percentage_anywhere - percentage_start
  ) %>%
  filter(expected)

plot_primer_difference <- primer_position_difference %>%
  ggplot(
    aes(
      x = sample,
      y = percentage_not_at_start,
      fill = read
    )
  ) +
  geom_col(
    width = 0.65,
    position = position_dodge(width = 0.7)
  ) +
  geom_text(
    aes(
      label = sprintf("%.3f%%", percentage_not_at_start)
    ),
    vjust = -0.4,
    size = 3.3
  ) +
  facet_wrap(~read, nrow = 1) +
  scale_fill_manual(
    values = c(
      "R1" = "#0072B2",
      "R2" = "#D55E00"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = "Sample",
    y = "Expected primer detected away from position 1 (%)",
    fill = "Read",
    title = "Expected primers not anchored at the 5' start",
    subtitle = "Difference between whole-read and position-1 detection"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

plot_primer_difference

ggsave(
  file.path(my_plots, "primer_position_difference.pdf"),
  plot_primer_difference,
  width = 8,
  height = 4.5
)

# creat cutadapt directoties for primers trimming
primerfree_dir <- file.path(common_files, "primer_removed")
cutadapt_log_dir <- file.path(my_root, "logs", "cutadapt")

dir.create(
  primerfree_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  cutadapt_log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

#Build one row per sample with paired R1/R2 input files
paired_files <- file_info %>%
  select(sample, read, file) %>%
  pivot_wider(
    names_from = read,
    values_from = file
  ) %>%
  arrange(sample)

sample.names <- paired_files$sample
raw_reads_fwd_paths <- paired_files$R1
raw_reads_rev_paths <- paired_files$R2

#Define cutadapt output files
primerfree_reads_fwd_paths <- file.path(
  primerfree_dir,
  paste0(sample.names, "_R1_primerfree.fastq.gz")
)

primerfree_reads_rev_paths <- file.path(
  primerfree_dir,
  paste0(sample.names, "_R2_primerfree.fastq.gz")
)

cut_logs <- file.path(
  cutadapt_log_dir,
  paste0(sample.names, ".log")
)

#check cutadapt availability in environment
cutadapt <- Sys.which("cutadapt")

if (!nzchar(cutadapt)) {
  cutadapt <- paste0(
    "C:/Users/riri4/OneDrive/Bureau/Stage_Bernard/",
    "16S_NCCR/cutadapt.exe"
  )
}

if (!file.exists(cutadapt)) {
  stop(
    "Cutadapt executable was not found: ",
    cutadapt
  )
}
cutadapt_version <- system2(
  cutadapt,
  args = "--version",
  stdout = TRUE,
  stderr = TRUE
)

print(cutadapt_version)

#define cutadapt arguments 
cutadapt_args <- c(
  "-g", fwd.primer,
  "-a", rev.primer.rev,
  "-G", rev.primer,
  "-A", fwd.primer.rev,
  "--times", "2",
  "--error-rate", "0.1",
  "--minimum-length", "10",
  "--overlap", "12",
  "--discard-untrimmed",
  "--match-read-wildcards"
)

cutadapt_status <- integer(length(sample.names))

# loop over files
for (i in seq_along(sample.names)) {
  
  message(
    "Running Cutadapt for ",
    sample.names[i],
    " (",
    i,
    "/",
    length(sample.names),
    ")"
  )
  
  cutadapt_status[i] <- system2(
    command = cutadapt,
    args = c(
      cutadapt_args,
      "-o", primerfree_reads_fwd_paths[i],
      "-p", primerfree_reads_rev_paths[i],
      raw_reads_fwd_paths[i],
      raw_reads_rev_paths[i]
    ),
    stdout = cut_logs[i],
    stderr = cut_logs[i]
  )
}

cutadapt_run_summary <- tibble(
  sample = sample.names,
  status = cutadapt_status,
  output_R1_exists = file.exists(primerfree_reads_fwd_paths),
  output_R2_exists = file.exists(primerfree_reads_rev_paths),
  log_exists = file.exists(cut_logs)
)
#a status of 0 means cutadapt completed successfuly 
print(cutadapt_run_summary)

#verify trimming 
primerfree_file_info <- tibble(
  sample = rep(sample.names, each = 2),
  read = rep(c("R1", "R2"), times = length(sample.names)),
  file = as.vector(
    rbind(
      primerfree_reads_fwd_paths,
      primerfree_reads_rev_paths
    )
  )
)
primer_detection_after <- primerfree_file_info %>%
  mutate(
    detection = map2(
      file,
      read,
      ~count_primers_in_fastq(.x, .y)
    )
  ) %>%
  select(sample, read, detection) %>%
  unnest(detection) %>%
  arrange(sample, read, primer)

print(primer_detection_after, n = Inf)

write_tsv(
  primer_detection_after,
  file.path(common_files, "primer_detection_after_long.tsv")
)

# Directory containing Cutadapt output files
primerfree_dir <- file.path(common_files, "primer_removed")

if (!dir.exists(primerfree_dir)) {
  stop("primer_removed directory does not exist: ", primerfree_dir)
}

# Get primer-free file paths
primerfree_reads_fwd_paths <- sort(
  list.files(
    primerfree_dir,
    pattern = "_R1_primerfree\\.fastq\\.gz$",
    full.names = TRUE
  )
)

primerfree_reads_rev_paths <- sort(
  list.files(
    primerfree_dir,
    pattern = "_R2_primerfree\\.fastq\\.gz$",
    full.names = TRUE
  )
)

#pair files
primerfree_file_info <- tibble(
  file = c(
    primerfree_reads_fwd_paths,
    primerfree_reads_rev_paths
  ),
  read = c(
    rep("R1", length(primerfree_reads_fwd_paths)),
    rep("R2", length(primerfree_reads_rev_paths))
  )
) %>%
  mutate(
    filename = basename(file),
    sample = str_remove(
      filename,
      "_R[12]_primerfree\\.fastq\\.gz$"
    )
  )
primerfree_paired <- primerfree_file_info %>%
  select(sample, read, file) %>%
  pivot_wider(
    names_from = read,
    values_from = file
  ) %>%
  arrange(sample)

sample.names <- primerfree_paired$sample

primerfree_reads_fwd_paths <- primerfree_paired$R1
primerfree_reads_rev_paths <- primerfree_paired$R2

names(primerfree_reads_fwd_paths) <- sample.names
names(primerfree_reads_rev_paths) <- sample.names

#trimmed reads quality profiles (plot)
quals_primerfree_R1 <- plotQualityProfile(
  primerfree_reads_fwd_paths,
  aggregate = TRUE
) +
  ggtitle("Primer-free R1 reads") +
  labs(
    x = "Read position",
    y = "Quality score"
  ) +
  scale_x_continuous(
    breaks = seq(0, 300, 25)
  ) +
  theme_bw()

quals_primerfree_R2 <- plotQualityProfile(
  primerfree_reads_rev_paths,
  aggregate = TRUE
) +
  ggtitle("Primer-free R2 reads") +
  labs(
    x = "Read position",
    y = "Quality score"
  ) +
  scale_x_continuous(
    breaks = seq(0, 300, 25)
  ) +
  theme_bw()

ggsave(
  file.path(my_plots, "quals_primerfree_R1.pdf"),
  plot = quals_primerfree_R1,
  device = "pdf",
  width = 8,
  height = 6
)

ggsave(
  file.path(my_plots, "quals_primerfree_R2.pdf"),
  plot = quals_primerfree_R2,
  device = "pdf",
  width = 8,
  height = 6
)

#Inspect primer-free reads before DADA2 quality filtering
#Quality profiles for each forward-read sample
quality_fwd <- plotQualityProfile(
  primerfree_reads_fwd_paths,
  aggregate = FALSE
) +
  ggtitle("Primer-free forward reads (R1)")

quality_fwd

ggsave(
  filename = file.path(
    my_plots,
    "quality_fwd.pdf"
  ),
  plot = quality_fwd,
  width = 8,
  height = 6,
  units = "in"
)

# Quality profiles for each reverse-read sample
quality_rev <- plotQualityProfile(
  primerfree_reads_rev_paths,
  aggregate = FALSE
) +
  ggtitle("Primer-free reverse reads (R2)")

quality_rev

ggsave(
  filename = file.path(
    my_plots,
    "quality_rev.pdf"
  ),
  plot = quality_rev,
  width = 8,
  height = 6,
  units = "in"
)

#or one plot per samples 
for (i in seq_along(sample.names)) {
  
  p_R1 <- plotQualityProfile(
    primerfree_reads_fwd_paths[i]
  ) +
    ggtitle(paste(sample.names[i], "- R1"))
  
  p_R2 <- plotQualityProfile(
    primerfree_reads_rev_paths[i]
  ) +
    ggtitle(paste(sample.names[i], "- R2"))
  
  ggsave(
    file.path(
      my_plots,
      paste0(sample.names[i], "_quality_R1.pdf")
    ),
    p_R1,
    width = 8,
    height = 6
  )
  
  ggsave(
    file.path(
      my_plots,
      paste0(sample.names[i], "_quality_R2.pdf")
    ),
    p_R2,
    width = 8,
    height = 6
  )
}

#trimming and filtering reads
filtered_dir <- file.path(common_files, "filtered_reads")

dir.create(
  filtered_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

filtered_reads_fwd_paths <- file.path(
  filtered_dir,
  paste0(sample.names, "_R1_filtered.fastq.gz")
)

filtered_reads_rev_paths <- file.path(
  filtered_dir,
  paste0(sample.names, "_R2_filtered.fastq.gz")
)

names(filtered_reads_fwd_paths) <- sample.names
names(filtered_reads_rev_paths) <- sample.names

#options tell filterAndTrim() how long the reads should remain, which low-quality reads to discard, 
#truncLen = c(260, 230) the reads should remain 260 and 230 bp after trimming 
#maxN = 0 Discard reads containing any ambiguous N
#maxEE = c(2, 5) Allow at most 2 expected errors in R1 and 5 in R2
#truncQ = 2 Cut a read at the first base with quality ≤2
#rm.phix = TRUE Remove reads matching the PhiX control genome
#compress = TRUE Write gzipped FASTQ output
#multithread = FALSE Use one processing thread
filter_summary <- filterAndTrim(
  fwd = primerfree_reads_fwd_paths,
  filt = filtered_reads_fwd_paths,
  rev = primerfree_reads_rev_paths,
  filt.rev = filtered_reads_rev_paths,
  truncLen = c(260, 230),
  maxN = 0,
  maxEE = c(2, 5),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE,
  verbose = TRUE
)

filter_summary

#check filtering retention
filter_summary_df <- as.data.frame(filter_summary) %>%
  rownames_to_column("sample") %>%
  as_tibble() %>%
  rename(
    reads_in = reads.in,
    reads_out = reads.out
  ) %>%
  mutate(
    retained_percentage = 100 * reads_out / reads_in
  )

print(filter_summary_df)

write_tsv(
  filter_summary_df,
  file.path(
    common_files,
    "filtering_retention.tsv"
  )
)

#check quality profiles after trimming and filtering 
quality_filtered_R1 <- plotQualityProfile(
  filtered_reads_fwd_paths,
  aggregate = TRUE
) +
  ggtitle("Filtered R1 reads") +
  labs(
    x = "Read position",
    y = "Quality score"
  ) +
  theme_bw()

quality_filtered_R2 <- plotQualityProfile(
  filtered_reads_rev_paths,
  aggregate = TRUE
) +
  ggtitle("Filtered R2 reads") +
  labs(
    x = "Read position",
    y = "Quality score"
  ) +
  theme_bw()

ggsave(
  filename = file.path(
    my_plots,
    "quality_filtered_R1.pdf"
  ),
  plot = quality_filtered_R1,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(
    my_plots,
    "quality_filtered_R2.pdf"
  ),
  plot = quality_filtered_R2,
  width = 8,
  height = 6
)

#check number of reads at each processing stage 
# Count FASTQ records and organize them in long format
count_stage_reads <- function(
    paths,
    samples,
    mate,
    stage
) {
  
  if (length(paths) != length(samples)) {
    stop(
      "The numbers of files and sample names differ for stage: ",
      stage,
      ", mate: ",
      mate
    )
  }
  
  if (!all(file.exists(paths))) {
    stop(
      "Some FASTQ files do not exist for stage: ",
      stage,
      ", mate: ",
      mate
    )
  }
  
  counts <- vapply(
    paths,
    function(path) {
      as.numeric(
        ShortRead::countFastq(path)$records
      )
    },
    numeric(1)
  )
  
  tibble(
    sample = samples,
    mate = mate,
    stage = stage,
    n_reads = counts
  )
}

read_counts <- bind_rows(
  
  # Raw reads
  count_stage_reads(
    raw_reads_fwd_paths,
    sample.names,
    mate = "R1",
    stage = "Raw"
  ),
  
  count_stage_reads(
    raw_reads_rev_paths,
    sample.names,
    mate = "R2",
    stage = "Raw"
  ),
  
  # Primer-free reads after Cutadapt
  count_stage_reads(
    primerfree_reads_fwd_paths,
    sample.names,
    mate = "R1",
    stage = "Primer removed"
  ),
  
  count_stage_reads(
    primerfree_reads_rev_paths,
    sample.names,
    mate = "R2",
    stage = "Primer removed"
  ),
  
  # Reads after DADA2 filterAndTrim
  count_stage_reads(
    filtered_reads_fwd_paths,
    sample.names,
    mate = "R1",
    stage = "Filtered"
  ),
  
  count_stage_reads(
    filtered_reads_rev_paths,
    sample.names,
    mate = "R2",
    stage = "Filtered"
  )
) %>%
  mutate(
    stage = factor(
      stage,
      levels = c(
        "Raw",
        "Primer removed",
        "Filtered"
      )
    ),
    mate = factor(
      mate,
      levels = c("R1", "R2")
    )
  ) %>%
  arrange(
    sample,
    stage,
    mate
  )

print(read_counts, n = Inf)

write_tsv(
  read_counts,
  file.path(
    common_files,
    "read_counts_preprocessing.tsv"
  )
)

read_count_pair_check <- read_counts %>%
  pivot_wider(
    names_from = mate,
    values_from = n_reads
  ) %>%
  mutate(
    difference = R1 - R2
  )

print(read_count_pair_check, n = Inf)

if (any(read_count_pair_check$difference != 0)) {
  warning(
    "Some processing stages have different R1 and R2 counts."
  )
}

#plot reads count 
read_counts <- read_tsv(
  file.path(
    common_files,
    "read_counts_preprocessing.tsv"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    stage = factor(
      stage,
      levels = c(
        "Raw",
        "Primer removed",
        "Filtered"
      )
    )
  )

preprocessing_dotplot <- ggplot(
  read_counts,
  aes(
    x = stage,
    y = n_reads,
    color = mate,
    group = interaction(sample, mate)
  )
) +
  geom_line(
    alpha = 0.45,
    linewidth = 0.7,
    position = position_dodge(width = 0.15)
  ) +
  geom_point(
    size = 3,
    position = position_dodge(width = 0.15)
  ) +
  facet_wrap(
    ~sample,
    nrow = 1
  ) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  scale_color_manual(
    values = c(
      "R1" = "#0072B2",
      "R2" = "#D55E00"
    )
  ) +
  labs(
    x = "Preprocessing stage",
    y = "Number of read pairs",
    color = "Mate",
    title = "Read retention during preprocessing",
    subtitle = "Raw reads, primer removal and quality filtering"
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

preprocessing_dotplot

ggsave(
  filename = file.path(
    my_plots,
    "preprocessing_read_counts.pdf"
  ),
  plot = preprocessing_dotplot,
  device = "pdf",
  width = 10,
  height = 4.5
)

#read retention percentage 
read_retention <- read_counts %>%
  group_by(sample, mate) %>%
  arrange(stage, .by_group = TRUE) %>%
  mutate(
    retained_from_previous =
      100 * n_reads / lag(n_reads),
    retained_from_raw =
      100 * n_reads / first(n_reads)
  ) %>%
  ungroup()

print(read_retention, n = Inf)

write_tsv(
  read_retention,
  file.path(
    common_files,
    "read_retention_preprocessing.tsv"
  )
)
