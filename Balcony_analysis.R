###############################################################################
# Inoculant persistence and community legacy in the tomato phyllosphere
# Analysis code for Batarseh et al. (FEMS Microbiology Letters)
#
# Reproduces the statistics, Figures 1-5 and S2-S10, and Supplementary
# Tables S1-S13 reported in the paper.
#
# How to run: put the input files in data/ (see README.md), set the working
# directory to the repository folder, and run the whole script, e.g.
#   source("Balcony_analysis.R")
# Outputs are written to output/ (figures, tables, figure legends, session info).
#
# Inputs (data/):
#   balcony_asv_table.rds                    ASV count table (66 plants + 2 extraction
#                                            blanks) and SILVA v138.1 taxonomy from DADA2
#   2024-10-15_Bacteria_Persistence.xlsx     leaflet plating at T1-T3: CFU, focal-strain
#                                            detection and fungal presence (sheet "Samples_Only")
#   2026-10-28_Balcony_TopLeaf_scoring.xlsx  focal-strain detection on the top leaf
#   2024-10-15_PlantTraits.xlsx              plant traits at harvest
#   isolate_16S_all_copies.fasta             16S rRNA gene copies from SynCom genome assemblies
#
# Section 1 (optional, run_dada2 = TRUE) rebuilds balcony_asv_table.rds from
# primer-trimmed reads (NCBI BioProject PRJNA1556283).
###############################################################################


# =============================================================================
# 0. SETTINGS
# =============================================================================
data_dir     <- "data"
out_dir      <- "output"
asv_rds      <- file.path(data_dir, "balcony_asv_table.rds")
persist_xlsx <- file.path(data_dir, "2024-10-15_Bacteria_Persistence.xlsx")
topleaf_xlsx <- file.path(data_dir, "2026-10-28_Balcony_TopLeaf_scoring.xlsx")
traits_xlsx  <- file.path(data_dir, "2024-10-15_PlantTraits.xlsx")
iso_fasta    <- file.path(data_dir, "isolate_16S_all_copies.fasta")

run_dada2    <- FALSE                       # TRUE to rebuild the ASV table from reads
fastq_dir    <- "trimmed"                   # primer-trimmed FASTQ files (cutadapt)
silva_train  <- file.path(data_dir, "silva_nr99_v138.1_train_set.fa.gz")

topleaf_date <- "2024-09-08"
# Plants excluded for low read depth after removal of insect-associated genera
low_depth_plants <- c("R_SC01_02", "R_SC02_01", "R_SC03_01", "R_SC04_01",
                      "R_SC09_05", "R_SC25_04", "R_SC35_05")
nperm   <- 999
seed    <- 123
fig_dpi <- 400


# =============================================================================
# 1. (OPTIONAL) DADA2: PRIMER-TRIMMED READS -> ASV TABLE
# =============================================================================
# Primers (799F 5'-AACMGGATTAGATACCCKG-3', 1193R 5'-ACGTCATCCCCACCTTCC-3') were
# removed from the 5' end of each read with cutadapt v4.9 before this step
# (equivalent command in README.md). Reads from the same MiSeq run
# that belong to another study are not part of this dataset, so a rerun learns
# error rates from these samples only and may give very slightly different ASVs.
if (run_dada2) {
  fnFs <- sort(list.files(fastq_dir, pattern = "_R1(_001)?\\.fastq\\.gz$", full.names = TRUE))
  fnRs <- sort(list.files(fastq_dir, pattern = "_R2(_001)?\\.fastq\\.gz$", full.names = TRUE))
  smp  <- sub("_S[0-9]+_R1(_001)?\\.fastq\\.gz$", "", sub("^trimmed_", "", basename(fnFs)))
  filtFs <- setNames(file.path(fastq_dir, "filt", basename(fnFs)), smp)
  filtRs <- setNames(file.path(fastq_dir, "filt", basename(fnRs)), smp)

  filt_out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(270, 260),
                                   maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                                   compress = TRUE, multithread = TRUE)
  ok <- file.exists(filtFs)
  errF <- dada2::learnErrors(filtFs[ok], multithread = TRUE)
  errR <- dada2::learnErrors(filtRs[ok], multithread = TRUE)
  dadaFs  <- dada2::dada(filtFs[ok], err = errF, multithread = TRUE)
  dadaRs  <- dada2::dada(filtRs[ok], err = errR, multithread = TRUE)
  mergers <- dada2::mergePairs(dadaFs, filtFs[ok], dadaRs, filtRs[ok])
  seqtab  <- dada2::makeSequenceTable(mergers)
  seqtab_nochim <- dada2::removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
  taxa <- dada2::assignTaxonomy(seqtab_nochim, silva_train, multithread = TRUE, tryRC = TRUE)

  getN  <- function(x) sum(dada2::getUniques(x))
  track <- data.frame(sample_id = smp[ok], input = filt_out[ok, 1], filtered = filt_out[ok, 2],
                      denoisedF = sapply(dadaFs, getN), denoisedR = sapply(dadaRs, getN),
                      merged = sapply(mergers, getN), nonchim = rowSums(seqtab_nochim))
  dir.create(out_dir, showWarnings = FALSE)
  write.csv(track, file.path(out_dir, "dada2_read_tracking.csv"), row.names = FALSE)
  saveRDS(list(counts = seqtab_nochim, taxonomy = taxa), asv_rds)
}


# =============================================================================
# 2. PACKAGES, OUTPUT FOLDERS AND PLOTTING HELPERS
# =============================================================================
needed  <- c("phyloseq", "vegan", "permute", "readxl", "dplyr", "tidyr",
             "ggplot2", "patchwork", "pairwiseAdonis", "openxlsx")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "),
                          "\n(pairwiseAdonis: devtools::install_github('pmartinezarbizu/pairwiseAdonis/pairwiseAdonis'))")

suppressPackageStartupMessages({
  library(phyloseq); library(vegan); library(permute); library(readxl)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(pairwiseAdonis)
})
options(dplyr.summarise.inform = FALSE)

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = fig_dpi, bg = "white")
}
tabs <- list()                                   # results collected for the supplementary workbook
keep_tab <- function(x, name) tabs[[name]] <<- as.data.frame(x)

# Colorblind-safe palettes: Okabe-Ito (groups) and Paul Tol muted (inocula, genera)
okabe <- c(grey = "#999999", orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
           yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00", pink = "#CC79A7")
tol_muted <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
               "#DDCC77", "#CC6677", "#882255", "#AA4499")
group_cols <- c("Control"               = okabe[["grey"]],
                "Not persistent"        = okabe[["vermillion"]],
                "Persistent"            = okabe[["blue"]],
                "Not detected T3"       = okabe[["vermillion"]],
                "Detected T3"           = okabe[["blue"]],
                "Non-persister species" = okabe[["vermillion"]],
                "Persister species"     = okabe[["blue"]])

inoc_levels <- c("MgCl2", "SC01", "SC02", "SC03", "SC04", "SC09", "SC13", "SC25", "SC35", "SCB8")
inoc_cols   <- setNames(c("#BBBBBB", tol_muted), inoc_levels)
inoc_species <- c(MgCl2 = "MgCl2 buffer control", SC01 = "Frigoribacterium endophyticum",
                  SC02 = "Pantoea agglomerans", SC03 = "Curtobacterium herbarum",
                  SC04 = "Erwinia tasmaniensis", SC09 = "Massilia aurea",
                  SC13 = "Exiguobacterium sibiricum", SC25 = "Pseudomonas rhodesiae",
                  SC35 = "Bacillus thuringiensis", SCB8 = "SynCom of all 8 strains")
inoc_key <- function(drop_control = FALSE) {
  k <- if (drop_control) inoc_species[-1] else inoc_species
  paste(sprintf("%s = %s", names(k), k), collapse = "; ")
}
inoc_from_id <- function(id) ifelse(grepl("^MgCl2", id, ignore.case = TRUE), "MgCl2",
                                    sub("^R_(SC[0-9B]+)_.*", "\\1", id))

# Axis and facet labels: italic species name + strain code, e.g. E. sibiricum (SC13)
inoc_abbrev <- c(SC01 = "F. endophyticum", SC02 = "P. agglomerans", SC03 = "C. herbarum",
                 SC04 = "E. tasmaniensis", SC09 = "M. aurea", SC13 = "E. sibiricum",
                 SC25 = "P. rhodesiae", SC35 = "B. thuringiensis")
inoc_pm <- setNames(vapply(inoc_levels, function(k) {
  if (k == "MgCl2") return("MgCl[2]~'control'")
  if (k == "SCB8")  return("plain('SynCom')")
  sprintf("italic('%s')~'(%s)'", inoc_abbrev[[k]], k)
}, character(1)), inoc_levels)
inoc_axis_labels <- function(x) parse(text = inoc_pm[as.character(x)])
inoc_lab     <- function(x) factor(inoc_pm[as.character(x)], levels = inoc_pm)
scale_x_inoc <- function(...) scale_x_discrete(labels = inoc_axis_labels, ...)

theme_set(theme_classic())
th_box   <- theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
th_pcoa  <- theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
                  legend.position = "bottom", legend.text = element_text(size = 12))
th_bar   <- theme(axis.text.x = element_text(size = 7), legend.position = "right",
                  legend.text = element_text(size = 8), strip.text = element_text(face = "bold", size = 9))
th_strip <- theme(strip.text = element_text(face = "bold", size = 9))

# Figures have no titles; titles and legends (with this run's statistics) go to figure_legends.md
fig_legends <- list()
add_legend <- function(fig, title, text) fig_legends[[fig]] <<- c(paste0("**", title, "**"), "", text)
set.seed(seed)


# =============================================================================
# 3. 16S DATA: BUILD PHYLOSEQ OBJECT AND FILTER READS
# =============================================================================
asv    <- readRDS(asv_rds)
counts <- asv$counts[grepl("^R_SC|^MgCl2|^Blank", rownames(asv$counts), ignore.case = TRUE), , drop = FALSE]
smeta  <- data.frame(sample_id = rownames(counts), row.names = rownames(counts)) %>%
  mutate(sample_type = ifelse(grepl("^Blank", sample_id), "extraction_blank", "plant"),
         inoculant   = ifelse(sample_type == "plant", inoc_from_id(sample_id), NA_character_),
         replicate   = as.numeric(sub(".*_(\\d+)$", "\\1", sample_id)))
ps_all <- phyloseq(otu_table(counts, taxa_are_rows = FALSE), sample_data(smeta),
                   tax_table(asv$taxonomy))

blanks <- prune_samples(smeta$sample_type == "extraction_blank", ps_all)
plants <- setdiff(smeta$sample_id[smeta$sample_type == "plant"], low_depth_plants)
ps_raw <- prune_samples(plants, ps_all)
stopifnot(nsamples(ps_raw) == 59, nsamples(blanks) == 2)

# Remove ASVs found in extraction blanks (contaminants), sap-feeding insect
# endosymbionts and bee-gut core genera (not leaf colonists)
insect_genera <- c("Buchnera", "Arsenophonus", "Candidatus Regiella", "Candidatus Hamiltonella",
                   "Candidatus Tremblaya", "Candidatus Portiera", "Wolbachia", "Rickettsia", "Sodalis")
bee_genera    <- c("Bombilactobacillus", "Gilliamella", "Snodgrassella", "Bifidobacterium",
                   "Bartonella", "Apilactobacillus")
genus_raw <- as(tax_table(ps_raw), "matrix")[, "Genus"]
is_insect <- !is.na(genus_raw) & genus_raw %in% insect_genera
is_bee    <- !is.na(genus_raw) & genus_raw %in% bee_genera
is_blank  <- taxa_sums(blanks) > 0

drop_zero <- function(x) prune_taxa(taxa_sums(x) > 0, x)
ps_full  <- drop_zero(prune_taxa(!is_blank, ps_raw))                         # blank ASVs removed only
ps_u     <- drop_zero(prune_taxa(!(is_blank | is_insect | is_bee), ps_raw))  # main dataset, unrarefied
ps_u_bee <- drop_zero(prune_taxa(!(is_blank | is_insect), ps_raw))           # bee-gut genera retained

rare_depth <- min(sample_sums(ps_u))
ps     <- rarefy_even_depth(ps_u, sample.size = rare_depth, rngseed = seed,
                            replace = FALSE, trimOTUs = TRUE, verbose = FALSE)   # main dataset, rarefied
ps_bee <- rarefy_even_depth(ps_u_bee, sample.size = min(sample_sums(ps_u_bee)), rngseed = seed,
                            replace = FALSE, trimOTUs = TRUE, verbose = FALSE)

raw_cnt <- as(otu_table(ps_raw), "matrix")
filter_summary <- data.frame(
  sample_id       = rownames(raw_cnt),
  raw_reads       = rowSums(raw_cnt),
  insect_reads    = rowSums(raw_cnt[, is_insect, drop = FALSE]),
  bee_gut_reads   = rowSums(raw_cnt[, is_bee, drop = FALSE]),
  blank_ASV_reads = rowSums(raw_cnt[, is_blank, drop = FALSE]),
  retained_reads  = sample_sums(ps_u)[rownames(raw_cnt)]) %>%
  mutate(pct_insect = 100 * insect_reads / raw_reads, pct_bee_gut = 100 * bee_gut_reads / raw_reads,
         pct_blank_ASVs = 100 * blank_ASV_reads / raw_reads)
keep_tab(filter_summary, "read_filtering")
filter_notes <- c(
  sprintf("ASVs detected in %d extraction blanks: %d (removed); median %.2f%% (max %.2f%%) of reads per plant.",
          nsamples(blanks), sum(is_blank), median(filter_summary$pct_blank_ASVs), max(filter_summary$pct_blank_ASVs)),
  sprintf("Insect endosymbiont genera removed (%s): median %.1f%% of reads per plant (max %.1f%%).",
          paste(insect_genera, collapse = ", "), median(filter_summary$pct_insect), max(filter_summary$pct_insect)),
  sprintf("Bee-gut genera removed (%s): median %.1f%% of reads per plant (max %.1f%%).",
          paste(bee_genera, collapse = ", "), median(filter_summary$pct_bee_gut), max(filter_summary$pct_bee_gut)),
  sprintf("Remaining reads rarefied to %d per plant (n = %d plants).", rare_depth, nsamples(ps)))

# ---- Fig S10: rarefaction curves (filtered, unrarefied reads) -----------------
rc_mat <- as(otu_table(ps_u), "matrix")
if (taxa_are_rows(ps_u)) rc_mat <- t(rc_mat)
rc_df <- rarecurve(rc_mat, step = 100, tidy = TRUE) %>%
  rename(sample_id = Site, reads = Sample, ASVs = Species) %>%
  mutate(inoculant = factor(inoc_from_id(sample_id), levels = inoc_levels))
slope_100 <- 100 * rareslope(rc_mat, sample = rare_depth)   # new ASVs per 100 extra reads at the rarefaction depth
figS10 <- ggplot(rc_df, aes(reads, ASVs, group = sample_id, colour = inoculant)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_vline(xintercept = rare_depth, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = inoc_cols, labels = inoc_axis_labels, name = "Inoculum") +
  xlab("Reads sampled") + ylab("ASVs observed")
save_fig(figS10, "FigS10_rarefaction_curves", w = 8, h = 5)
add_legend("FigS10", "Rarefaction curves",
  paste0("ASV richness against the number of reads sampled for each plant (n = ", nsamples(ps_u), " plants), ",
         "after removal of blank-associated ASVs, insect endosymbiont and bee-gut genera, and before rarefaction. ",
         "The dashed line marks the rarefaction depth (", rare_depth, " reads). At this depth plants gained a median of ",
         sprintf("%.2f", median(slope_100)), " new ASVs per 100 additional reads (range ",
         sprintf("%.2f", min(slope_100)), "-", sprintf("%.2f", max(slope_100)), "). Seven plants with low read counts ",
         "were excluded before this step. Colours show the inoculum: ", inoc_key(), "."))


# =============================================================================
# 4. PLATING DATA AND PERSISTENCE GROUPS
# =============================================================================
pers_raw <- read_excel(persist_xlsx, sheet = "Samples_Only", .name_repair = "minimal")
names(pers_raw)[15] <- "CFU_g_KBHA"     # two columns share the name CFU_gram_leaf_tissue
names(pers_raw)[21] <- "CFU_g_rif"

pers_plant <- pers_raw %>%
  transmute(Isolate, Block, Sampling_Timepoint,
            CFU_g_KBHA = as.numeric(CFU_g_KBHA), CFU_g_rif = as.numeric(CFU_g_rif),
            detected = as.numeric(Focal_Bacteria_Persistence)) %>%
  group_by(Isolate) %>%
  summarise(Block        = dplyr::first(Block),
            T1           = detected[Sampling_Timepoint == "T1"],
            T2           = detected[Sampling_Timepoint == "T2"],
            T3           = detected[Sampling_Timepoint == "T3"],
            n_detect     = T1 + T2 + T3,
            logrif_mean  = mean(log10(CFU_g_rif + 1)),
            logKBHA_mean = mean(log10(CFU_g_KBHA + 1)))

# Top (youngest) leaf, which developed after inoculation; not used in the persistence groups
topleaf <- read_excel(topleaf_xlsx, sheet = "TopLeaf_scoring") %>%
  transmute(Isolate, TopLeaf = suppressWarnings(as.numeric(TopLeaf_focal_detected)))
pers_plant <- left_join(pers_plant, topleaf, by = "Isolate")
stopifnot(all(sample_names(ps) %in% pers_plant$Isolate))

# Persistence: plant level (detected at >=2 of 3 timepoints), detected at T3, and species level
persister_species <- c("SC01", "SC02", "SC03", "SC04", "SC25", "SCB8")
meta <- data.frame(sample_data(ps)) %>%
  select(sample_id, inoculant, replicate) %>%
  left_join(pers_plant, by = c("sample_id" = "Isolate")) %>%
  mutate(
    inoculant = factor(inoculant, levels = inoc_levels),
    persist_plant = factor(case_when(inoculant == "MgCl2" ~ "Control",
                                     n_detect >= 2        ~ "Persistent",
                                     TRUE                 ~ "Not persistent"),
                           levels = c("Control", "Not persistent", "Persistent")),
    persist_T3 = factor(case_when(inoculant == "MgCl2" ~ "Control",
                                  T3 == 1              ~ "Detected T3",
                                  TRUE                 ~ "Not detected T3"),
                        levels = c("Control", "Not detected T3", "Detected T3")),
    persist_species = factor(case_when(inoculant == "MgCl2"             ~ "Control",
                                       inoculant %in% persister_species ~ "Persister species",
                                       TRUE                             ~ "Non-persister species"),
                             levels = c("Control", "Non-persister species", "Persister species")),
    Block = factor(Block))
meta$species <- unname(inoc_species[as.character(meta$inoculant)])
rownames(meta) <- meta$sample_id

sample_data(ps)      <- sample_data(meta)
sample_data(ps_u)    <- sample_data(meta[sample_names(ps_u), ])
sample_data(ps_full) <- sample_data(meta[sample_names(ps_full), ])
sample_data(ps_bee)  <- sample_data(meta[sample_names(ps_bee), ])

schemes <- c(persist_plant   = "Plant level (detected at >=2 of 3 timepoints)",
             persist_T3      = "Detected at T3",
             persist_species = "Species level")
keep_tab(as.data.frame.matrix(table(meta$inoculant, meta$persist_plant)) %>%
           tibble::rownames_to_column("inoculant"), "group_counts")

# ---- Top leaf vs T3 lower leaflet (all inoculated plants) --------------------
tl_all <- pers_plant %>%
  mutate(inoculant = factor(inoc_from_id(Isolate), levels = inoc_levels)) %>%
  filter(inoculant != "MgCl2") %>% droplevels()
tl_by_inoc <- tl_all %>%
  group_by(inoculant) %>%
  summarise(n_plants = n(), T3_detected = sum(T3), TopLeaf_scored = sum(!is.na(TopLeaf)),
            TopLeaf_detected = sum(TopLeaf, na.rm = TRUE),
            T3_pos_and_TopLeaf_pos = sum(T3 == 1 & TopLeaf %in% 1),
            T3_pos_and_TopLeaf_neg = sum(T3 == 1 & TopLeaf %in% 0)) %>%
  ungroup()
tl_by_inoc <- bind_rows(tl_by_inoc %>% mutate(inoculant = as.character(inoculant)),
                        tl_by_inoc %>% summarise(inoculant = "All inoculated",
                                                 across(where(is.numeric), sum)))
keep_tab(tl_by_inoc, "topleaf_by_inoculant")

tl_pair  <- tl_all %>% filter(!is.na(TopLeaf))
n_lost   <- sum(tl_pair$T3 == 1 & tl_pair$TopLeaf == 0)
n_gained <- sum(tl_pair$T3 == 0 & tl_pair$TopLeaf == 1)
mcn_p    <- binom.test(n_lost, n_lost + n_gained)$p.value      # exact McNemar test
tl_test  <- data.frame(
  test = "Focal detection, T3 lower leaflet vs top leaf (exact McNemar, inoculated plants scored at both)",
  n_plants = nrow(tl_pair), T3_detected = sum(tl_pair$T3), TopLeaf_detected = sum(tl_pair$TopLeaf),
  detected_T3_only = n_lost, detected_TopLeaf_only = n_gained, p = mcn_p,
  n_topleaf_unscored = sum(is.na(tl_all$TopLeaf)), topleaf_date = topleaf_date)
keep_tab(tl_test, "topleaf_test")

# ---- Fig 1: focal-strain detection heatmap ----------------------------------
tp_labels <- c(T1 = "T1", T2 = "T2", T3 = "T3", TopLeaf = "Top leaf")
fig1 <- meta %>%
  filter(inoculant != "MgCl2") %>%
  select(sample_id, inoculant, T1, T2, T3, TopLeaf) %>%
  pivot_longer(c(T1, T2, T3, TopLeaf), names_to = "timepoint", values_to = "detected") %>%
  mutate(inoculant = droplevels(inoc_lab(inoculant)),
         timepoint = factor(tp_labels[timepoint], levels = tp_labels),
         detected  = factor(ifelse(is.na(detected), "NA", as.character(detected)), levels = c("0", "1", "NA")),
         plant     = sub("^R_SC[0-9B]+_", "", sample_id)) %>%
  ggplot(aes(timepoint, plant, fill = detected)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = function(d) d[d$detected == "NA", ], fill = "white", colour = "grey55",
            linewidth = 0.4, width = 0.9, height = 0.9, show.legend = FALSE) +
  geom_text(data = function(d) d[d$detected == "NA", ], label = "–", colour = "grey45",
            size = 4, show.legend = FALSE) +
  geom_vline(xintercept = 3.5, linetype = 2, colour = "grey40") +
  facet_wrap(~inoculant, nrow = 2, scales = "free_y", labeller = label_parsed) +
  scale_fill_manual(values = c("0" = "#EEEEEE", "1" = okabe[["blue"]], "NA" = "#FFFFFF"),
                    labels = c("0" = "Not detected", "1" = "Detected", "NA" = "Not scored"),
                    drop = FALSE, name = "Focal strain\non KBHA-rif") +
  guides(fill = guide_legend(override.aes = list(colour = c(NA, NA, "grey55")))) +
  xlab("Lower leaflet (T1-T3) | new growth") + ylab("Replicate plant") +
  th_strip + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(fig1, "Fig1_focal_strain_detection_heatmap", w = 11, h = 5.5)
tl_tot <- tl_by_inoc %>% filter(inoculant == "All inoculated")
add_legend("Fig1", "Focal strains persisted on inoculated leaflets but rarely reached new leaves",
  paste0("Detection of the rifampicin-marked focal strain on KBHA-rif plates from one lower leaflet per plant at ",
         "T1 (2024-08-10), T2 (2024-08-25) and T3 (2024-09-08), and from the top (youngest) leaf (", topleaf_date,
         "), which grew after the seedlings were dip-inoculated. Detection is at least one colony with the expected ",
         "morphology. Only plants retained in the 16S dataset are shown (n = ", sum(meta$inoculant != "MgCl2"),
         " inoculated plants). Outlined cells with a dash are top leaves that could not be scored. Of the ",
         tl_test$n_plants, " inoculated plants scored on both leaves (all ", tl_tot$n_plants, " inoculated plants, ",
         tl_test$n_topleaf_unscored, " top leaves unscored), the focal strain was detected on ", tl_test$T3_detected,
         " T3 leaflets and ", tl_test$TopLeaf_detected, " top leaves (", tl_test$detected_T3_only, " T3 only, ",
         tl_test$detected_TopLeaf_only, " top leaf only; exact McNemar p = ", signif(mcn_p, 2), "). ",
         "16S samples were collected from whole-shoot sonicates at harvest on 10 September 2024."))


# =============================================================================
# 5. COMMUNITY COMPOSITION BY PERSISTENCE GROUP
# =============================================================================
bray   <- phyloseq::distance(ps, method = "bray")
bray_m <- as.matrix(bray)
ord    <- ordinate(ps, method = "PCoA", distance = bray)

adonis_row <- function(a) data.frame(F = a$F[1], R2 = a$R2[1], p = a$`Pr(>F)`[1])
beta_tests <- function(dist, data, grp) {
  set.seed(seed)
  f   <- as.formula(paste("dist ~", grp))
  ov  <- adonis2(f, data = data, permutations = nperm)
  pw  <- pairwise.adonis2(f, data = data, permutations = nperm)
  pw  <- pw[names(pw) != "parent_call"]
  dsp <- permutest(betadisper(dist, data[[grp]], type = "centroid"), permutations = nperm)
  bind_rows(
    data.frame(scheme = grp, test = "PERMANOVA overall", comparison = "all groups", adonis_row(ov)),
    bind_rows(lapply(names(pw), function(n)
      data.frame(scheme = grp, test = "PERMANOVA pairwise", comparison = gsub("_vs_", " vs ", n),
                 adonis_row(pw[[n]])))),
    data.frame(scheme = grp, test = "Dispersion (betadisper)", comparison = "all groups",
               F = dsp$tab$F[1], R2 = NA, p = dsp$tab$`Pr(>F)`[1]))
}
permanova_tab <- bind_rows(lapply(names(schemes), function(g) beta_tests(bray, meta, g)))
keep_tab(permanova_tab, "permanova_by_scheme")

# ---- Fig 2 and Fig S5: PCoA ---------------------------------------------------
eig     <- round(100 * ord$values$Relative_eig[1:2], 1)
pcoa_df <- data.frame(ord$vectors[, 1:2], meta[rownames(ord$vectors), ])
p_lab   <- permanova_tab %>% filter(test == "PERMANOVA overall")
plot_pcoa <- function(df, grp, eig) {
  ggplot(df, aes(Axis.1, Axis.2, colour = .data[[grp]])) +
    geom_point(size = 4) +
    stat_ellipse(level = 0.95) +
    scale_colour_manual(values = group_cols) +
    labs(x = paste0("PCoA1 (", eig[1], "%)"), y = paste0("PCoA2 (", eig[2], "%)"), colour = NULL) +
    th_pcoa
}
fig2 <- plot_pcoa(pcoa_df, "persist_plant", eig) + plot_pcoa(pcoa_df, "persist_species", eig) +
  plot_annotation(tag_levels = "A")
save_fig(fig2, "Fig2_PCoA_by_persistence", w = 12, h = 6)
save_fig(plot_pcoa(pcoa_df, "persist_T3", eig), "FigS5_PCoA_by_T3_detection", w = 6.5, h = 6)
pl <- function(g) { x <- p_lab %>% filter(scheme == g); sprintf("R² = %.3f, p = %.3f", x$R2, x$p) }
add_legend("Fig2", "Leaf bacterial community composition does not differ by focal-strain persistence",
  paste0("Principal coordinates analysis of Bray-Curtis dissimilarities (filtered, rarefied to ", rare_depth,
         " reads; n = ", nrow(meta), " plants). Plants are grouped by (A) plant-level persistence, focal strain ",
         "detected at >=2 of 3 timepoints (PERMANOVA ", pl("persist_plant"), "); (B) species-level persistence, ",
         "with SC09, SC13 and SC35 as non-persister species (", pl("persist_species"), "). Ellipses are 95% ",
         "confidence ellipses. ", nperm, " permutations."))
add_legend("FigS5", "Community composition grouped by focal-strain detection at T3",
  paste0("Same ordination as Fig 2, grouped by detection of the focal strain at T3 (PERMANOVA ",
         pl("persist_T3"), ")."))


# =============================================================================
# 6. ARE PLANTS THAT LOST THEIR INOCULANT CLOSER TO THE BUFFER CONTROLS?
# =============================================================================
ctl_ids <- meta$sample_id[meta$inoculant == "MgCl2"]
meta$dist_to_ctl <- vapply(meta$sample_id, function(s) mean(bray_m[s, setdiff(ctl_ids, s)]), numeric(1))
inoc_only <- droplevels(filter(meta, inoculant != "MgCl2"))

within_bc <- function(grp) {
  bind_rows(lapply(levels(meta[[grp]]), function(g) {
    ids <- meta$sample_id[meta[[grp]] == g]
    m <- bray_m[ids, ids]
    data.frame(scheme = grp, group = g, n = length(ids),
               mean_within_group_BC = mean(m[upper.tri(m)]),
               mean_BC_to_controls  = mean(meta$dist_to_ctl[meta[[grp]] == g]))
  }))
}
dist_ctl_summary <- bind_rows(lapply(names(schemes), within_bc))
keep_tab(dist_ctl_summary, "dissimilarity_summary")

w1 <- wilcox.test(dist_to_ctl ~ persist_plant,   data = inoc_only, exact = FALSE)
w2 <- wilcox.test(dist_to_ctl ~ persist_species, data = inoc_only, exact = FALSE)
c1 <- suppressWarnings(cor.test(inoc_only$n_detect,    inoc_only$dist_to_ctl, method = "spearman"))
c2 <- suppressWarnings(cor.test(inoc_only$logrif_mean, inoc_only$dist_to_ctl, method = "spearman"))

ps_inoc   <- subset_samples(ps, inoculant != "MgCl2")
meta_inoc <- data.frame(sample_data(ps_inoc))
bray_inoc <- phyloseq::distance(ps_inoc, method = "bray")
set.seed(seed)
a_within  <- adonis2(bray_inoc ~ persist_plant, data = meta_inoc,
                     permutations = how(blocks = meta_inoc$inoculant, nperm = nperm))
set.seed(seed)
a_inoc    <- adonis2(bray_inoc ~ inoculant, data = meta_inoc, permutations = nperm)

dist_ctl_tests <- data.frame(
  test = c("Distance to controls: Persistent vs Not persistent (Wilcoxon)",
           "Distance to controls: Persister vs Non-persister species (Wilcoxon)",
           "Distance to controls ~ number of detections (Spearman, inoculated)",
           "Distance to controls ~ mean log10 rif CFU (Spearman, inoculated)",
           "PERMANOVA persistence, permuted within inoculant (inoculated only)",
           "PERMANOVA inoculant identity (inoculated only)"),
  statistic = c(w1$statistic, w2$statistic, c1$estimate, c2$estimate, a_within$F[1], a_inoc$F[1]),
  R2 = c(NA, NA, NA, NA, a_within$R2[1], a_inoc$R2[1]),
  p  = c(w1$p.value, w2$p.value, c1$p.value, c2$p.value, a_within$`Pr(>F)`[1], a_inoc$`Pr(>F)`[1]))
keep_tab(dist_ctl_tests, "distance_tests")

# ---- Fig 3 --------------------------------------------------------------------
fig3 <- ggplot(meta, aes(x = persist_plant, y = dist_to_ctl, fill = persist_plant)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, size = 2) +
  scale_fill_manual(values = group_cols) +
  xlab("Focal strain persistence") + ylab("Mean Bray-Curtis dissimilarity to MgCl2 plants") +
  th_box
save_fig(fig3, "Fig3_distance_to_controls", w = 6, h = 5)
dsum <- dist_ctl_summary %>% filter(scheme == "persist_plant")
add_legend("Fig3", "Plants that lost the focal strain are not more similar to buffer controls",
  paste0("Mean Bray-Curtis dissimilarity of each plant to the MgCl2 control plants (controls compared ",
         "with the other controls). Mean distance: Control ", round(dsum$mean_BC_to_controls[1], 3),
         ", Not persistent ", round(dsum$mean_BC_to_controls[2], 3), ", Persistent ",
         round(dsum$mean_BC_to_controls[3], 3), ". Persistent vs Not persistent, Wilcoxon p = ",
         signif(w1$p.value, 2), ". Mean within-group dissimilarity: ",
         paste(sprintf("%s %.3f", dsum$group, dsum$mean_within_group_BC), collapse = ", "), "."))


# =============================================================================
# 7. EACH INOCULUM VS CONTROLS; ALPHA DIVERSITY AND DOMINANCE
# =============================================================================
alpha <- suppressWarnings(estimate_richness(ps, measures = c("Observed", "Shannon"))) %>%
  mutate(sample_id = rownames(.)) %>%
  left_join(meta, by = "sample_id")
rel <- transform_sample_counts(ps, function(x) x / sum(x))
alpha$max_ASV_share <- apply(as(otu_table(rel), "matrix"), 1, max)[alpha$sample_id]

each_vs_ctl <- bind_rows(lapply(setdiff(inoc_levels, "MgCl2"), function(sc) {
  ids <- meta$sample_id[meta$inoculant %in% c(sc, "MgCl2")]
  ps_sub <- prune_samples(ids, ps)
  m <- data.frame(sample_data(ps_sub)) %>% mutate(inoculant = droplevels(inoculant))
  d <- phyloseq::distance(ps_sub, "bray")
  set.seed(seed)
  a  <- adonis2(d ~ inoculant, data = m, permutations = nperm)
  bd <- permutest(betadisper(d, m$inoculant, type = "centroid"), permutations = nperm)
  sh  <- alpha$Shannon[match(ids, alpha$sample_id)]
  grp <- m$inoculant[match(ids, m$sample_id)]
  data.frame(inoculant = sc, n = sum(grp == sc),
             PERMANOVA_R2 = a$R2[1], PERMANOVA_p = a$`Pr(>F)`[1],
             dispersion_p = bd$tab$`Pr(>F)`[1],
             median_Shannon = median(sh[grp == sc]),
             median_Shannon_control = median(sh[grp == "MgCl2"]),
             Shannon_p = wilcox.test(sh[grp == sc], sh[grp == "MgCl2"], exact = FALSE)$p.value)
})) %>% mutate(PERMANOVA_p_BH = p.adjust(PERMANOVA_p, "BH"),
               Shannon_p_BH   = p.adjust(Shannon_p, "BH"))
keep_tab(each_vs_ctl, "each_vs_control")

ids  <- meta$sample_id[meta$inoculant %in% c("SC13", "SC35", "MgCl2")]
mp   <- meta[ids, ] %>% mutate(grp = ifelse(inoculant == "MgCl2", "Control", "SC13+SC35"))
dp   <- as.dist(bray_m[ids, ids])
set.seed(seed)
a_pool  <- adonis2(dp ~ grp, data = mp, permutations = nperm)
bd_pool <- betadisper(dp, mp$grp, type = "centroid")
keep_tab(data.frame(comparison = "SC13 + SC35 vs MgCl2",
                    PERMANOVA_R2 = a_pool$R2[1], PERMANOVA_p = a_pool$`Pr(>F)`[1],
                    dispersion_control   = mean(bd_pool$distances[mp$grp == "Control"]),
                    dispersion_SC13_SC35 = mean(bd_pool$distances[mp$grp == "SC13+SC35"]),
                    dispersion_p = permutest(bd_pool, permutations = nperm)$tab$`Pr(>F)`[1]),
         "pooled_vs_control")

alpha_tests <- bind_rows(lapply(names(schemes), function(g) {
  kw <- kruskal.test(alpha$Shannon ~ alpha[[g]])
  pw <- pairwise.wilcox.test(alpha$Shannon, alpha[[g]], p.adjust.method = "BH", exact = FALSE)$p.value
  pw_long <- as.data.frame(as.table(pw)) %>% filter(!is.na(Freq)) %>%
    transmute(scheme = g, test = "Shannon pairwise Wilcoxon (BH)",
              comparison = paste(Var1, "vs", Var2), p = Freq)
  bind_rows(data.frame(scheme = g, test = "Shannon Kruskal-Wallis", comparison = "all groups",
                       p = kw$p.value), pw_long)
}))
alpha_tests <- bind_rows(alpha_tests,
  data.frame(scheme = "inoculant", test = "Shannon Kruskal-Wallis", comparison = "all inoculants",
             p = kruskal.test(Shannon ~ inoculant, data = alpha)$p.value),
  data.frame(scheme = "persist_plant", test = "Observed ASVs Kruskal-Wallis", comparison = "all groups",
             p = kruskal.test(Observed ~ persist_plant, data = alpha)$p.value))
keep_tab(alpha_tests, "alpha_tests")

# ---- Fig 4 --------------------------------------------------------------------
fig4a <- ggplot(alpha, aes(x = inoculant, y = Shannon, fill = inoculant)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, size = 2) +
  scale_fill_manual(values = inoc_cols) +
  scale_x_inoc() + xlab("Inoculum") + ylab("Shannon Index") + th_box
fig4b <- ggplot(alpha, aes(x = inoculant, y = 100 * max_ASV_share, fill = inoculant)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, size = 2) +
  scale_fill_manual(values = inoc_cols) +
  scale_x_inoc() + xlab("Inoculum") + ylab("Most abundant ASV (% of reads)") + th_box
save_fig(fig4a / fig4b + plot_annotation(tag_levels = "A"),
         "Fig4_Shannon_and_dominance_by_inoculant", w = 7.5, h = 8.5)
add_legend("Fig4", "SC13 and SC35 plants have more diverse, less dominated communities",
  paste0("(A) Shannon diversity of ASVs and (B) share of reads from the single most abundant ASV per plant ",
         "(filtered and rarefied). Non-persister species: SC09, SC13, SC35. Versus MgCl2 controls (Wilcoxon): ",
         paste(sprintf("%s p = %.3f", each_vs_ctl$inoculant, each_vs_ctl$Shannon_p), collapse = "; "),
         ". Inoculants: ", inoc_key(), "."))


# =============================================================================
# 8. GENUS-LEVEL COMPOSITION
# =============================================================================
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
rel_g    <- transform_sample_counts(ps_genus, function(x) x / sum(x))
tt_g     <- as.data.frame(as(tax_table(rel_g), "matrix"))
gen_mat  <- as.data.frame(as(otu_table(rel_g), "matrix"))
colnames(gen_mat) <- make.unique(ifelse(is.na(tt_g$Genus), "Unassigned", tt_g$Genus))

meta$dominant_genus <- colnames(gen_mat)[apply(gen_mat[meta$sample_id, ], 1, which.max)]
meta$dominant_share <- apply(gen_mat[meta$sample_id, ], 1, max)

gen_rank   <- names(sort(colMeans(gen_mat), decreasing = TRUE))
top9       <- head(gen_rank[!grepl("^Unassigned", gen_rank)], 9)
genus_cols <- c(setNames(tol_muted, top9), Other = "#DDDDDD")

bar_df <- gen_mat %>%
  mutate(sample_id = rownames(.)) %>%
  pivot_longer(-sample_id, names_to = "Genus", values_to = "Abundance") %>%
  mutate(Genus = ifelse(Genus %in% top9, Genus, "Other")) %>%
  group_by(sample_id, Genus) %>% summarise(Abundance = sum(Abundance)) %>%
  left_join(meta %>% select(sample_id, inoculant), by = "sample_id") %>%
  mutate(Genus = factor(Genus, levels = c(top9, "Other")))

# ---- Fig S2 and Fig S3 ------------------------------------------------------
figS2 <- ggplot(bar_df, aes(x = sample_id, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~inoc_lab(inoculant), scales = "free_x", labeller = label_parsed) +
  scale_fill_manual(values = genus_cols) +
  scale_x_discrete(labels = function(x) sub("^.*_", "", x)) +
  ylab("Relative Abundance") + xlab("Replicate plant") +
  guides(fill = guide_legend(ncol = 1)) + th_bar
save_fig(figS2, "FigS2_genus_composition_per_plant", w = 10, h = 7)
add_legend("FigS2", "Genus-level composition of each plant",
  paste0("Relative abundance of the nine most abundant genera across all plants; all remaining genera are ",
         "grouped as Other (filtered and rarefied). Panels are inoculants: ", inoc_key(), "."))

figS3 <- bar_df %>% group_by(inoculant, Genus) %>%
  summarise(Abundance = sum(Abundance) / n_distinct(sample_id)) %>%
  ggplot(aes(x = inoculant, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = genus_cols) +
  ylab("Mean Relative Abundance") + scale_x_inoc() + xlab("Inoculum") +
  guides(fill = guide_legend(ncol = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right", legend.text = element_text(size = 8))
save_fig(figS3, "FigS3_genus_composition_mean_by_inoculant", w = 7.5, h = 5.5)
add_legend("FigS3", "Mean genus-level composition by inoculant",
  "Mean relative abundance of the nine most abundant genera per inoculant (same data as Fig S2).")

# ---- Fig S4: PCoA coloured by dominant genus ---------------------------------
dom_top <- intersect(names(genus_cols), unique(meta$dominant_genus))
pcoa_df$dominant_genus <- factor(ifelse(meta[rownames(pcoa_df), "dominant_genus"] %in% dom_top,
                                        meta[rownames(pcoa_df), "dominant_genus"], "Other"),
                                 levels = c(dom_top, "Other"))
figS4 <- ggplot(pcoa_df, aes(Axis.1, Axis.2, fill = dominant_genus, shape = persist_plant)) +
  geom_point(size = 4, colour = "grey20") +
  scale_fill_manual(values = genus_cols[c(dom_top, "Other")], name = "Dominant genus") +
  scale_shape_manual(values = c(Control = 21, `Not persistent` = 24, Persistent = 22), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21), ncol = 2)) +
  labs(x = paste0("PCoA1 (", eig[1], "%)"), y = paste0("PCoA2 (", eig[2], "%)")) +
  th_pcoa + theme(legend.box = "vertical")
save_fig(figS4, "FigS4_PCoA_by_dominant_genus", w = 8, h = 7)
add_legend("FigS4", "PCoA clusters track which genus dominates each plant, not persistence",
  paste0("Same ordination as Fig 2, with points filled by the most abundant genus in each plant and shaped ",
         "by plant-level persistence. Clusters correspond to Curtobacterium-dominated, Serratia-dominated and ",
         "more even communities; balcony block does not explain them (see Table S10)."))

# ---- Genus differences: SC13, SC35 and pooled vs controls (exploratory) -----
genus_da <- function(test_ids, label) {
  ctl  <- meta$sample_id[meta$inoculant == "MgCl2"]
  keep <- colMeans(gen_mat[c(test_ids, ctl), ]) > 0.005
  data.frame(comparison = label, genus = colnames(gen_mat)[keep],
             mean_test    = colMeans(gen_mat[test_ids, keep]),
             mean_control = colMeans(gen_mat[ctl, keep]),
             p = apply(gen_mat[, keep], 2, function(x) wilcox.test(x[test_ids], x[ctl], exact = FALSE)$p.value)) %>%
    mutate(log2FC = log2((mean_test + 1e-4) / (mean_control + 1e-4)),
           p_BH = p.adjust(p, "BH")) %>% arrange(p)
}
keep_tab(bind_rows(
  genus_da(meta$sample_id[meta$inoculant %in% c("SC13", "SC35")], "SC13+SC35 vs control"),
  genus_da(meta$sample_id[meta$inoculant == "SC13"], "SC13 vs control"),
  genus_da(meta$sample_id[meta$inoculant == "SC35"], "SC35 vs control")), "genus_differences")


# =============================================================================
# 9. TRACKING FOCAL STRAINS: ASVs IDENTICAL TO GENOME 16S COPIES
# =============================================================================
fa  <- readLines(iso_fasta)
hdr <- grep("^>", fa)
iso <- vapply(seq_along(hdr), function(i) {
  end <- if (i < length(hdr)) hdr[i + 1] - 1 else length(fa)
  toupper(paste(fa[(hdr[i] + 1):end], collapse = ""))
}, character(1))
names(iso) <- sub("^>(\\S+).*", "\\1", fa[hdr])          # e.g. SC35_copy3
iso_id <- sub("_copy.*", "", names(iso))
rc <- function(s) paste(rev(strsplit(chartr("ACGT", "TGCA", s), "")[[1]]), collapse = "")

asv_seqs  <- taxa_names(ps_u)
focal_map <- bind_rows(lapply(asv_seqs, function(a) {
  a_rc <- rc(a)
  hit  <- vapply(iso, function(s) grepl(a, s, fixed = TRUE) || grepl(a_rc, s, fixed = TRUE), logical(1))
  if (!any(hit)) return(NULL)
  data.frame(asv = a, isolate = unique(iso_id[hit]))
})) %>%
  group_by(asv) %>% mutate(shared_with = paste(sort(isolate), collapse = ";")) %>% ungroup() %>%
  left_join(as.data.frame(as(tax_table(ps_u), "matrix")) %>% mutate(asv = rownames(.)) %>%
              select(asv, Family, Genus), by = "asv")
keep_tab(focal_map %>% mutate(asv_length = nchar(asv)), "asv_matches")

inoculated <- c("SC01", "SC02", "SC03", "SC04", "SC09", "SC13", "SC25", "SC35")
rel_u <- as(otu_table(transform_sample_counts(ps_u, function(x) x / sum(x))), "matrix")
cnt_u <- as(otu_table(ps_u), "matrix")

focal_ra <- bind_rows(lapply(inoculated, function(sc) {
  a <- focal_map$asv[focal_map$isolate == sc]
  data.frame(isolate = sc, sample_id = rownames(rel_u),
             focal_reads = if (length(a)) rowSums(cnt_u[, a, drop = FALSE]) else 0,
             focal_RA    = if (length(a)) rowSums(rel_u[, a, drop = FALSE]) else 0)
})) %>%
  left_join(meta %>% select(sample_id, inoculant, T1, T2, T3, TopLeaf), by = "sample_id") %>%
  mutate(group = case_when(as.character(inoculant) == isolate ~ "Inoculated with this strain",
                           inoculant == "MgCl2"                ~ "MgCl2 control",
                           inoculant == "SCB8"                 ~ "SynCom",
                           TRUE                                ~ "Other single inoculant"))

focal_summary <- focal_ra %>%
  group_by(isolate, group) %>%
  summarise(n = n(), prevalence = mean(focal_reads > 0), mean_RA_pct = 100 * mean(focal_RA)) %>%
  ungroup() %>%
  left_join(bind_rows(lapply(inoculated, function(sc) {
    x <- focal_ra %>% filter(isolate == sc)
    data.frame(isolate = sc,
               p_inoculated_vs_control = tryCatch(suppressWarnings(
                 wilcox.test(x$focal_RA[x$group == "Inoculated with this strain"],
                             x$focal_RA[x$group == "MgCl2 control"], exact = FALSE)$p.value),
                 error = function(e) NA_real_))
  })), by = "isolate") %>%
  left_join(focal_map %>% filter(isolate %in% inoculated) %>% group_by(isolate) %>%
              summarise(n_ASVs = n_distinct(asv), shared_with = paste(unique(shared_with), collapse = " | ")),
            by = "isolate")
keep_tab(focal_summary, "focal_summary")
keep_tab(focal_ra %>% filter(group == "Inoculated with this strain") %>%
           transmute(isolate, sample_id, culture_T1T2T3 = paste0(T1, T2, T3),
                     culture_TopLeaf = TopLeaf, focal_reads, focal_RA_pct = 100 * focal_RA),
         "focal_vs_culture")

# ---- Fig 5 --------------------------------------------------------------------
focal_group_cols <- c("Inoculated with this strain" = okabe[["vermillion"]],
                      "MgCl2 control"               = okabe[["grey"]],
                      "Other single inoculant"      = okabe[["sky"]],
                      "SynCom"                      = okabe[["green"]])
fig5 <- focal_ra %>%
  filter(isolate != "SC03") %>%     # SC03's ASV is identical to the dominant environmental Curtobacterium ASV
  mutate(group = factor(group, levels = names(focal_group_cols))) %>%
  ggplot(aes(group, 100 * focal_RA + 0.001, fill = group)) +
  geom_jitter(shape = 21, colour = "black", width = 0.2, height = 0, size = 2) +
  facet_wrap(~inoc_lab(isolate), nrow = 1, labeller = label_parsed) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10), labels = c("0", "0.01", "0.1", "1", "10")) +
  scale_fill_manual(values = focal_group_cols, name = NULL) +
  xlab(NULL) + ylab("Focal ASV relative abundance (%)") +
  th_strip + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
save_fig(fig5, "Fig5_focal_ASV_tracking", w = 13, h = 4.5)
fs <- focal_summary %>% filter(group %in% c("Inoculated with this strain", "MgCl2 control")) %>%
  select(isolate, group, prevalence, mean_RA_pct) %>%
  pivot_wider(names_from = group, values_from = c(prevalence, mean_RA_pct))
add_legend("Fig5", "Exact-match focal ASVs in leaf communities at harvest (10 September 2024)",
  paste0("Relative abundance (filtered, unrarefied; log scale, 0 plotted at the axis floor) ",
         "of ASVs identical to a 16S rRNA gene copy from each focal strain's genome assembly. SC03 is omitted because ",
         "its ASV is identical to the dominant environmental Curtobacterium ASV; SC09 had no exact-match ASV. ",
         "Prevalence in inoculated plants vs controls: ",
         paste(sprintf("%s %.0f%% vs %.0f%%", fs$isolate,
                       100 * fs[["prevalence_Inoculated with this strain"]],
                       100 * fs[["prevalence_MgCl2 control"]]), collapse = "; "),
         ". Inoculants: ", inoc_key(drop_control = TRUE), "."))


# =============================================================================
# 10. SENSITIVITY ANALYSES AND ALTERNATIVE READ FILTERING
# =============================================================================
perm_one <- function(dist, data, rhs, label, blocks = NULL) {
  set.seed(seed)
  perms <- if (is.null(blocks)) nperm else how(blocks = blocks, nperm = nperm)
  a <- adonis2(as.formula(paste("dist ~", rhs)), data = data, permutations = perms)
  data.frame(test = label, adonis_row(a))
}
jac    <- phyloseq::distance(ps, method = "jaccard", binary = TRUE)
meta_u <- meta[sample_names(ps_u), ]
hel    <- vegdist(decostand(as(otu_table(ps_u), "matrix"), "hellinger"), "euclidean")
ps_no1335 <- subset_samples(ps, !inoculant %in% c("SC13", "SC35"))
m_no1335  <- data.frame(sample_data(ps_no1335))
ps_sc09   <- subset_samples(ps, inoculant %in% c("SC09", "MgCl2"))
m_sc09    <- data.frame(sample_data(ps_sc09)) %>% mutate(inoculant = droplevels(inoculant))

sensitivity <- bind_rows(
  perm_one(jac, meta, "persist_plant",     "Jaccard (presence/absence): plant-level persistence"),
  perm_one(jac, meta, "persist_species",   "Jaccard (presence/absence): species-level persistence"),
  perm_one(hel, meta_u, "persist_plant",   "Hellinger, unrarefied: plant-level persistence"),
  perm_one(hel, meta_u, "persist_species", "Hellinger, unrarefied: species-level persistence"),
  perm_one(phyloseq::distance(ps_no1335, "bray"), m_no1335, "persist_species",
           "Bray-Curtis: species-level persistence without SC13 and SC35"),
  perm_one(phyloseq::distance(ps_sc09, "bray"), m_sc09, "inoculant", "Bray-Curtis: SC09 vs control"),
  perm_one(bray, meta, "Block", "Bray-Curtis: balcony block"),
  perm_one(bray, meta, "persist_plant", "Bray-Curtis: plant-level persistence, permuted within block",
           blocks = meta$Block))

meta$insect_fraction <- with(filter_summary[match(meta$sample_id, filter_summary$sample_id), ],
                             (insect_reads + bee_gut_reads) / raw_reads)
ps_full_rare <- rarefy_even_depth(ps_full, sample.size = min(sample_sums(ps_full)),
                                  rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
alt_tests <- function(px, label) {
  d <- phyloseq::distance(px, method = "bray")
  m <- meta[sample_names(px), ]
  sub_test <- function(sc) {
    ids <- m$sample_id[m$inoculant %in% c(sc, "MgCl2")]
    mm  <- m[ids, ] %>% mutate(inoculant = droplevels(inoculant))
    perm_one(as.dist(as.matrix(d)[ids, ids]), mm, "inoculant", paste0(label, ": ", sc, " vs control"))
  }
  sh <- suppressWarnings(estimate_richness(px, measures = "Shannon"))
  sh$sample_id <- rownames(sh)
  sh <- left_join(sh, m %>% select(sample_id, persist_species), by = "sample_id")
  bind_rows(
    perm_one(d, m, "persist_plant",   paste0(label, ": plant-level persistence")),
    perm_one(d, m, "persist_T3",      paste0(label, ": detected at T3")),
    perm_one(d, m, "persist_species", paste0(label, ": species-level persistence")),
    sub_test("SC13"), sub_test("SC35"), sub_test("SC09"),
    data.frame(test = paste0(label, ": Shannon ~ species-level persistence (Kruskal-Wallis)"),
               F = NA, R2 = NA, p = kruskal.test(Shannon ~ persist_species, data = sh)$p.value))
}
alt_filtering <- bind_rows(
  data.frame(test = "Insect + bee-gut read fraction ~ plant-level persistence (Kruskal-Wallis)",
             F = NA, R2 = NA, p = kruskal.test(insect_fraction ~ persist_plant, data = meta)$p.value),
  data.frame(test = "Insect + bee-gut read fraction ~ inoculant (Kruskal-Wallis)",
             F = NA, R2 = NA, p = kruskal.test(insect_fraction ~ inoculant, data = meta)$p.value),
  alt_tests(ps_bee, "Bee-gut genera retained"),
  alt_tests(ps_full_rare, "Full community"))
keep_tab(bind_rows(sensitivity   %>% mutate(analysis = "Sensitivity (main filtering)"),
                   alt_filtering %>% mutate(analysis = "Alternative filtering")) %>%
           select(analysis, test, F, R2, p), "sensitivity")

# ---- Fig S8 -------------------------------------------------------------------
figS8a <- ggplot(meta, aes(x = inoculant, y = 100 * insect_fraction, fill = inoculant)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, size = 2) +
  scale_fill_manual(values = inoc_cols) +
  scale_x_inoc() + xlab("Inoculum") + ylab("Insect-associated reads (%)") + th_box
ord_full  <- ordinate(ps_full_rare, method = "PCoA", distance = phyloseq::distance(ps_full_rare, method = "bray"))
eig_full  <- round(100 * ord_full$values$Relative_eig[1:2], 1)
pcoa_full <- data.frame(ord_full$vectors[, 1:2], meta[rownames(ord_full$vectors), ])
save_fig(figS8a + plot_pcoa(pcoa_full, "persist_plant", eig_full) +
           plot_layout(widths = c(1, 1.2)) + plot_annotation(tag_levels = "A"),
         "FigS8_full_community_incl_insect_genera", w = 13, h = 6)
ft <- setNames(alt_filtering$p, alt_filtering$test)
add_legend("FigS8", "Insect-associated reads do not differ by treatment; full-community results are consistent",
  paste0("(A) Percentage of reads per plant from insect endosymbiont and bee-gut genera (",
         paste(c(insect_genera, bee_genera), collapse = ", "), "); Kruskal-Wallis by inoculant p = ",
         signif(ft[2], 2), ", by plant-level persistence p = ", signif(ft[1], 2), "; ",
         sum(meta$insect_fraction > 0.5), " of ", nrow(meta), " plants had >50% insect-associated reads. ",
         "(B) PCoA of Bray-Curtis dissimilarities on the full community (only blank ASVs removed) rarefied to ",
         min(sample_sums(ps_full)), " reads. Full-community PERMANOVA: plant-level persistence p = ",
         signif(ft["Full community: plant-level persistence"], 2), "; species-level p = ",
         signif(ft["Full community: species-level persistence"], 2), "; SC13 vs control p = ",
         signif(ft["Full community: SC13 vs control"], 2), "; SC35 vs control p = ",
         signif(ft["Full community: SC35 vs control"], 2), ". See Table S10 for the bee-gut-retained version."))


# =============================================================================
# 11. PLANT TRAITS AND FUNGAL PRESENCE ON PLATES
# =============================================================================
traits <- read_excel(traits_xlsx, sheet = 1, na = c("", "NA")) %>%
  rename(Open_flowers = `Yellow/open_flowers`) %>%
  mutate(across(c(Plant_Height_cm, Plant_Weight, Stem_Weight, Plant_And_Stem_Weight,
                  Flower_clusters, Open_flowers), ~ suppressWarnings(as.numeric(.x))),
         inoculant = factor(inoc_from_id(Bacterial_Inoculum_Replicate), levels = inoc_levels),
         Block = factor(Block))
trait_vars <- c(Plant_Height_cm = "Plant height (cm)", Plant_And_Stem_Weight = "Plant + stem fresh weight (g)",
                Flower_clusters = "Flower clusters (count)", Open_flowers = "Open flowers (count)")
trait_tests <- bind_rows(lapply(names(trait_vars), function(v) {
  d <- traits %>% filter(!is.na(.data[[v]]))
  a <- anova(lm(as.formula(paste(v, "~ inoculant + Block")), data = d))
  data.frame(trait = trait_vars[[v]], n = nrow(d), mean = mean(d[[v]]), sd = sd(d[[v]]),
             ANOVA_inoculant_F = a["inoculant", "F value"], ANOVA_inoculant_p = a["inoculant", "Pr(>F)"],
             ANOVA_block_p = a["Block", "Pr(>F)"],
             KruskalWallis_inoculant_p = kruskal.test(as.formula(paste(v, "~ inoculant")), data = d)$p.value)
}))
keep_tab(trait_tests, "plant_traits")

figS6 <- traits %>%
  pivot_longer(all_of(names(trait_vars)), names_to = "trait", values_to = "value") %>%
  mutate(trait = factor(trait_vars[trait], levels = trait_vars)) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = inoculant, y = value, fill = inoculant)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, size = 1.5) +
  facet_wrap(~trait, scales = "free_y") +
  scale_fill_manual(values = inoc_cols) +
  scale_x_inoc() + xlab("Inoculum") + ylab(NULL) + th_box + th_strip
save_fig(figS6, "FigS6_plant_traits", w = 9.5, h = 8)
tt <- setNames(trait_tests$ANOVA_inoculant_p, trait_tests$trait)
add_legend("FigS6", "Inoculation did not affect tomato growth or flowering",
  paste0("Plant height, combined plant and stem fresh weight, flower clusters and open flowers at the end of the ",
         "experiment (n = ", nrow(traits), " plants). Two-way ANOVA (inoculant + block), inoculant effect: ",
         paste(sprintf("%s p = %.2f", names(tt), tt), collapse = "; "), "."))

fung <- pers_raw %>%
  transmute(Isolate, Sampling_Timepoint, fungal = as.numeric(Fungal_PA),
            focal = as.numeric(Focal_Bacteria_Persistence)) %>%
  mutate(inoculant = factor(inoc_from_id(Isolate), levels = inoc_levels))
keep_tab(fung %>% group_by(inoculant, Sampling_Timepoint) %>%
           summarise(n = n(), fungal_prevalence = mean(fungal)) %>% ungroup(), "fungal_prevalence")

fung_plant <- fung %>% group_by(Isolate, inoculant) %>%
  summarise(n_fungal = sum(fungal), n_focal = sum(focal)) %>% ungroup()
fung_inoc  <- fung %>% filter(inoculant != "MgCl2") %>% droplevels()
fung_noSyn <- fung_inoc %>% filter(inoculant != "SCB8")
cmh <- function(d) mantelhaen.test(table(factor(d$focal, levels = c(1, 0)),
                                         factor(d$fungal, levels = c(1, 0)), d$Sampling_Timepoint))
cmh_all   <- cmh(fung_inoc)
cmh_noSyn <- cmh(fung_noSyn)
rho_fp <- suppressWarnings(cor.test(fung_plant$n_focal[fung_plant$inoculant != "MgCl2"],
                                    fung_plant$n_fungal[fung_plant$inoculant != "MgCl2"], method = "spearman"))
fung_tests <- data.frame(
  test = c("Fungal detections per plant (0-3) ~ inoculant (Kruskal-Wallis)",
           "Focal detection vs fungal presence, stratified by timepoint (CMH), inoculated plants",
           "Same, excluding SynCom plants",
           "Plant-level: number of focal detections vs number of fungal detections (Spearman)"),
  estimate = c(NA, cmh_all$estimate, cmh_noSyn$estimate, rho_fp$estimate),
  estimate_type = c(NA, "common odds ratio", "common odds ratio", "rho"),
  p = c(kruskal.test(n_fungal ~ inoculant, data = fung_plant)$p.value,
        cmh_all$p.value, cmh_noSyn$p.value, rho_fp$p.value))
keep_tab(fung_tests, "fungal_tests")

# ---- Fig S7 -------------------------------------------------------------------
figS7 <- fung_inoc %>%
  mutate(focal = factor(ifelse(focal == 1, "Focal strain detected", "Focal strain not detected"),
                        levels = c("Focal strain not detected", "Focal strain detected"))) %>%
  group_by(Sampling_Timepoint, focal) %>%
  summarise(n = n(), prev = mean(fungal), se = sqrt(prev * (1 - prev) / n)) %>%
  ggplot(aes(x = Sampling_Timepoint, y = 100 * prev, fill = focal)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85, colour = "grey20") +
  geom_errorbar(aes(ymin = 100 * pmax(prev - se, 0), ymax = 100 * pmin(prev + se, 1)),
                position = position_dodge(0.8), width = 0.2) +
  geom_text(aes(y = 2, label = paste0("n=", n)), position = position_dodge(0.8), vjust = 0, size = 3) +
  scale_fill_manual(values = c("Focal strain not detected" = okabe[["vermillion"]],
                               "Focal strain detected" = okabe[["blue"]]), name = NULL) +
  xlab("Sampling timepoint") + ylab("Leaflets with fungi on plates (%)") +
  theme(legend.position = "bottom")
save_fig(figS7, "FigS7_fungal_presence_vs_focal_detection", w = 6, h = 5)
add_legend("FigS7", "Fungi were less often recovered from leaflets where the focal strain was detected",
  paste0("Percentage of leaflet samples from inoculated plants with fungal growth on the plates, split by whether ",
         "the focal strain was re-isolated from the same sample (mean +/- binomial SE). Cochran-Mantel-Haenszel ",
         "test stratified by timepoint: common OR = ", sprintf("%.2f", cmh_all$estimate), ", p = ",
         signif(cmh_all$p.value, 2), "; excluding SynCom plants OR = ", sprintf("%.2f", cmh_noSyn$estimate),
         ", p = ", signif(cmh_noSyn$p.value, 2), ". Fungal detections per plant did not differ by inoculant ",
         "(Kruskal-Wallis p = ", signif(fung_tests$p[1], 2), ")."))


# =============================================================================
# 12. CULTURABLE BACTERIAL DENSITIES OVER TIME (all colonies, all 66 plants)
# =============================================================================
cfu <- pers_raw %>%
  transmute(Isolate, Sampling_Timepoint,
            KBHA     = suppressWarnings(as.numeric(CFU_g_KBHA)),
            KBHA_rif = suppressWarnings(as.numeric(CFU_g_rif)),
            focal    = as.numeric(Focal_Bacteria_Persistence)) %>%
  mutate(inoculant = factor(inoc_from_id(Isolate), levels = inoc_levels)) %>%
  pivot_longer(c(KBHA, KBHA_rif), names_to = "medium", values_to = "CFU_g") %>%
  mutate(medium = factor(medium, levels = c("KBHA", "KBHA_rif"),
                         labels = c("'KBHA (all colonies)'", "'KBHA-rif (all colonies)'")))
keep_tab(cfu %>%
           group_by(inoculant, Sampling_Timepoint, medium) %>%
           summarise(n = sum(!is.na(CFU_g)), median_CFU_g = median(CFU_g, na.rm = TRUE),
                     min_CFU_g = min(CFU_g, na.rm = TRUE), max_CFU_g = max(CFU_g, na.rm = TRUE),
                     n_zero = sum(CFU_g == 0, na.rm = TRUE)) %>%
           ungroup() %>% mutate(medium = gsub("'", "", medium)), "cfu_summary")

figS9 <- ggplot(cfu, aes(Sampling_Timepoint, CFU_g + 1)) +
  geom_line(aes(group = Isolate), colour = "grey65", linewidth = 0.4) +
  geom_point(aes(fill = factor(focal, levels = c(0, 1))), shape = 21, size = 1.8, colour = "grey20") +
  facet_grid(medium ~ inoc_lab(inoculant), labeller = label_parsed) +
  scale_y_log10(breaks = 10^c(0, 2, 4, 6, 8, 10),
                labels = c(expression(0), expression(10^2), expression(10^4),
                           expression(10^6), expression(10^8), expression(10^10))) +
  scale_fill_manual(values = c("0" = "#FFFFFF", "1" = okabe[["blue"]]),
                    labels = c("0" = "Not detected", "1" = "Detected"),
                    name = "Focal strain on KBHA-rif", drop = FALSE) +
  xlab("Sampling timepoint") + ylab("CFU per g leaf tissue") +
  th_strip + theme(strip.text.x = element_text(size = 7), legend.position = "bottom")
save_fig(figS9, "FigS9_CFU_per_plant_over_time", w = 14, h = 6)
add_legend("FigS9", "Culturable bacterial densities on inoculated leaflets over time",
  paste0("Colony-forming units per gram of leaf tissue from one lower leaflet per plant at T1 (2024-08-10), ",
         "T2 (2024-08-25) and T3 (2024-09-08), on KBHA (top) and KBHA-rifampicin (bottom). All colonies were ",
         "counted regardless of morphology, so counts include non-focal bacteria. Lines join samples from the ",
         "same plant (n = ", n_distinct(cfu$Isolate), " plants); filled points mark samples where the focal strain ",
         "was recovered. Zero counts are plotted at the axis floor. Inoculants: ", inoc_key(), "."))


# =============================================================================
# 13. SUPPLEMENTARY TABLES S1-S13, FIGURE LEGENDS, SESSION INFO
# =============================================================================
# Per-plant table: metadata, persistence, distance to controls, dominant genus and alpha diversity
keep_tab(meta %>%
           left_join(alpha %>% select(sample_id, Observed, Shannon, max_ASV_share), by = "sample_id") %>%
           relocate(sample_id, inoculant, species, replicate, Block),
         "plant_table")

# Table number -> blocks (results name, label shown above the block; NA = no label)
supp_layout <- list(
  S1  = list(c("cfu_summary", NA)),
  S2  = list(c("plant_table", NA)),
  S3  = list(c("read_filtering", NA)),
  S4  = list(c("asv_matches", NA)),
  S5  = list(c("plant_traits", NA)),
  S6  = list(c("topleaf_by_inoculant", "A. Focal-strain detection on the T3 lower leaflet and the top leaf, by inoculant"),
             c("topleaf_test",         "B. Exact McNemar test, T3 lower leaflet vs top leaf")),
  S7  = list(c("group_counts",          "A. Number of plants per inoculant in each plant-level persistence group"),
             c("permanova_by_scheme",   "B. PERMANOVA (overall and pairwise) and dispersion tests by persistence scheme"),
             c("dissimilarity_summary", "C. Within-group dissimilarity and dissimilarity to MgCl2 controls by persistence group"),
             c("distance_tests",        "D. Distance-to-control tests, correlations, and PERMANOVA within and among inoculants")),
  S8  = list(c("each_vs_control",   "A. Each inoculant vs MgCl2 controls: composition, dispersion and Shannon diversity"),
             c("pooled_vs_control", "B. SC13 and SC35 pooled vs MgCl2 controls"),
             c("alpha_tests",       "C. Alpha diversity tests by persistence scheme and inoculant")),
  S9  = list(c("genus_differences", NA)),
  S10 = list(c("sensitivity", NA)),
  S11 = list(c("focal_summary", NA)),
  S12 = list(c("focal_vs_culture", NA)),
  S13 = list(c("fungal_prevalence", "A. Prevalence of fungal growth on plates by inoculant and timepoint"),
             c("fungal_tests",      "B. Tests of fungal presence")))
stopifnot(all(unlist(lapply(supp_layout, function(b) sapply(b, `[`, 1))) %in% names(tabs)))

scheme_labels <- c(persist_plant = "Plant level (detected at ≥2 of 3 timepoints)", persist_T3 = "Detected at T3",
                   persist_species = "Species level (SC09, SC13, SC35 = non-persisters)", inoculant = "Inoculum identity")
species_lookup <- c(inoc_species, "All inoculated" = "All inoculated plants")
prep_tab <- function(x) {
  if ("scheme" %in% names(x)) x$scheme <- ifelse(x$scheme %in% names(scheme_labels), scheme_labels[x$scheme], x$scheme)
  for (col in c("inoculant", "isolate")) {
    if (col %in% names(x) && !"species" %in% names(x)) {
      sp <- unname(species_lookup[as.character(x[[col]])]); sp[is.na(sp)] <- ""
      x <- tibble::add_column(x, species = sp, .after = col)
    }
  }
  x %>% mutate(across(where(is.factor), as.character),
               across(where(is.double), ~ ifelse(abs(.x) >= 1000, round(.x), signif(.x, 4))))
}

wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 11)
st_head  <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9E2EC")
st_label <- openxlsx::createStyle(textDecoration = c("bold", "italic"))
for (tn in names(supp_layout)) {
  sh <- paste("Table", tn); openxlsx::addWorksheet(wb, sh)
  row <- 1; widths <- numeric(0)
  for (b in seq_along(supp_layout[[tn]])) {
    blk <- supp_layout[[tn]][[b]]; dat <- prep_tab(tabs[[blk[1]]])
    csv_name <- paste0(tn, if (length(supp_layout[[tn]]) > 1) LETTERS[b] else "", "_", blk[1], ".csv")
    write.csv(dat, file.path(tab_dir, csv_name), row.names = FALSE, na = "")
    if (!is.na(blk[2])) {
      openxlsx::writeData(wb, sh, blk[2], startRow = row)
      openxlsx::addStyle(wb, sh, st_label, rows = row, cols = 1)
      row <- row + 1
    }
    openxlsx::writeData(wb, sh, dat, startRow = row, headerStyle = st_head, keepNA = FALSE)
    w <- pmin(pmax(nchar(names(dat)) + 2,
                   vapply(dat, function(v) max(nchar(as.character(v)), 0, na.rm = TRUE) + 2, numeric(1))), 60)
    widths <- c(pmax(widths[seq_along(w)], w, na.rm = TRUE), widths[-seq_along(w)])
    row <- row + nrow(dat) + 2
  }
  openxlsx::setColWidths(wb, sh, cols = seq_along(widths), widths = widths)
  if (length(supp_layout[[tn]]) == 1) openxlsx::freezePane(wb, sh, firstActiveRow = 2)
}
openxlsx::saveWorkbook(wb, file.path(out_dir, "Supplementary_Tables_S1-S13.xlsx"), overwrite = TRUE)

legend_order <- c("Fig1", "Fig2", "Fig3", "Fig4", "Fig5", paste0("FigS", 2:10))
writeLines(c("# Figure titles and legends", "",
             paste0("Generated ", Sys.Date(), " (", nperm, " permutations). Statistics are filled in from this run."), "",
             "## Data processing", "", paste("-", filter_notes), "",
             unlist(lapply(intersect(legend_order, names(fig_legends)), function(f)
               c(paste0("## ", f), "", fig_legends[[f]], "")))),
           file.path(out_dir, "figure_legends.md"))
saveRDS(ps,   file.path(out_dir, "ps_balcony_filtered_rarefied.rds"))
saveRDS(ps_u, file.path(out_dir, "ps_balcony_filtered_unrarefied.rds"))
writeLines(c(paste("Run:", Sys.time()), paste("Permutations:", nperm),
             paste("Rarefaction depth:", rare_depth), "", capture.output(sessionInfo())),
           file.path(out_dir, "session_info.txt"))
message("Done. Figures, tables and legends written to ", normalizePath(out_dir))
