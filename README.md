# balcony-phyllosphere
Tracking single-strain bacterial inoculants and 16S community assembly on outdoor-grown tomato leaves (R, DADA2, phyloseq)

# Inoculant persistence and community legacy in the tomato phyllosphere

Analysis code for Batarseh et al., *FEMS Microbiology Letters* (submitted).

Gnotobiotically germinated tomato seedlings (cv. Moneymaker) were dip-inoculated with one of eight rifampicin-resistant leaf bacteria, all eight together (SynCom), or buffer (MgCl₂), and grown outdoors. Focal strains were tracked by plating leaflets at three timepoints and the top leaf, and whole-shoot bacterial communities were profiled at harvest by 16S rRNA gene (V5–V7) amplicon sequencing.

## Contents

| File | Description |
|---|---|
| `Balcony_analysis.R` | Single script: (optional) DADA2, read filtering, all statistics, Figures 1–5 and S2–S9, Supplementary Tables S1–S13 |
| `data/balcony_asv_table.rds` | ASV count table (66 plants + 2 extraction blanks) and SILVA v138.1 taxonomy |
| `data/2024-10-15_Bacteria_Persistence.xlsx` | Leaflet plating at T1–T3: CFU g⁻¹, focal-strain detection, fungal presence |
| `data/2026-10-28_Balcony_TopLeaf_scoring.xlsx` | Focal-strain detection on the top (youngest) leaf |
| `data/2024-10-15_PlantTraits.xlsx` | Plant height, fresh weight and flower counts at harvest |
| `data/isolate_16S_all_copies.fasta` | All 16S rRNA gene copies from the SynCom hybrid genome assemblies |

Raw reads: NCBI BioProject [PRJNA1556283](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1556283).

## Running the analysis

1. Install R (≥ 4.4) and the packages: `phyloseq` (Bioconductor), `vegan`, `permute`, `readxl`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `openxlsx`, and `pairwiseAdonis` (`devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")`).
2. Open R in this folder and run `source("Balcony_analysis.R")`.
3. Results are written to `output/`: figures (PDF and PNG), `Supplementary_Tables_S1-S13.xlsx`, one CSV per table, `figure_legends.md` (legends with the statistics from the run) and `session_info.txt`.

To rebuild the ASV table from the SRA reads:

1. Remove the anchored 5′ primers with cutadapt v4.9 (run in the folder with the raw FASTQ files), then move the `trimmed_*` files into `trimmed/`:

   ```bash
   for f in *_R1_001.fastq.gz; do
       base=$(basename "$f" _R1_001.fastq.gz)
       cutadapt \
           -g ^AACMGGATTAGATACCCKG \
           -G ^ACGTCATCCCCACCTTCC \
           -o trimmed_${base}_R1.fastq.gz \
           -p trimmed_${base}_R2.fastq.gz \
           ${base}_R1_001.fastq.gz ${base}_R2_001.fastq.gz \
           --minimum-length 200 \
           --discard-untrimmed \
           --quiet
   done
   ```

2. Download the SILVA v138.1 training set (`silva_nr99_v138.1_train_set.fa.gz`) to `data/`.
3. Set `run_dada2 <- TRUE` in `Balcony_analysis.R` and run it.

The published ASV table was produced with DADA2 v1.22.0 (in R v4.1.2) from a sequencing run shared with another study, so a rerun on these samples alone, or with a newer DADA2 version, can give very slightly different ASVs. All downstream analyses were run in R v4.6.1 (see `output/session_info.txt` after a run).

## Citation

If you use this code, please cite the paper.

## License

MIT
