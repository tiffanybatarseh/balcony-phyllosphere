# balcony-phyllosphere
Tracking single-strain bacterial inoculants and 16S community assembly on outdoor-grown tomato leaves (R, DADA2, phyloseq)

# Inoculant persistence and community legacy in the tomato phyllosphere

Gnotobiotically germinated tomato seedlings (cv. Moneymaker) were dip-inoculated with one of eight rifampicin-resistant leaf bacteria, all eight together (SynCom), or buffer (MgCl₂), and grown outdoors. Focal strains were tracked by plating leaflets at three timepoints and the top leaf, and whole-shoot bacterial communities were profiled at harvest by 16S rRNA gene (V5–V7) amplicon sequencing.

## Contents

| File | Description |
|---|---|
| `Balcony_analysis.R` | Single script: (optional) DADA2, read filtering, all statistics, Figures 1–5 and S2–S10, Supplementary Tables S1–S13 |
| `data/balcony_asv_table.rds` | ASV count table (66 plants + 2 extraction blanks) and SILVA v138.1 taxonomy |
| `data/2024-10-15_Bacteria_Persistence.xlsx` | Leaflet plating at T1–T3: CFU g⁻¹, focal-strain detection, fungal presence |
| `data/2026-10-28_Balcony_TopLeaf_scoring.xlsx` | Focal-strain detection on the top (youngest) leaf |
| `data/2024-10-15_PlantTraits.xlsx` | Plant height, fresh weight and flower counts at harvest |
| `data/isolate_16S_all_copies.fasta` | All 16S rRNA gene copies from the SynCom hybrid genome assemblies |

Raw reads: NCBI BioProject [PRJNA1556283](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1556283).

## Running the analysis

To rebuild the ASV table from the SRA reads:

1. Remove the 5′ primers with cutadapt v4.9. The original command was not retained; the command below reproduces the trimmed files (primers removed from the 5′ end of each read, no quality trimming, reads without a primer kept):

   ```bash
   for r1 in *_R1_001.fastq.gz; do
     s=${r1%_R1_001.fastq.gz}
     cutadapt -g AACMGGATTAGATACCCKG -G ACGTCATCCCCACCTTCC \
              -o trimmed/trimmed_${s}_R1.fastq.gz -p trimmed/trimmed_${s}_R2.fastq.gz \
              ${s}_R1_001.fastq.gz ${s}_R2_001.fastq.gz
   done
   ```

2. Download the SILVA v138.1 training set (`silva_nr99_v138.1_train_set.fa.gz`) to `data/`.
3. Set `run_dada2 <- TRUE` in `Balcony_analysis.R` and run it.

## Citation

If you use this code, please cite the paper.

## License

MIT
