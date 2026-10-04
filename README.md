# H5Nx synonymous genomic architecture — analysis code and data

This repository contains the custom scripts and derived data used in:

> Teng Long, Rachel H.H. Ching, Kenrie P.Y. Hui, J.S. Malik Peiris, Michael C.W. Chan.
> **Synonymous Genomic Architecture of H5Nx Influenza Virus Reveals Potential of
> Systemic Dissemination and Tropism.** *Virus Evolution* (submitted).

## Overview

The study analyzes 2,472 complete H5Nx influenza A virus genomes (H5N1, H5N2, H5N5,
H5N6, H5N8 and minor subtypes) circulating in East Asia, combining phylogenetic
reconstruction with codon usage analyses (correspondence analysis, RSCU, CAI, ENC,
neutrality plots) to characterize the synonymous genomic architecture of H5Nx viruses.

## Repository structure

```
H5Nx_codon_architecture/
├── README.md
├── 1_data/
│   ├── CAI_results/                   # per-gene CAI values (metadata_CAI.csv and relatives)
│   ├── RSCU_output/                   # RSCU output per segment x subtype
│   ├── CodonW_output/                 # Raw data of CodonW output files of each subtypes and segments
│   ├── ENC_data/                      # CodonW ENC output
│   ├── HA_structure/                  # Analysis of arginine abundance in H5 HA
│   ├── Genomic_diversity/             # Files for nucleotide diversity analysis         
│   ├── Neutrality_results/            # GC123 tables from 3_calculate_codon_GC123.py
│   └── Trees/                         # HA and NA ML trees (.nwk, IQ-TREE, 1,000 UFBoot)
└── 2_scripts/
    ├── 01_data_prep/
    │   ├── NCBI_sequence_prep/               # only for NCBI downloaded data, N/A for GISAID originated data
    │   │   ├── 1_split_meta_fasta.py         # Split + rename headers to >isolate|assembly (GCA_ only)
    │   │   ├── 2_filter_set_assemblies.py    # Redundant safety check: drop any residual 'set:' entries
    │   │   └── 3_rna2dna.py                  # Transfer rna to dna
    │   ├── neutrality/
    │   │   ├── 1_Fasta_remove_terminal_codon.py         # Remove the terminal codons for analysis
    │   │   ├── 2_Fasta_remove_base.py                   # Remove ATG(T), TGG(W), ATT(I), ATC(I), ATA(I)
    │   │   └── 3_calculate_codon_GC123.py               # Calculate and generate GC123 tables   
    │   └── structure_analysis/         
    │       └── analyze_h5_subtypes.py                   # Calculate the arginine abundance in H5 HA  
    └── 02_figures/
        ├── fig1_phylo_diversity.R                  # Figure 1A-B
        ├── fig2_correspondence_analysis.R          # Figure 2A-C
        ├── fig3_rscu_ha.R                          # Figure 3A-B
        ├── fig4_cai.R                              # Figure 4A-B
        ├── fig5_enc.R                              # Figure 5 and Figure S2A
        ├── fig6_na.R                               # Figure 6 and Figure S2B
        └── figS1_cai.R                             # Figure S1
        
```

## Requirements

- Python >= 3.9, `biopython >= 1.80`
- R >= 4.5.2, packages: `ggplot2`, `ggpubr`, `rstatix`, `dplyr`, `ggtree`, `phytools`, `phangorn`, `pheatmap`, `plotly`
- External tools: `MAFFT v7.49`, `IQ-TREE v3.0.1`, `CodonW v1.4.4`, `Dedupe` (Geneious Prime v2025.2), `CAIcal` (<https://ppuigbo.me/programs/CAIcal/>)


## Data availability notes

- **Sequence data are NOT redistributed in this repository** in compliance with the
  GISAID Terms of Use. All analyzed sequences are available from the GISAID EpiFlu
  Database and the NCBI Influenza Virus Resource under the accession numbers listed
  in Supplementary Table S1 of the paper.
- We gratefully acknowledge all data contributors — the originating and submitting
  laboratories responsible for obtaining the specimens and generating the sequences
  shared via the GISAID and NCBI Initiative.
- Files under `data/` are **derived analysis results** (RSCU matrices, CAI values,
  GC statistics, tree files), which are provided here to ensure full reproducibility
  of every figure and table.


## Contact

Teng Long — <longt@hku.hk>
