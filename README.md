# LRRFIP1 Analysis Code

## Repository structure

```
burden_analysis/          Rare variant burden testing (UK Biobank WGS)
  regenie_pipeline.sh       Regenie SKAT-O pipeline (step 1 + step 2)
  build_annotations.py      Annotation, mask, and set-list construction
  run_vep.sh                VEP variant annotation (Docker)
  figures/                  Manuscript figure scripts (forest plot, VEP, directionality)

isoform_shortread/        Short-read isoform quantification (Salmon v1.10.3)
  fibroblast_composite.py   Main figure: fibroblast LRRFIP1 isoform fractions + PCA
  ipsc_composite.py         Supplementary figure: iPSC isoform fractions + PCA
  supp_table_ipsc.py        Supplementary table: iPSC isoform data
  linked_genes_figures.py   Linked genes expression analysis
  top100_and_heatmap.py     Top 100 expressed genes and heatmap
  figure_data_tables.py     Data tables underlying main figures
  regenie_pertest_table.py  Supplementary Regenie per-test breakdown
  NR3C1_expression.py       NR3C1 expression analysis
  pairwise_delta_shortread.py  Supplementary: short-read pairwise isoform delta
  pairwise_delta_combined.py   Main figure: combined short- and long-read pairwise delta

isoform_longread/         Long-read isoform analysis (PacBio IsoSeq)
  pigeon_run.sh             Pigeon transcript classification (GENCODE v49)
  run_isoquant_flnc.sh      IsoQuant FLNC transcript assignment
  pairwise_delta.py         Pairwise isoform fraction delta (long-read only)
  isoquant_lrrfip1_composite.py   IsoQuant LRRFIP1 composite figure
  flnc_isoform_distribution.py    FLNC isoform distribution
  lrrfip1_isoseq_isoform_distribution.py  IsoSeq isoform distribution
  lrrfip1_6_reference_isoforms_mapping.py  Reference isoform mapping

motif_analysis/           Transcription factor binding site disruption (motifbreakR)
  run_motifbreakr_hocomoco.R  motifbreakR analysis with HOCOMOCO motifs (per-SNP PDFs + CSVs)
  motif_summary_tables.R      Post-processing: boolean matrix + counts tables (TSV + Excel)
```

## Dependencies

### Python (>= 3.10)
pandas, numpy, matplotlib, seaborn, scikit-learn, scipy, openpyxl

### R (Bioconductor)
motifbreakR, BSgenome.Hsapiens.UCSC.hg38, TFBSTools, seqLogo, GenomicRanges, data.table, gridExtra, openxlsx

### External tools
- [Regenie](https://rgcgithub.github.io/regenie/) v3.4+
- [Salmon](https://combine-lab.github.io/salmon/) v1.10.3
- [Ensembl VEP](https://www.ensembl.org/vep) (Docker image)
- [PacBio IsoSeq / Pigeon](https://isoseq.how/)
- [IsoQuant](https://github.com/ablab/IsoQuant) v3.4+

### Reference annotation
- GENCODE v49 (GRCh38)

## Data availability

[Citation of the paper here]

## Notes

- All scripts use relative paths (`./data/...`). Adjust the `data/` directory to point to your local data.
- Salmon quantification used a GENCODE v49 full-transcriptome index with decoy sequences.
- Isoform numbering (1-6) refers to GENCODE v49 LRRFIP1 transcripts (ENST IDs listed in supplementary tables), not RefSeq nomenclature.
- The Regenie pipeline was executed on the UK Biobank Research Analysis Platform (DNAnexus).

## Contact

Carla Bello — carla.bello@medgen.uzh.ch
