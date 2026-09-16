cd ./data/pacbio/pigeon_inputs

# Step 1: Prepare v49 annotation + reference files
pigeon prepare \
  gencode.v49.annotation.gtf \
  human_GRCh38_no_alt_analysis_set.fasta \
  refTSS_v3.3_human_coordinate.hg38.sorted.bed \
  "intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv"

# Step 2: Prepare collapsed GFF
pigeon prepare collapse_isoforms.gff

# Step 3: Classify against v49
pigeon classify \
  collapse_isoforms.sorted.gff \
  gencode.v49.annotation.sorted.gtf \
  human_GRCh38_no_alt_analysis_set.fasta \
  --cage-peak refTSS_v3.3_human_coordinate.hg38.sorted.bed \
  --poly-a polyA.list.txt \
  --fl collapse_isoforms.flnc_count.txt \
  -d ../outputs

# Step 4: Filter
pigeon filter \
  ../outputs/collapse_isoforms.sorted_classification.txt \
  --isoforms collapse_isoforms.sorted.gff \
  -d ../outputs
