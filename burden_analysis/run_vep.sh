docker run \
    -v /path/to/project/vep_analysis:/vep/input \
    -v /path/to/project/vep_analysis/vep_cache:/opt/vep/.vep \
    ensemblorg/ensembl-vep \
    vep \
    --input_file /vep/input/6genes_sites_clean.vcf.gz \
    --output_file /vep/input/vep_6genes_results.txt \
    --cache \
    --offline \
    --assembly GRCh38 \
    --everything \
    --tab \
    --fork 6 \
    --force_overwrite \
    --no_stats
