#!/bin/bash

# Define the list of SRP IDs and group Comparison
SRP_LIST="SRP_list.txt"
Comparison="hcc_v_nonhcc"

# Loop through each SRP ID and run the rMATS preparation script in the background
while read -r SRP; do
  nohup Rscript scripts/rmats_prep_post.R "$SRP" paired Comparison > "${SRP}.log" 2>&1 &
done < "$SRP_LIST"

# Run the other rmats processing scripts
Rscript scripts/impute_reformat.R -s SRP_list.txt -m metadata/metadata.txt -c Comparison -g data/Homo_sapiens.GRCh38.103.gtf
Rscript scripts/JC.R Comparison
Rscript scripts/add_gene_names.R Comparison
Rscript scripts/rmats_tables.R Comparison
