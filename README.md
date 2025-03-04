# HCC_PBMC_ASEs

This repository contains scripts used in the manuscript "Alternative Splicing Immune Signature for Detecting HBV-Related Hepatocellular Carcinoma in Peripheral Blood Mononuclear Cells." Raw data can be found in the manuscript itself or in the following SRP datasets: SRP094502, SRP162958, SRP201023, SRP312121, SRP328165, SRP330954, SRP332585, and SRP412160.

Directory Structure:

•	Figures_Analysis: Scripts for analysis, correlations, regressions, ontology analysis, exon annotation, gene expression meta-analysis, and visualization.
•	rMATS: Scripts for alternative splicing analysis.
o	Before running, ensure you generate a text file listing all SRP datasets used. Additionally, create text files control.txt and condition.txt for each SRP dataset, containing BAM file locations for control and condition samples.
o	Execute the bash script run_rmats.sh for processing.
•	RF: Scripts and visualizations for predictive modeling.
