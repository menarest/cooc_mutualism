[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10931679.svg)](https://doi.org/10.5281/zenodo.10931679)

# Supplementary code for analyses of Menares et al. (in review), "Co-occurrence patterns do not predict mutualistic interactions between plant and butterfly species".

This repository contains the analysis done for Menares et al. (in review), which compares the accuracy of existing co-occurrence methods for estimating pairwise species interactions between plants and butterfly species from observational occurrence (P/A and abundance) data against species interactions estimated experimentally from long-term flower visitation data. 

# Requirements: 

To run all analyses from scratch including data wrangling and prepararion steps, you will need to download the data from BExIS using the IDs: 

- plants: 23586 (2)
- butterflies: 12526 (2)
- species taxonomy table: 31733 (8)
- species interactions:
  - Imago ALB: 31734 (8)
  - Larva ALB: 31735 (4)
  - Imago SCH: 31736 (3)
  - Larva SCH: 31737 (3)
- flower availability: 4302 (2) and 4964 (2)
- flower visitation: 10160 (3)

> [!NOTE]
> check that the datasets are named the same as in the scripts so you can run the analyses with no problem. 

If you prefer to run only the accuracy and data validation analyses use this dataset from BExIS: 31738 (3) #TODO add DOI once available

# Steps: 
1. run the [data_tidy script](scripts/wrangling/data_tidying.R) to clean and tidy the datasets
2. continue with the pairs files in order (pairs_1.., 2, 3, and 4)
3. run the files in the [scripts/analysis folder](scripts/analysis) in this order:
   - accuracy_co-occurrence
   - modeling_trophint_cooccur
   - troph_data_validation
   - matrices_cooccur

This should generate all files and plots of the paper and place them into the already folders.
