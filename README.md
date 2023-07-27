# data-and-code-for-fps

there are more data that oversized, please check on this site for the data[data link](https://drive.google.com/drive/folders/1Qu7oPUNQkYPxbC_eWwcLygtuGywIMhM1?usp=sharing)

# Dataset Description


## aver_pheno_snpspec_xy.rds

**Description:** This dataset includes averaged phenotypic and spectral data for individual trees in the plantation of 2021. It is stored in the RDS format.

## cor_allmonth.rds

**Description:** This dataset provides It's a list file containing gs models made with different data, for each month, and then run 100 times  in the study. It is stored in the RDS format.

## fam_location_all_trees.csv

**Description:** This CSV file contains family and location information for all individual trees in the plantation. It includes relevant attributes such as family name, site, and block.

## las_list.rds

**Description:** The 'las_list.rds' dataset comprises a list of LiDAR point cloud data for the entire plantation for the whole year of 2021. It is stored in the RDS format.

## rastimagelist.rds

**Description:** The 'rastimagelist.rds' dataset consists of a list of raster images for various spectral bands or vegetation indices. It is stored in the RDS format.
## code File Descriptions
1. **GWAS CODE.R**

   **Description:** The `GWAS CODE.R` file contains code for performing Genome-Wide Association Studies (GWAS) on the provided dataset using the VIs as PBWAS. GWAS is a statistical method used to identify genetic variants associated with particular traits or phenotypes. The code involves data preprocessing and possibly statistical analysis of significant genetic markers. The output may include p-values, effect sizes, and significant genetic loci associated with specific traits.

2. **data extraction.R**

   **Description:** The `data extraction.R` file includes code for extracting and preprocessing specific data from the provided datasets. It involve filtering, cleaning, and aggregating data to create subsets or summaries required for further analysis. This script is crucial for preparing the data for subsequent tasks such as PBWAS or correlation analysis.

3. **correlation of PS selection function.R**

   **Description:** The `correlation of PS selection function.R` file contains code for computing PS using the GS methods.
   

For any code and data related issues, please contact aj7105@gmail.com 


