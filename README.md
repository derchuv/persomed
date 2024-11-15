# Immune Modules to guide Diagnosis and Personalized Treatment of Inflammatory Skin Diseases

This repository contains R scripts for generating graphs for the paper **"Awaiting for DOI number"**.

## Repository Structure

- `data/`: Contains raw and processed data.
- `scripts/`: Contains R scripts for data processing and analysis.
- `output/`: Contains generated graphs and tables.

## How to Run

1. Clone the repository:
   ```bash
   git clone git@github.com:derchuv/persomed.git
   ```
1. Move to the repository:
   ```bash
   cd persomed/data
   ```
1. Get RCC files, store and name them correctly:
   ```bash
   wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE280nnn/GSE280220/suppl/GSE280220_RAW.tar
   mkdir rcc
   tar -xvf GSE280220_RAW.tar -C rcc
   rm GSE280200_RAW.tar
   cd rcc
   for file in *.gz; do
      gzip -d "$file"
      new_filename="${file%.gz}"
      mv "$new_filename" "${new_filename#*_}"
   done
   ```
1. Run the main script that creates all pdf's used for the figure.
1. Check the results in the output folder.