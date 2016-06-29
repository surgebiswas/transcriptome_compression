# Transcriptome compression
### Analysis code for results reported in:
Surojit Biswas, Konstantin Kerner, Paulo Jose Pereira Lima Teixeira, Jeffery L. Dangl, Vladimir Jojic, Philip A. Wigge. "Tradict enables high fidelity reconstruction of the euykaryotic transcriptome from 100 marker genes." Submitted (2016). BioRxiv: http://biorxiv.org/content/early/2016/06/21/060111.

With the exception of data, this repository contains everything required to perform downloads from the SRA and recapitulate analyses performed in our paper. The repository is composed of several directories. Each is listed below with a  description of its purpose and contents. 

# Repository directories
1. **data_download** - Contains the custom scripts necessary to download and quantify raw RNA-Seq reads from NCBI's Short Read Archive (SRA). This repository contains it's own README that describes how to perform the data downloads. An overview of the download workflow is also given in the Materials and Methods section in the Supplemental Information of our paper. 
2. **util** - Contains accessory code used by analysis scripts. Code in this directory is commented as needed to explain function and usage. 
3. **analysis** -- Contains all code used to perform analysis and generate figures. The main analysis 'master scripts' are in the `Athaliana` and `Mmusculus` subfolders and have file names of the following form: `analyze_X_transcriptome_compression.m`.  These master scripts contain blocks of code, demarcated by `if-else` blocks, that each perform an analysis of the data or produce sub-panels of a figure. These scripts are commented to provide descriptions of what they are doing, as are all dependency scripts. 


Send questions/comments to surojitbiswas@g.harvard.edu.
