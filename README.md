# APM_Sarcoma
Analysis of APM (antigen presenting and processing machinery) presence in Sarcoma

# Introduction

Sarcomas, a diverse range of malignant tumors, pose challenges due to limited treatment options and unfavorable outcomes. The antigen processing and presenting machinery (APM) is pivotal in creating tumor-specific antigens for immune recognition and targeting. Yet, the APM status in sarcomas remains poorly understood.

## Material and Methods

## Results
Our study uncovered APM defects in all the investigated sarcomas, with HLA Class I β2-microglobulin and HLA Class II being the least affected. Noteworthy, histotype-specific defects include LMP10 in LMS, MLPS, and DDLPS, and TAP2 in UPS, GIST, and CH. Spatial variations in APM expression were observed, with high-grade areas exhibiting distinct patterns. Patient-level heterogeneity also was noted, with certain individuals prone to APM component deficiencies. Importantly, loss of any APM component predicted distant metastasis in LMS and DDLPS, and overall survival in LMS.

## Conclusion
Study reveals pronounced defects in Antigen Presenting Machinery (APM) components within sarcomas, showcasing notable variations across histotypes and tumoral regions. Key alterations involve HLA Class I subunit β2-microglobulin, HLA Class I subunit α (HC10), and MHC I transporting unit TAP2. Notably, the absence of APM components serves as a prognostic indicator for distant metastases and overall survival in LMS and DDLPS. This underscores the clinical significance of APM dysregulation in these sarcoma subtypes, providing crucial insights into molecular mechanisms and informing personalized therapeutic strategies.

## Repository Structure

`Survival Analysis.R` This file contains the Kaplan Meier Survival Analysis. It produces two Figures. They corresponds to Figure 5A and Figure 5B in the paper. 

`Overall Survival Hazard Ratio.R` R file to run the OS_Grade.stan model. It will produce the Supplementary Figure 4

`Distant Metastasis Hazard Ratio.R` R file to run the DFS_Grade.stan model. It will produce the Supplementary Figure 3

`APM_Presence_Spt_Differences_Interpatient_Variability.R` First analysis part of the paper. It will produce the Figure 1,2,3,4.

`OS_Grade.stan` This file contains the model for Overall Survival according to APM's presence on Histotypes

`DFS_Grade.stan` This file contains the model for Disease Free Survival according to APM's presence on Histotypes

## Dataset

The datasets used in this project can be downloaded from Zenodo at the following link: ZENODO

## Requirements

The following softwares were used for developing and running this project:

* R version 4.1.2
* Stan version 2.21.0


## License

This project is licensed under the Creative Commons Attribution 4.0 International.
