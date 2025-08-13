# WIM_HB

Welcome to the GitHub repository for 
**When your heart isn't in it anymore: Cardiac correlates of task disengagement**, a combined EEG-ECG-pupillometry study of brain-body interaction during spontaneous cognition.

This repository contains the scripts and functions required to reproduce the analyses reported in the associated [preprint](https://www.biorxiv.org/content/10.1101/2024.06.21.599851v2).
You can access the data files on which these scripts operate from the accompanying OSF repositories: [Melbourne data](https://osf.io/ey3ca/); [Paris data](https://osf.io/v9xsw/).

## Directions
In order to reproduce the analyses reported in the manuscript, first download/clone this repository into a suitable directory and navigate to the `localdef_WIM_HB.m` script.
You will need to modify this script in order to define local paths to the data (which you will need to download from the OSF, as well as the requisite *MATLAB* toolboxes.

Scripts for EEG, ECG, and eye-tracker preprocessing are located in the `preproc` subdirectory.
The pipeline can be run in its entirety via the high-level script `WIM_HB_preproc_pipeline.m`.
Preprocessed files can then be called by relevant scripts in the `analysis` subdirectory to emulate the signal processing and mixed-effects modelling reported in the manuscript.

All `MATLAB` analysis scripts can be invoked by running `WIM_HB_analysis_pipeline.m`.
Note, the mixed-effects analysis can be run independently of the preprocessing pipeline by downloading the relevant intermediary datafiles from the `analysis` component of the [Paris repository](https://osf.io/pc74r/).

The results of these analysis procesdures can be visualised by running the `plotHEP.m` and `plotGCMI.m` scripts (`analysis\plots` subdirectory).

## Citation
If the materials archived here are useful for your own research, please cite this repository including the appropriate [release version](#current-release) information (year and doi; see below for details):

> Corcoran, A.W., Le Coz, A., Hohwy, J., & Andrillon, T. {*year*}. WIM_HB {*version*} [Software]. Retrieved from https://github.com/corcorana/WIM_HB. {*doi*}

Please also cite the accompanying manuscript.

## License
This software is freely available for redistribution and/or modification under the terms of the GNU General Public Licence.
It is distributed WITHOUT WARRANTY; without even the implied warranty of merchantability or fitness for a particular purpose. 
See the [GNU General Public License](https://github.com/corcorana/SWS_NVS_code/blob/main/LICENSE) for more details.

## System information
The preprocessing and analysis pipeline archived here was built and tested on a 64-bit system running Microsoft Windows 10 Enterprise Version 10.0 (Build 19045).
*MATLAB* and *R* software version details, and attached toolboxes/packages etc., are listed below.


**MATLAB version 9.12.0.1884302 (R2022a)**

**Java Version**: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

**R version 4.2.1 (2022-06-23 ucrt)**

**Platform:** x86_64-w64-mingw32/x64 (64-bit)

**attached base packages:**

[1] stats     graphics  grDevices utils     datasets  methods   base     

**other attached packages:**

 [1] cowplot_1.1.1     ggh4x_0.2.8       ggeffects_1.1.2   MASS_7.3-58.1    
 [5] emmeans_1.7.5     itsadug_2.4.1     plotfunctions_1.4 mgcViz_0.1.9     
 [9] qgam_1.3.4        mgcv_1.8-40       nlme_3.1-158      forcats_0.5.1    
[13] stringr_1.5.0     dplyr_1.1.2       purrr_1.0.1       readr_2.1.2      
[17] tidyr_1.3.0       tibble_3.2.1      ggplot2_3.5.1     tidyverse_1.3.2  
[21] pacman_0.5.1     

**loaded via a namespace (and not attached):**

 [1] matrixStats_0.62.0  fs_1.5.2            lubridate_1.8.0     doParallel_1.0.17  
 [5] RColorBrewer_1.1-3  httr_1.4.3          tools_4.2.1         backports_1.4.1    
 [9] utf8_1.2.2          R6_2.5.1            KernSmooth_2.23-20  DBI_1.1.3          
[13] colorspace_2.0-3    withr_2.5.0         tidyselect_1.2.0    gridExtra_2.3      
[17] GGally_2.1.2        compiler_4.2.1      cli_3.5.0           rvest_1.0.2        
[21] xml2_1.3.3          scales_1.3.0        mvtnorm_1.1-3       digest_0.6.29      
[25] minqa_1.2.4         pkgconfig_2.0.3     htmltools_0.5.3     lme4_1.1-30        
[29] dbplyr_2.2.1        fastmap_1.1.0       rlang_1.1.1         readxl_1.4.0       
[33] rstudioapi_0.13     shiny_1.7.2         generics_0.1.3      jsonlite_1.8.0     
[37] googlesheets4_1.0.0 magrittr_2.0.3      Matrix_1.4-1        Rcpp_1.0.9         
[41] munsell_0.5.0       fansi_1.0.3         viridis_0.6.2       lifecycle_1.0.3    
[45] stringi_1.7.8       plyr_1.8.7          grid_4.2.1          parallel_4.2.1     
[49] promises_1.2.0.1    crayon_1.5.2        miniUI_0.1.1.1      lattice_0.20-45    
[53] haven_2.5.0         splines_4.2.1       hms_1.1.1           pillar_1.9.0       
[57] boot_1.3-28         estimability_1.4    codetools_0.2-18    reprex_2.0.1       
[61] glue_1.6.2          modelr_0.1.8        vctrs_0.6.2         nloptr_2.0.3       
[65] tzdb_0.3.0          httpuv_1.6.5        foreach_1.5.2       cellranger_1.1.0   
[69] gtable_0.3.1        reshape_0.8.9       assertthat_0.2.1    mime_0.12          
[73] xtable_1.8-4        broom_1.0.0         coda_0.19-4         later_1.3.0        
[77] googledrive_2.0.0   viridisLite_0.4.1   gargle_1.2.0        iterators_1.0.14   
[81] gamm4_0.2-6         ellipsis_0.3.2     