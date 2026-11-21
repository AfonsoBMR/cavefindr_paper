(in submission process) Rodrigues, A. & Fazenda, P. (2026) Detecting cave entrances from open airborne LiDAR data: the cavefindr R package and a case study in mainland Portugal

Script "validation.R" is the code necessary to reproduce the analysis and figures in the paper. It was run on cavefindr version 1.0.0, which can be installed using 

```r 
remotes::install_github("AfonsoBMR/cavefindr@f8029c74ea403c24cc372d2ac3df04244b2d2518")
```

Folder "data" contains all supporting data except the LiDAR tiles, which are too heavy to include here and can be sourced from DGT using the following link:  https://cdd.dgterritorio.gov.pt/dgt-fe/catalogos/colecoes?filter=LiDAR+-+Modelos+Digitais+do+Relevo

Table "verified_coordinates.csv" are the results from field verification.

**R session parameters used to run the code:**

R version 4.4.3 (2025-02-28 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8    LC_MONETARY=English_United Kingdom.utf8
[4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.utf8    

time zone: Europe/Lisbon
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] future_1.82.0

loaded via a namespace (and not attached):
 [1] gtable_0.3.6         ggplot2_4.0.4        raster_3.6-32        devtools_2.5.2       lattice_0.23-1       vctrs_0.7.3         
 [7] tools_4.4.3          generics_0.1.4       parallel_4.4.3       tibble_3.3.1         proxy_0.4-29         pkgconfig_2.0.3     
[13] KernSmooth_2.23-27   data.table_1.18.6.1  checkmate_2.3.4      RColorBrewer_1.1-3   S7_0.2.2             lifecycle_1.0.5     
[19] compiler_4.4.3       farver_2.1.2         progress_1.2.3       terra_1.9-50         codetools_0.2-20     stars_0.7-3         
[25] usethis_3.2.2        class_7.3-24         lazyeval_0.2.3       pillar_1.11.1        crayon_1.5.3         ellipsis_0.3.3      
[31] classInt_0.4-11      cachem_1.1.0         lidR_4.3.3-9000      sessioninfo_1.2.4    abind_1.4-8          parallelly_1.50.0   
[37] digest_0.6.39        tidyselect_1.2.1     listenv_1.1.0        sf_1.1-3             dplyr_1.2.1          purrr_1.2.2         
[43] cowplot_1.2.0        rprojroot_2.1.1      fastmap_1.2.0        grid_4.4.3           here_1.0.2           cli_3.6.6           
[49] magrittr_2.0.5       randomForest_4.7-1.2 pkgbuild_1.4.8       e1071_1.7-17         prettyunits_1.2.0    scales_1.4.0        
[55] backports_1.5.1      sp_2.2-4             globals_0.20.0       rlas_1.9.5           otel_0.3.0           hms_1.1.4           
[61] memoise_2.0.1        cavefindr_1.0.0      rlang_1.3.0          Rcpp_1.1.2           glue_1.8.1           DBI_1.3.0           
[67] pkgload_1.5.3        rstudioapi_0.19.0    R6_2.6.1             fs_2.1.0             units_1.0-1         