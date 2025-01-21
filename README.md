# Scripts, data and figures of Guery *et al.* 2024

This repository contains the pipeline, data and figures included in:

Guery, M.-A., Ceesay, S., Drammeh, S., Jaiteh, F. K., d’Alessandro, U., Bousema, T., Conway, D. J., & Claessens, A. (2024). Household clustering and seasonal genetic variation of Plasmodium falciparum at the community-level in The Gambia. *eLife*, *13*. [https://doi.org/10.7554/eLife.103047.1](https://doi.org/10.7554/eLife.103047.1).

# Pipeline

The details of this pipeline are described in the methods section of the corresponding article.

The script `barcode-genome_pipeline.R` is the main  script controlling each step of the pipeline. Beware that some steps require the output of previous ones.

The script `make-tables.R` uses the output of the pipeline to produce the tables attached with the article.

## Environment

The pipeline was run with R version 4.3.3 on Ubuntu 22.

```R
─ Session info ────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.3.3 (2024-02-29)
 os       Ubuntu 22.04.5 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language (EN)
 collate  C.UTF-8
 ctype    C.UTF-8
 tz       Europe/Paris
 date     2025-01-21
 rstudio  2024.12.0+467 Kousa Dogwood (desktop)
 pandoc   2.9.2.1 @ /usr/bin/pandoc

─ Packages ────────────────────────────────────────────────────────────────────────────────────────────
 package      * version    date (UTC) lib source
 ape          * 5.7-1      2023-03-13 [1] CRAN (R 4.3.2)
 binom        * 1.1-1.1    2022-05-02 [1] CRAN (R 4.3.3)
 BiocGenerics * 0.48.1     2023-11-01 [1] Bioconductor
 cli            3.6.3      2024-06-21 [1] CRAN (R 4.3.3)
 colorspace     2.1-1      2024-07-26 [1] CRAN (R 4.3.3)
 cowplot      * 1.1.3      2024-01-22 [1] CRAN (R 4.3.2)
 digest         0.6.36     2024-06-23 [1] CRAN (R 4.3.3)
 dplyr        * 1.1.4      2023-11-17 [1] CRAN (R 4.3.2)
 fansi          1.0.6      2023-12-08 [1] CRAN (R 4.3.2)
 farver         2.1.2      2024-05-13 [1] CRAN (R 4.3.3)
 generics       0.1.3      2022-07-05 [1] CRAN (R 4.3.2)
 ggforce      * 0.4.1      2022-10-04 [1] CRAN (R 4.3.2)
 ggnewscale   * 0.5.0      2024-07-19 [1] CRAN (R 4.3.3)
 ggpattern    * 1.0.1      2022-11-09 [1] CRAN (R 4.3.2)
 ggplot2      * 3.5.0      2024-02-23 [1] CRAN (R 4.3.3)
 glue           1.7.0      2024-01-09 [1] CRAN (R 4.3.3)
 gtable         0.3.5      2024-04-22 [1] CRAN (R 4.3.3)
 hms            1.1.3      2023-03-21 [1] CRAN (R 4.3.2)
 igraph       * 2.0.3      2024-03-13 [1] CRAN (R 4.3.3)
 IRanges      * 2.36.0     2023-10-24 [1] Bioconductor
 lattice        0.22-5     2023-10-24 [1] CRAN (R 4.3.2)
 lifecycle      1.0.4      2023-11-07 [1] CRAN (R 4.3.2)
 lubridate    * 1.9.3      2023-09-27 [1] CRAN (R 4.3.2)
 magrittr       2.0.3      2022-03-30 [1] CRAN (R 4.3.2)
 MASS           7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.2)
 munsell        0.5.1      2024-04-01 [1] CRAN (R 4.3.3)
 nlme           3.1-164    2023-11-27 [1] CRAN (R 4.3.2)
 pillar         1.9.0      2023-03-22 [1] CRAN (R 4.3.2)
 pkgconfig      2.0.3      2019-09-22 [1] CRAN (R 4.3.2)
 plyr           1.8.9      2023-10-02 [1] CRAN (R 4.3.2)
 polyclip       1.10-6     2023-09-27 [1] CRAN (R 4.3.2)
 purrr          1.0.2      2023-08-10 [1] CRAN (R 4.3.2)
 R6             2.5.1      2021-08-19 [1] CRAN (R 4.3.2)
 RColorBrewer * 1.1-3      2022-04-03 [1] CRAN (R 4.3.2)
 Rcpp           1.0.12     2024-01-09 [1] CRAN (R 4.3.2)
 readr        * 2.1.5      2024-01-10 [1] CRAN (R 4.3.2)
 reshape2     * 1.4.4      2020-04-09 [1] CRAN (R 4.3.2)
 rlang          1.1.4      2024-06-04 [1] CRAN (R 4.3.3)
 rstudioapi     0.15.0     2023-07-07 [1] CRAN (R 4.3.2)
 S4Vectors    * 0.40.2     2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
 scales       * 1.3.0      2023-11-28 [1] CRAN (R 4.3.2)
 sessioninfo    1.2.2      2021-12-06 [1] CRAN (R 4.3.2)
 stringi        1.8.4      2024-05-06 [1] CRAN (R 4.3.3)
 stringr      * 1.5.1      2023-11-14 [1] CRAN (R 4.3.2)
 tibble         3.2.1      2023-03-20 [1] CRAN (R 4.3.2)
 tidyr        * 1.3.1      2024-01-24 [1] CRAN (R 4.3.2)
 tidyselect     1.2.1      2024-03-11 [1] CRAN (R 4.3.3)
 timechange     0.3.0      2024-01-18 [1] CRAN (R 4.3.2)
 tweenr         2.0.2      2022-09-06 [1] CRAN (R 4.3.2)
 tzdb           0.4.0      2023-05-12 [1] CRAN (R 4.3.2)
 utf8           1.2.4      2023-10-22 [1] CRAN (R 4.3.2)
 vctrs          0.6.5      2023-12-01 [1] CRAN (R 4.3.2)
 withr          3.0.0      2024-01-16 [1] CRAN (R 4.3.2)
```

# Contact

Marc-Antoine Guery, postdoc at Antoine Claessens lab, LPHI, University of Montpellier.

Refer to corresponding author of the article for additional information.

