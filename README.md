# R2CPCT - Toolkit for the analysis of WGS and RNA-Seq data from the CPCT-02 and WIDE study.

![Author](https://img.shields.io/badge/Author-Job%20van%20Riet-orange.svg)
![license](https://img.shields.io/badge/license-%20GNU%20GPLv3-yellow.svg)

This package contains various *R* functions which are used to import, process and analyze WGS and (matching) RNA-Seq samples derived from the CPCT/Hartwig Medical Foundation collaboration in context of the CPCT-02, DRUP and WIDE studies.

## Installation

The latest development version (which requires R ≥v4.0) can easily be installed from this repository using devtools:
```R
# Install using devtools.
devtools::install_github(repo = "J0bbie/R2CPCT")

# Load library
library(R2CPCT)
```

### Additional annotation of genomic variants with VEP.

Prior to import, files need to be annotated with VEP using the workflow and scripts available in https://github.com/J0bbie/VariantAnnotation_VEP.
