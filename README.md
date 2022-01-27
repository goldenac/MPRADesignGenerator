# MPRADesignGenerator

## Installation ##

Use devtools to install MPRADesignGenerator. 
```
install.packages("devtools")
devtools::install_github("goldenac/MPRADesignGenerator")

library("MPRADesignTools")
```
MPRADesignGenerator relies on the **Biostrings** and **BSgenome.Hsapiens.UCSC.hg38** packages from Bioconductor. Install these packages using the following commands.
```
BiocManager::install("Biostrings")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library("Biostrings")
library("BSgenome.Hsapiens.UCSC.hg38")
```

## Generating Your File ##
generate(variant_input_file, tag_input_file, scrambled_input_file)

Test MPRA

## About MPRADesignGenerator ##

- Purpose
- How it works

## Troubleshooting ##

- Problems loading packages
- File does not exist
