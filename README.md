# MPRADesignGenerator


## Installation ##

Use devtools to install MPRADesignGenerator. If you do not have devtools installed, use the following:
```
install.packages("devtools")
```
Then install and load the package:
```
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

- Required files and format
- Using the generate function: generate(variant_input_file, tag_input_file, scrambled_input_file)
- Test: Provide test files, command to run test files, and explain what should be output


## About MPRADesignGenerator ##

- Purpose
- What is output
- How it works
    * How are alt sequences created for each variant type
    * How are restriction enzyme sites handled


## Troubleshooting ##

- Problems loading packages
- File does not exist
