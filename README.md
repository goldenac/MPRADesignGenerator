# MPRADesignGenerator

## Installation ##

Use devtools to install MPRADesignGenerator. 
```
install.packages("devtools")
devtools::install_github("goldenac/MPRADesignGenerator")

library("MPRADesignTools")
```
MPRADesignGenerator relies on the **Biostrings** and **BSgenome.Hsapiens.UCSC.hg38** packages from Bioconductor. Make sure these packages are installed using the following commands.
```
BiocManager::install("Biostrings")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library("Biostrings")
library("BSgenome.Hsapiens.UCSC.hg38")
```

generate(variant_input_file, tag_input_file, scrambled_input_file)
