# HK1FigureNotebook

[![DOI](https://zenodo.org/badge/472765261.svg)](https://zenodo.org/badge/latestdoi/472765261)

This julia project contains all source data and code to analyse data and regenerate genomics figures for:


Wakeling, M. N., Owens, N., Hopkinson, J. R., Johnson, M. B., Houghton, J., Dastamani, A., Flaxman, C. S., Wyatt, R. C., Hewat, T. I., Hopkins, J. J., Laver, T. W., van Heugten, R., Weedon, M. N., De Franco, E., Patel, K. A., Ellard, S., Morgan, N. G., Cheesman, E., Banerjee, I., Hattersley, A. T., … Flanagan, S. E. (2022). Non-coding variants disrupting a tissue-specific regulatory element in HK1 cause congenital hyperinsulinism. *Nature Genetics*, https://doi.org/10.1038/s41588-022-01204-x




Note repository contains source data for all gene expression and quantifications of ChIP-seq and ATAC-seq data at particular locu. Repository does not contain alignment files, but code is supplied to regenerate quantifications if suitable alignment files are present.

## Prerequistes
Julia >= 1.6, all julia packages and their versions are specified in the included `Project.toml` and `Manifest.toml`. Additionally, two python packages are used via https://github.com/JuliaPy/PyCall.jl:

  1. https://matplotlib.org/ as used by https://github.com/exeter-tfs/MotifScanner.jl to aid draw transcription factor motif sequence logos.
  2. https://github.com/open2c/cooler - uses the the cooler api to access HiC data in cooler format, this is only required should you wish to regenerate source data from cooler file.

If any errors are encountered generating figures ensure that https://github.com/JuliaPy/PyCall.jl and https://matplotlib.org/ is installed in the version of Python used to PyCall.


## Installation
```bash
git clone https://github.com/owensnick/HK1FigureNotebook.jl
cd HK1FigureNotebook.jl
julia
```
Within julia activiate the local 
```julia
] # to enter into Pkg mode
activate .
instantiate ## for first time installation
```
To regenerate figures either use jupyter notebook within `notebooks` directory or use script as follows:
```julia
 include("notebooks/hk1_genomics_figures.jl")
 ```
This will generate a `figures` folder and will generate all figure panels in `svg` and `png` format.