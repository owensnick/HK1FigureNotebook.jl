# HK1FigureNotebook

This julia project contains all source data and code to analyse data and regenerate figures for:



A novel disease mechanism leading to the expression of a disallowed gene in the pancreatic beta-cell identified by non-coding, regulatory mutations controlling HK1
Matthew N. Wakeling, Nick D. L. Owens, Jessica R. Hopkinson, Matthew B. Johnson, Jayne A.L. Houghton, Antonia Dastamani, Christine S. Flaxman, Rebecca C. Wyatt, Thomas I. Hewat, Jasmin J. Hopkins, Thomas W. Laver, Rachel Van Heugten, Michael N. Weedon, Elisa De Franco, Kashyap A. Patel, Sian Ellard, Noel G. Morgan, Edmund Cheesman, Indraneel Banerjee, Andrew T. Hattersley, Mark J. Dunne, International Congenital Hyperinsulinism Consortium, Sarah J. Richardson, Sarah E. Flanagan

medRxiv 2021.12.03.21267240; doi: https://doi.org/10.1101/2021.12.03.21267240


Note repository contains source data for all gene expression and quantifications of ChIP-seq and ATAC-seq data at particular locu. Repository does not contain alignment files, but code is supplied to regenerate quantifications if suitable alignment files are present.

## Prerequistes
All julia packages and their versions are specified in the included `Project.toml` and `Manifest.toml`. Additionally, two python packages are used via https://github.com/JuliaPy/PyCall.jl:

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
