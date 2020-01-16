# Umi4C Package #

Welcome to umi4cPackage Github repository.

### What is umi4cPackage? ###
umi4cPackage is a processing and analysis pipeline for UMI-4C experiment. 

## Installation
### Requirements:
The pipeline is designed to run on a standard linux machine. The basic requirements are: 

- _Perl_
- _Bowtie2_: <http://bowtie-bio.sourceforge.net/bowtie2>.
- Reference index for _Bowtie2_.
- R packages:
    * _devtools_.
    * _misha_.
    * _zoo_.


#### Installation
```
#!r
devtools::install_github("tanaylab/umi4cpackage", vignette = TRUE)
library(umi4cPackage)
```

#### Using the package
Please refer to the package vignette for usage and a complete step by step example.
```
#!r
vignette('umi4c-usage', package='umi4cPackage') 
```
