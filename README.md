# Umi4C Package #

Welcome to umi4cPackage Bitbucket repository.

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


#### Installing misha package:
Since _misha_ package is not in CRAN, we will need to install it directly as follow:  
Download the package from: [here.](http://compgenomics.weizmann.ac.il/tanay/?page_id=617).   
Install the package:
```
!#r
install.packages("misha_3.2.7.tar.gz", repos=NULL) #Installs misha package from file
```


#### Importing UMI-4C package:
Download and install *umi4cPackage*: 
```
#!r
devtools::install_bitbucket("tanaylab/umi4cpackage")
library(umi4cPackage)
```

#### Using the package
Please refer to the package vignette for usage and complete step by step example.
'''
#!r
vignette('umi4c-usage', package=umi4cPackage)
'''