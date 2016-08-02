## ---- eval=FALSE---------------------------------------------------------
#  install.packages("misha_3.2.7.tar.gz", repos=NULL) #Installs misha package from file

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_bitbucket("tanaylab/umi4cpackage")
#  library(umi4cPackage)

## ---- eval=FALSE---------------------------------------------------------
#  library("umi4cPackage")
#  # Use build 'hg19' for DpnII enzyme
#  p4cBuildRequirements(output_dir = "my_u4c_proj", genome = "hg19", reseq = "GATC")

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  p4cDumpConfFiles(conf_dir = "my_u4c_proj/conf")

## ---- eval=FALSE---------------------------------------------------------
#  p4cLoadConfFiles(conf_dir = "my_u4c_proj/conf")

## ---- eval=FALSE---------------------------------------------------------
#  p4cCreate4CseqTrack(sample_ids = 1)

## ---- eval=FALSE---------------------------------------------------------
#  CMK_fc <- p4cNewProfile("umi4C_example_CMK_ANK1_TSS", scope_5=200000, scope_3=200000)

## ---- eval=FALSE---------------------------------------------------------
#  plot(CMK_fc)

## ---- eval=FALSE---------------------------------------------------------
#  plot(CMK_fc, png_fn="figs/CMK_ANK1_TSS.png")

## ---- eval=FALSE---------------------------------------------------------
#  p4cCreate4CseqTrack(sample_ids = c(2,3))

## ---- message=FALSE, eval=FALSE------------------------------------------
#  # Generate a second p4cProfile on the same scope
#  DND41_fc <- p4cNewProfile("umi4C_example_DND41_ANK1_TSS", scope_5=200000, scope_3=200000)
#  
#  # plot a comparative plot
#  plot(CMK_fc, DND41_fc)

## ---- message=FALSE, eval=FALSE------------------------------------------
#  fold_change <- p4cIntervalsMean(CMK_fc, DND41_fc, start=41665000, end=41675000)
#  knitr::kable(fold_change, align = 'l', digits=2)

