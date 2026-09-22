#!/bin/R

args <- commandArgs(trailingOnly = TRUE)

COV  <- args[1]
PERC <- args[2]

#Set up workspace
setwd("/home/nhowe/epi-pop/Rmd")
getwd() #Display the current working directory

# Loading libraries
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("methylKit")

# Ran first time around, now is unnecessary
# install.packages("readxl")
# install.packages("tidyverse")
# install.packages("janitor")
library(methylKit)
library(readxl) ## read excel file 
library(plyr) ## needs to be loaded before dplyr 
library(tidyverse)
library(janitor) # for clean_names()

# report sessioninfo for package versions upon run
sessioninfo::session_info()

## Importing data
# loading metadata
metadata <- read_tsv("csv_outputs/POP_epigenetic_samples_metadata_AgeSex_20260729.tsv") %>% 
  dplyr::select(ABLG, SPECIMENID, AGE, SEX, LENGTH, WEIGHT) %>%
  mutate(ABLG = as.character(ABLG))

file_list <- list.files(path="/share2/nhowe/epi-pop/coverage",
                pattern = "*deduplicated.sorted.CpG_report.merged_CpG_evidence.cov.gz$", full.names=TRUE)

if(length(file_list)==0){
  print("Not calling in any files. Quitting")
  q(save = "no", status = 0)
}

print(paste(length(file_list),"will be included in the dataframe."))
unlist(file_list)

# File paths and sample IDs
# sample_ids <- sapply(file_list, function(x) {
#   str_split_i(basename(x), "_", 1)
# })
sample_ids <- unlist(file_list) %>% basename() %>% str_split_i("_", 1)

print("sample_ids:")
sample_ids[1:10]

# Pass numeric continuous ages into treatment
samples <- metadata %>%
  filter(ABLG %in% as.vector(sample_ids)) %>%
  mutate(ABLG = factor(ABLG, levels = as.vector(sample_ids))) %>%
  arrange(ABLG)
print("head(samples) arranged by ABLG:")
head(samples)

# Import coverage files
cov <- methRead(
  location = as.list(file_list),
  sample.id = as.list(sample_ids),
  assembly = "SebastesalutusPOP6.ragtag",
  # treatment = samples$AGE,
  treatment = rep(0, length(sample_ids)),
  pipeline = "bismarkCoverage",
  context = "CpG"
)

saveRDS(cov,"rdata/methylkit_methRead.rds")
print("methylkit methRead cov rds saved")

# hi.perc discards bases w/in top 99.9th percentile of coverage in each sample
filtered.cov <- filterByCoverage(cov,
                              lo.count=COV, # cov is an arg[1]
                              lo.perc=NULL,
                              hi.count=NULL,hi.perc=99.9)

saveRDS(filtered.cov, str_c("rdata/methylkit_methRead.cov",COV,"x.rds"))
print("filtered methylkit methRead cov rds saved")

getCoverageStats(filtered.cov, plot = FALSE)

norm.cov <- normalizeCoverage(filtered.cov)
saveRDS(norm.cov, str_c("rdata/methylkit_methRead.cov",COV,"x.norm.rds"))
print("filtered methylkit methRead cov rds saved")

getCoverageStats(norm.cov, plot = FALSE)

meth <- methylKit::unite(norm.cov, min.per.group = 1L, destrand=FALSE)
head(meth)

saveRDS(meth, str_c("rdata/methylkit_methRead.cov",COV,"x.norm.unite.rds"))
print("filtered, normalized, united methylkit methRead cov rds saved")

# sample clustering
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
ggsave("rdata/clusterplot1.jpeg")

diffMeth <- calculateDiffMeth(meth)

saveRDS(diffMeth, str_c("rdata/methylkit_methRead.cov",COV,"x.diffMeth.rds"))
print("filtered diff methylated methylkit rds saved")

diffMeth25p <- getMethylDiff(diffMeth,difference=25,qvalue=0.01)

saveRDS(diffMeth25p, str_c("rdata/methylkit_methRead.cov",COV,"x.diffMeth25p.rds"))
print("filtered diff methylated 25 p methylkit rds saved")

print("Diff meth per chr:")
diffMethPerChr(diffMeth25p,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
