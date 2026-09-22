#!/bin/R

args <- commandArgs(trailingOnly = TRUE)

COV  <- args[1]
PERC <- args[2]

#Set up workspace
setwd("/home/nhowe/epi-pop/Rmd")
getwd() #Display the current working directory

# Ran first time around, now is unnecessary
# install.packages("readxl")
# install.packages("tidyverse")
# install.packages("janitor")

#Loading libraries
library(readxl) ## read excel file 
library(plyr) ## needs to be loaded before dplyr 
library(tidyverse)
library(janitor) # for clean_names()

# report sessioninfo for package versions upon run
sessioninfo::session_info()

## Importing data
# filenames <- list.files(path="/share2/nhowe/epi-pop/coverage/filtered/20X",
#                 pattern = ".*R1.*90perc.*\\.bed$", full.names=TRUE)
filenames <- list.files(path=str_c("/share2/nhowe/epi-pop/coverage/filtered/",COV,"X"),
                pattern = str_c(".*R1.*",PERC,"perc.*\\.bed$"), full.names=TRUE)

if(length(filenames)==0){
  print("Not calling in any files. Quitting")
  q(save = "no", status = 0)
}

print(paste(length(filenames),"will be included in the dataframe."))

# combine files into a dataframe
df <- filenames %>% 
  set_names(.) %>% map_dfr(read.table, .id = "sample", header=F)

# clean up sample names
df$sample <- basename(df$sample) %>% gsub(pattern = "_.*", replacement = "\\1", .)

# changing column names
colnames(df) <- c("ABLG", "scaffold", "start", "stop", "percent.meth", "meth", "unmeth")
head(df)

df <- df %>%
  mutate(locus = str_c(scaffold,"_",start))

invariant_loci <- df %>%
  group_by(locus) %>%
  dplyr::summarise(
    invar = (max(percent.meth, na.rm = TRUE) == min(percent.meth, na.rm = TRUE)),
    barelyVar90 = ifelse(min(percent.meth, na.rm = TRUE) > 90 | max(percent.meth, na.rm = TRUE) < 10, TRUE, FALSE)
  ) %>%
  filter(invar | barelyVar90) %>%
  pull(locus)

print("# of invariant or barely variant loci:")
nrow(invariant_loci)

df <- df %>%
  filter(grepl("^NC",scaffold)) %>%
  filter(!(locus %in% invariant_loci))

# loading metadata
#metadata <- read_tsv("csv_outputs/POP_epigenetic_samples_metadata_AgeSex_20260729.tsv") %>% 
#  dplyr::select(ABLG, SPECIMENID, AGE, SEX, LENGTH, WEIGHT) %>%
#  mutate(ABLG = as.character(ABLG))

#df <- left_join(df, metadata, by = "ABLG")

#print("After left join")
#head(df)

## adds one column for 21994932 x 12
#df <- df %>% unite(Loc, c("scaffold", "start"), sep=" ", remove=F)

#print("after unite")
print("head df:")
head(df)

print("number of loci:")
nrow(unique(df$locus))

# Save Rdata
saveRDS(df, file = str_c("rdata/bedfiles-",COV,"x-",PERC,"perc.variant.rds"))
print("rds saved")
