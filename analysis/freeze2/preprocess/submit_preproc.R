#! /usr/bin/env Rscript
library(getopt)

spec = matrix(c(
      'cohort', 	'c', 1, "character"
    ), byrow=TRUE, ncol=4)
opt = getopt(spec)

stopifnot( ! opt$cohrot %in% c("RUSH", "HBCC", "MSSM", "FULL", "AGING"))

library(rmarkdown)
library(tidyverse)

path = "/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/h5ad_final/"

h5ad_files = c(
RUSH = "RUSH_2023-09-12_16_13.h5ad", 
HBCC = "HBCC_2023-09-12_16_28.h5ad",
MSSM = "MSSM_2023-09-12_17_04.h5ad",
FULL = "FULL_2023-09-12_18_32.h5ad",
AGING = "AGING_2023-09-12_20_35.h5ad")

h5ad_files = sapply(h5ad_files, function(x) paste0(path, x))

create_job = function(x){
  # system("git pull origin master");

  tmp.dir = paste0("/sc/arion/scratch/hoffmg01/", x, "_", round(runif(1)*1e7))
  dir.create(tmp.dir)
  setwd(tmp.dir)

  # create tmp Rmd for each dataset
  fldr = "/sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/preprocess/"
  file.orig = paste0(fldr, "preprocess.Rmd")
  file = paste0(fldr, "preprocess_", x,".Rmd") 
  file.copy(file.orig, file)

  # x %>% 
  #   walk(function(x) render(file,
  #           params = list(DATASET = h5ad_files[x]),
  #           output_file = paste0(fold, "/preprocess_", x, ".html")))
}

run = function(x){

  # copy to process_{DATASET}.Rmd and execute
  fldr = "/sc/arion/projects/CommonMind/hoffman/NPS-AD/work/nps_ad/analysis/freeze2/preprocess/"
  file1 = paste0(fldr, "preprocess.Rmd")
  file2 = paste0(fldr, "preprocess_", x,".Rmd")

  file.copy(file1, file2)

  render(file2,
          params = list(DATASET = h5ad_files[x]),
          output_file = paste0(fldr, "/preprocess_", x, ".html"))
}

names(h5ad_files) %>%
  walk(create_job)
  
run( opt$cohort)