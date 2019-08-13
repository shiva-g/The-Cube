
library(optparse,quietly = T)
library(yaml,quietly = T)        
library(tidyverse, quietly = T)
library(dplyr, quietly = T)
library(ggplot2,quietly = T)
library(ggridges,quietly = T)
library(ggrepel, quietly = T)


capture <- commandArgs(trailingOnly = TRUE)

opt1 = list(make_option(c("--input"), type = "character", default = "input.yml", dest = "input"))


user_input <- function(name, argv) {
  return(name %in% names(argv))
}

argv <- parse_args(OptionParser(option_list = opt1))

if (user_input("input", argv)) {
  input = argv$input
  
  if(file.exists(input)){
    input.yaml <- read_yaml(input)
  }else{
    message('Input YAML not found.\n')
    break;
  }
} else {
  message('Cannot proceed without input yaml file. Please use "--input" flag .\n')
}

source("scripts/helper_file.R")
source("scripts/3d_arrays.R")

check <- file.exists(paste0(input.yaml$outputDir,'/results/ar_mat.rds'),paste0(input.yaml$outputDir,'/results/sr_mat.rds'),paste0(input.yaml$outputDir,'/results/ir_mat.rds'))

if(sum(check) != 3){
  message("All files (matrics) required to proceed, were not created.")
  break;
} else {
  source("scripts/hpo_associations.R")
  source("scripts/wilcoxon_test.R")
}


