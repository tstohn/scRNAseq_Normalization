#!/usr/bin/env Rscript

#parameters: 
# 1. <normalization method>
# 2. <name of dataset inside the dataset directory specified in settings.ini 
#     or 'ALL' to run for all datasets>

library(tidyverse)
library(here)
source(here("src/normalization", "functions.R"))

library(edgeR)

# NORNALIZATION FUNCTIONS

#outlier removal
is_value_outside_sd<-function(value_vector, sd)
{
  return(any(value_vector > 5*sd))
}
prefilter_dataset<-function(data)
{
  sample_ids_to_filter<-c()
  
  #get outlier: > 5 X sd
  outlier_sd<- data %>% 
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = sample_id, values_from = ab_count) %>%
    rowwise() %>% 
    mutate(sd = sd(across(where(is.numeric))))
  outlier_sd<-outlier_sd %>%
    purrr::map_df(is_value_outside_sd, outlier_sd$sd) %>%
    select(-sd, -ab_id) %>%
    pivot_longer(cols = everything()) %>%
    filter(value == TRUE)
  
  outlier_na<- data %>%
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = ab_id, values_from = ab_count) %>%
    rowwise() %>% 
    mutate(na=any(is.na(across(!starts_with("sample_id"))))) %>%
    select(sample_id, na) %>%
    filter(na==TRUE)
  
  sample_ids_to_filter<-append(sample_ids_to_filter, as.vector(outlier_sd$name))
  sample_ids_to_filter<-append(sample_ids_to_filter, as.vector(outlier_na$sample_id))
  
  data <- data[ !data$sample_id %in% sample_ids_to_filter, ]
  return(data)
}

run_tmm<-function(data)
{
  print("RUNNING TMM NORMALIZATION")
  countdata <- data %>%
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = ab_id, values_from = ab_count ) %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>% 
    t()
  
  normfactors <- edgeR::calcNormFactors(countdata,sumTrim = 0.05, logratioTrim = 0) %>% 
    enframe("sample_id", "normfactor")
  
  return(normfactors)
}

run_normalization<-function(dataset, method)
{
  datapath<-get_dataset_path()
  datapath<-paste(datapath, "/", dataset, sep="")
  output_folder<-"~/bin/NORMALIZED_DATASETS"
  data<-read_tsv(datapath)
  
  columns = get_dataset_specific_column_names(datapath)
  data <- data %>%
    rename(
      sample_id=columns[1],
      ab_id=columns[2],
      ab_count=columns[3],
      batch_id=columns[4],
      cluster_id=columns[5],
      ab_type=columns[6]
    ) %>%
    prefilter_dataset()
  
  print(paste("RUN NORMALIZATION[", method, "] FOR DATASET[", datapath, "]", sep = ""))
  if(method=="TMM")
  {
    x<-run_tmm(data)
  }
}

# MAIN FUNCTION
main<-function(method, dataset)
{
  datapath<-get_dataset_path()
  if(dataset=="ALL")
  {
    #no lookbehind in default R
    dataset_vector<-grep("(?<!ini)$", list.files(datapath, full.names = FALSE), perl=T, value=T)
    tmp<-sapply(dataset_vector, function(x) run_normalization(x, method))
  }
  else
  {
    run_normalization(dataset, method)
  }
}

# RSCRIPT
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=2)
{
  print("WRONG NUMBER OF ARGUMENTS:\'RScript NormalizationScript.R <method> <dataset>\'")
}
main(args[1], args[2])
