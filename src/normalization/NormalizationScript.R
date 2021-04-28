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
  ab_ids_to_filter<-c()
  
  #filter AB first, bcs according to AB outliers samples are removed
  #FILTER ANTI-BODIES with mean overall (treatment, sample, ...) count < 40
  outlier_lowUmi<- data %>%
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = sample_id, values_from = ab_count) %>%
    column_to_rownames("ab_id") %>%
    mutate(mean=rowMeans(., na.rm = TRUE)) %>%
    filter(mean < 40) %>%
    rownames_to_column(var="ab_id") %>%
    select(ab_id, mean)
  
  ab_ids_to_filter<-append(sample_ids_to_filter, as.vector(outlier_lowUmi$ab_id))
  data <- data[ !data$ab_id %in% ab_ids_to_filter, ]
  
  #FILTER SAMPLES
  outlier_sd<- data %>% 
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = sample_id, values_from = ab_count) %>%
    rowwise() %>% 
    mutate(sd = sd(across(starts_with("plate")), na.rm = TRUE))
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

  data <- data[ !data$sample_id %in% sample_ids_to_filter,]
  
  return(data)
}

remove_batch_effect<-function(data)
{
  matrix_for_batch_correction<-select(data, sample_id, ab_id, ab_count_normalized) %>%
    pivot_wider(names_from = sample_id, values_from = ab_count_normalized) %>%
    column_to_rownames("ab_id") %>%
    as.matrix() 
  sample_order_in_correction_matrix <- colnames(matrix_for_batch_correction)
  batch_matrix<- data %>% 
    select(sample_id, batch_id) %>%
    unique()
  batch_vector <- batch_matrix[match(sample_order_in_correction_matrix, batch_matrix$sample_id),] %>%
    pull(batch_id)
  
  batch_corrected_data <- limma::removeBatchEffect(matrix_for_batch_correction, batch = batch_vector) %>% 
    as.data.frame() %>% 
    rownames_to_column("ab_id") %>% 
    pivot_longer(-ab_id, names_to = "sample_id", values_to = "ab_count_normalized")
  
  data_norm<-data %>%
    select(-ab_count_normalized) %>% 
    left_join(batch_corrected_data)
  
  return(data_norm)
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
  
  libsize <- data %>% 
    group_by(sample_id) %>% 
    mutate(lib_size = sum(ab_count)) %>% 
    select(sample_id, lib_size) %>% 
    ungroup() %>% 
    distinct()
  normalized_data <- data %>% 
    left_join(normfactors, by = "sample_id") %>% 
    left_join(libsize, by = "sample_id") %>%  
    mutate(ab_count_normalized = ab_count/(lib_size*normfactor))
  
  return(normalized_data)
}

run_normalization<-function(dataset, method)
{
  datapath<-get_dataset_path()
  datapath<-paste(datapath, "/", dataset, sep="")
  output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), dataset)
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
    normalized_data<- data %>% 
      run_tmm() %>%
      remove_batch_effect()
    write_tsv(normalized_data, file=output_table)
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
