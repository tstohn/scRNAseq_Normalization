#!/usr/bin/env Rscript

#parameters: 
# 1. <normalization method>
# 2. <name of dataset inside the dataset directory specified in settings.ini 
#     or 'ALL' to run for all datasets>

library(tidyverse)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(here)
library(splitstackshape)
source(here("src/normalization", "functions.R"))

library(edgeR)

# NORNALIZATION FUNCTIONS

#outlier removal
is_value_outside_sd<-function(value_vector, sd, mean)
{
  return(any( (value_vector > mean+5*sd) | (value_vector < mean-5*sd) ))
}
prefilter_dataset<-function(data)
{
  sample_ids_to_filter<-c()
  ab_ids_to_filter<-c()
  
  #filter AB first, bcs according to AB outliers samples are removed
  #FILTER ANTI-BODIES with median over all AB(treatment, sample, ...) count < 40
  outlier_lowUmi_ab<- data %>%
    select(sample_id, ab_id, ab_count) %>%
    group_by(ab_id) %>%
    summarise(median = median(ab_count)) %>%
    filter(median < 40) %>%
    ungroup()
  
  ab_ids_to_filter<-append(sample_ids_to_filter, as.vector(outlier_lowUmi_ab$ab_id))
  data <- data[ !data$ab_id %in% ab_ids_to_filter, ]
  
  #FILTER SAMPLES
  outlier_sd<- data %>% 
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = sample_id, values_from = ab_count) %>%
    rowwise() %>% 
    mutate(sd = sd(across(starts_with("plate")), na.rm = TRUE),
           mean = rowMeans(across(starts_with("plate")), na.rm = TRUE))
  outlier_sd<-outlier_sd %>%
    purrr::map_df(is_value_outside_sd, outlier_sd$sd, outlier_sd$mean) %>%
    select(-sd, -ab_id) %>%
    pivot_longer(cols = everything()) %>%
    filter(value == TRUE)
  
  #simply based on first dataset, many low count samples were seen
  #maybe integrate a statistic approach
  outlier_lowUmi_sample<- data %>%
    select(sample_id, ab_id, ab_count) %>%
    group_by(sample_id) %>%
    summarise(sum = sum(ab_count)) %>%
    filter(sum < 25000) %>%
    ungroup()
  
  outlier_na<- data %>%
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = ab_id, values_from = ab_count) %>%
    rowwise() %>% 
    mutate(na=any(is.na(across(!starts_with("sample_id"))))) %>%
    select(sample_id, na) %>%
    filter(na==TRUE)

  sample_ids_to_filter<-append(sample_ids_to_filter, as.vector(outlier_sd$name))
  sample_ids_to_filter<-append(sample_ids_to_filter, as.vector(outlier_na$sample_id))
  sample_ids_to_filter<-append(sample_ids_to_filter, as.vector(outlier_lowUmi_sample$sample_id))
  
  print(paste0("Removing outliers: ", length(unique(sample_ids_to_filter))))
  data <- data[ !data$sample_id %in% sample_ids_to_filter,]
  
  return(data)
}

remove_batch_effect_manually <- function(data){
   data %>%
     group_by(ab_id) %>%
     dplyr::mutate(overall_median_ab_count = median(ab_count_normalized)) %>%
     group_by(ab_id, batch_id) %>%
     dplyr::mutate(plate_median_ab_count = median(ab_count_normalized)) %>%
     dplyr::mutate(ab_count_normalized = (ab_count_normalized/plate_median_ab_count) *
                     overall_median_ab_count) %>%
     select(-overall_median_ab_count, -plate_median_ab_count) %>%
     ungroup()
}

remove_batch_effect<-function(data, log_transform = FALSE)
{
  if(log_transform)
  {
    data <- data %>% mutate(ab_count_normalized = log(ab_count_normalized))
  }
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

sample_x_rows_of_sample<-function(sample, data, x)
{
  data[data$sample_id == sample, ] %>% 
    expandRows("ab_count") %>%
    sample_n(x, replace=FALSE) %>%
    group_by_all() %>%
    dplyr::summarise(ab_count_normalized = n()) %>%
    ungroup()
}

count_different_ab_for_sample<-function(sample, data)
{
  data[data$sample_id == sample, ] %>% 
    nrow()
}

run_subsampling<-function(data)
{
  #get number of antibodies per sample before subsampling 
  #(subsampling might introduce empty AB counts that must be removed)
  ab_number_before_subsampling <- purrr::map(unique(data$sample_id), count_different_ab_for_sample, data)
  assertthat::assert_that(length(unique(ab_number_before_subsampling)) == 1)
  ab_number_before_subsampling <- unique(ab_number_before_subsampling) %>% as.numeric()
  
  #minimum samplesize to use as subsampling size
  sample_size <- data %>%
    group_by(sample_id) %>%
    summarise(sum = sum(ab_count)) %>%
    summarise(min = min(sum)) %>%
    ungroup() %>%
    as.numeric()
    
  #subsample
  normalized_data <-
    purrr::map_dfr(unique(data$sample_id), sample_x_rows_of_sample, data, sample_size)
    
  #remove new zero AB counts
  sample_ids_to_filter<-c()
  outlier <- normalized_data %>% group_by(sample_id) %>%
    summarise(n = n()) %>%
    filter(n != ab_number_before_subsampling) %>%
    ungroup()
  sample_ids_to_filter<-append(sample_ids_to_filter, as.vector(outlier$sample_id))
  if(length(sample_ids_to_filter) > 0)
  {
    print(paste0("WARNING: Removing ", length(sample_ids_to_filter), " samples due to zero AB counts after subsampling"))
  }
  #quality check, arbitrarily set threshold to not remove more than 10% of samples due to subsampling
  assertthat::assert_that(length(sample_ids_to_filter) < 0.1*length(unique(data$sample_id)))
  normalized_data <- normalized_data[ !normalized_data$sample_id %in% sample_ids_to_filter,]
  
  return(normalized_data)
}

run_normalization<-function(dataset, method)
{
  datapath<-get_dataset_path()
  datapath<-paste(datapath, "/", dataset, sep="")
  data<-read_tsv(datapath)
  dir.create(paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset)), showWarnings = TRUE)
  
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
      remove_batch_effect(log_transform = TRUE)
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/TMM.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="SUBSAMPLING")
  {
    normalized_data<- data %>% 
      run_subsampling() %>%
      remove_batch_effect(log_transform = TRUE)
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/SUBSAMPLED.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="SCTRANSFORM")
  {

  }
  else if(method=="CLR_COMPOSITIONS")
  {

  }
  else
  {
    print(paste0("WARNING: No normalization metjod called: ", method))
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
