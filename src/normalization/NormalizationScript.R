#!/usr/bin/env Rscript

#parameters: 
# 1. <normalization method>
# 2. <name of dataset inside the dataset directory specified in settings.ini 
#     or 'ALL' to run for all datasets>

suppressPackageStartupMessages(deployrUtils::deployrPackage("tidyverse"))
options(tidyverse.quiet = TRUE)
deployrUtils::deployrPackage("dplyr", warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
deployrUtils::deployrPackage("here")
deployrUtils::deployrPackage("splitstackshape")
source(here("src/normalization", "functions.R"))
library("edgeR")
suppressPackageStartupMessages(deployrUtils::deployrPackage("compositions"))
library(Seurat)

# NORNALIZATION FUNCTIONS

#outlier removal
is_value_outside_sd<-function(value_vector, sd, mean)
{
  return(any( (value_vector > mean+5*sd) | (value_vector < mean-5*sd) ))
}

#TODO: works only for one dataset now:
# - remove samples of low total AB counts according to distribution of samples
# - remove AB wit low median count/ zero count more sophistically
prefilter_dataset<-function(data, thresholds)
{
  sample_ids_to_filter<-c()
  ab_ids_to_filter<-c()
  
  #filter AB first, bcs according to AB outliers samples are removed
  #FILTER ANTI-BODIES with median over all AB(treatment, sample, ...) count < 40
  outlier_lowUmi_ab<- data %>%
    select(sample_id, ab_id, ab_count) %>%
    group_by(ab_id) %>%
    summarise(median = median(ab_count)) %>%
    filter(median < as.numeric(thresholds[2])) %>%
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
    filter(sum < as.numeric(thresholds[1])) %>%
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
  
  print(paste0("Removing outliers: ", length(unique(sample_ids_to_filter)), " of in total: ", length(unique(data$sample_id)), " samples"))
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
  #TMM takes data in a matrix of features X samples
  countdata <- data %>%
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = ab_id, values_from = ab_count ) %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>% 
    t()
  
  normfactors <- edgeR::calcNormFactors(countdata,sumTrim = 0.05, logratioTrim = 0) %>% 
    enframe("sample_id", "normfactor")
  #on certain platform the output of normfactors is a tibble
  #with the rownames of countdata (Mac, Ubuntu)
  #however for e.g. OpenSuse its simply the rownumber
  normfactors$sample_id <- colnames(countdata)
  normfactors <- normfactors %>%  
    mutate(sample_id = as.character(sample_id))
  
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

run_clr_seurat<-function(data)
{
  print("RUNNING CLR SEURAT NORMALIZATION")
  #Seurat first creates a SeuratObject
  #it takes data in a matrix of features X samples
  countdata <- data %>%
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = ab_id, values_from = ab_count ) %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>% 
    t()
  #create seurat object
  seurat_data <- CreateSeuratObject(counts = countdata)
  #normalize
  normalized_data <- NormalizeData(seurat_data, normalization.method = "CLR")
  #regain a dataframe of features X samples
  normalized_data_dgCMatrix <- normalized_data@assays$RNA@data
  normalized_data_frame <- as.data.frame(as.matrix(normalized_data_dgCMatrix))
  
  final_normalized_data <- normalized_data_frame %>%
    tibble::rownames_to_column(var = "ab_id") %>%
    pivot_longer(-ab_id, names_to = "sample_id", values_to = "ab_count_normalized")
    
  #map the results to our origional data frame
  combined_data <- left_join(data, final_normalized_data, by = c("sample_id", "ab_id"))
  
  return(combined_data)
}

run_sctransform<-function(data, batchEffect) 
{
  if(batchEffect)
  {
    print("RECENTLY NO DEDICATED BATCH EFF REMOVAL FOR SCTRANSFORM")
  }
  print("RUNNING SCTRANSFORM NORMALIZATION FROM SEURAT")
  #Seurat first creates a SeuratObject
  #it takes data in a matrix of features X samples
  countdata <- data %>%
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = ab_id, values_from = ab_count ) %>%
    column_to_rownames("sample_id") %>%
    as.matrix() %>% 
    t()
  #create seurat object
  seurat_data <- CreateSeuratObject(counts = countdata)
  #normalize
  normalized_data <- SCTransform(seurat_data)
  #regain a dataframe of features X samples
  #SCTransformed values are in the new assay SCT
  normalized_data_dgCMatrix <- normalized_data@assays$SCT@data
  normalized_data_frame <- as.data.frame(as.matrix(normalized_data_dgCMatrix))
  
  final_normalized_data <- normalized_data_frame %>%
    tibble::rownames_to_column(var = "ab_id") %>%
    pivot_longer(-ab_id, names_to = "sample_id", values_to = "ab_count_normalized")
  
  #map the results to our origional data frame
  combined_data <- left_join(data, final_normalized_data, by = c("sample_id", "ab_id"))
  
  return(combined_data)

run_leave_one_out_tmm<-function(data)
{
  print("RUNNING LEAVE ONE OUT TMM NORMALIZATION")
  
  CbindedTMMOutput <- data
  colVector = c()
  i <- 0
  for (abId in unique(data$ab_id))
  {
    tmpData <- data[data$ab_id != abId,]
    tmpTMM <- run_tmm(tmpData)
    newColName <- paste0("LEAVEONEOUTTMM", as.character(i))
    colVector = append(colVector, newColName)
    tmpTMM <- tmpTMM %>%
      select(sample_id, ab_id, ab_count_normalized) %>%
      rename( !!newColName := ab_count_normalized)
    CbindedTMMOutput <- merge(x = CbindedTMMOutput, y = tmpTMM, by = c('sample_id','ab_id'), all.x = TRUE)
    i = i+1
  }

  normalized_data <- CbindedTMMOutput
  normalized_data$ab_count_normalized <- rowSums(normalized_data[,grepl( "LEAVEONEOUTTMM" , names( normalized_data ) )], na.rm=TRUE)
  normalized_data <- normalized_data %>%
    select(-all_of(colVector)) %>%
    mutate(ab_count_normalized = ab_count_normalized / (length(colVector)-1))
    
  return(normalized_data)
}

run_clr<-function(data)
{
  print("RUNNING CLR (COMPOSITION) NORMALIZATION")
  #CLR takes data in a matrix of samples X features
  countdata <- data %>%
    select(sample_id, ab_id, ab_count) %>%
    pivot_wider(names_from = ab_id, values_from = ab_count ) %>%
    column_to_rownames("sample_id") %>%
    as.matrix()
  
  clr_normalized_counts<-as.data.frame(compositions::clr(countdata)) %>%
    rownames_to_column(var = "sample_id") %>% 
    pivot_longer(-sample_id, names_to = "ab_id", values_to = "ab_count_normalized")
  
  normalized_data <- data %>%
    left_join(clr_normalized_counts, by=c("sample_id", "ab_id"))
  
  return(normalized_data)
}

sample_x_rows_of_sample<-function(sample, data, x, fill_with_zeros)
{
  subsampled_data <- data[data$sample_id == sample, ] %>% 
    expandRows("ab_count") %>%
    sample_n(x, replace=FALSE) %>%
    group_by_all() %>%
    dplyr::summarise(ab_count_normalized = n()) %>%
    ungroup()
  
  zero_rows = data[data$sample_id == sample, ] %>%
    select(-ab_count) %>%
    mutate(ab_count_normalized = 0) %>%
    unique()
  
  if(fill_with_zeros){
    subsampled_plus_zero <- subsampled_data %>%
      full_join(zero_rows) %>%
      group_by(across(c(-ab_count_normalized))) %>%
      dplyr::summarise(ab_count_normalized = sum(ab_count_normalized)) %>%
      ungroup()
  }
  else
  {
    return(subsampled_data)
  }
    
  return(subsampled_plus_zero)
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
  fill_nan_with_zero = TRUE
  normalized_data <-
    purrr::map_dfr(unique(data$sample_id), sample_x_rows_of_sample, data, sample_size, fill_nan_with_zero)

  #remove new zero AB counts
  if(!fill_nan_with_zero)
  {
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
  }
  
  return(normalized_data)
}

run_leave_one_out_subsampling<-function(data) 
{
  print("RUNNING LEAVE ONE OUT SUBSAMPLING NORMALIZATION")
  
  CbindedSubsamplingOutput <- data
  colVector = c()
  i <- 0
  for (abId in unique(data$ab_id))
  {
    tmpData <- data[data$ab_id != abId,]
    tmpSubs <- run_subsampling(tmpData)
    newColName <- paste0("LEAVEONEOUTSUBS", as.character(i))
    colVector = append(colVector, newColName)
    tmpSubs <- tmpSubs %>%
      select(sample_id, ab_id, ab_count_normalized) %>%
      rename( !!newColName := ab_count_normalized)
    CbindedSubsamplingOutput <- merge(x = CbindedSubsamplingOutput, y = tmpSubs, by = c('sample_id','ab_id'), all.x = TRUE)
    i = i+1
  }
  
  normalized_data <- CbindedSubsamplingOutput
  normalized_data$ab_count_normalized <- rowSums(normalized_data[,grepl( "LEAVEONEOUTSUBS" , names( normalized_data ) )], na.rm=TRUE)
  normalized_data <- normalized_data %>%
    select(-all_of(colVector)) %>%
    mutate(ab_count_normalized = ab_count_normalized / (length(colVector)-1))
  
  return(normalized_data)
}

run_libsize_normalization<-function(data)
{
  libsize <- data %>% 
    group_by(sample_id) %>% 
    mutate(lib_size = sum(ab_count)) %>% 
    select(sample_id, lib_size) %>% 
    ungroup() %>% 
    distinct()
  normalized_data <- data %>% 
    left_join(libsize, by = "sample_id") %>%  
    mutate(ab_count_normalized = ab_count/(lib_size))
  
  return(normalized_data)
}

run_normalization<-function(dataset, method)
{
  datapath<-get_dataset_path()
  datapath<-paste(datapath, "/", dataset, sep="")
  data<-read_tsv(datapath, col_types = cols())
  dir.create(paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset)), showWarnings = FALSE)
  
  columns = get_dataset_specific_column_names(datapath)
  dataset_processing_variables = get_prefilter_thresholds(datapath)
  
  data <- data %>%
    rename(
      sample_id=columns[1],
      ab_id=columns[2],
      ab_count=columns[3],
      batch_id=columns[4],
      cluster_id=columns[5],
      ab_type=columns[6]
    ) %>%
    prefilter_dataset(dataset_processing_variables)
  
  print(paste("RUN NORMALIZATION[", method, "] FOR DATASET[", datapath, "]", sep = ""))
  if(method=="TMM")
  {
    normalized_data<- data %>% 
      run_tmm() 
    if(dataset_processing_variables[3] == 1)
    {
      normalized_data <- remove_batch_effect(normalized_data, log_transform = FALSE)
    }
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/TMM.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="SUBSAMPLING")
  {
    normalized_data<- data %>% 
      run_subsampling()
    if(dataset_processing_variables[3] == 1)
    {
      normalized_data <- remove_batch_effect(normalized_data, log_transform = FALSE)
    }
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/SUBSAMPLED.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="LIBSIZE")
  {
    normalized_data<- data %>% 
      run_libsize_normalization() 
    if(dataset_processing_variables[3] == 1)
    {
      normalized_data <- remove_batch_effect(normalized_data, log_transform = FALSE)
    }
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/LIBSIZE.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="SCTRANSFORM")
  {
    #perform batch correction within SCTRANSFORM
    batchEffect = FALSE
    if(dataset_processing_variables[3] == 1)
    {
      batchEffect = TRUE
    }
    normalized_data<- data %>% 
      run_sctransform(batchEffect) 
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/LIBSIZE.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="CLR_COMPOSITIONS")
  {
    #CLR is alraedy a log transformation
    normalized_data<- data %>% 
      run_clr()
    if(dataset_processing_variables[3] == 1)
    {
      normalized_data <- remove_batch_effect(normalized_data, log_transform = FALSE)
    }
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/CLR.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="CLR_SEURAT")
  {
    #CLR is alraedy a log transformation
    normalized_data<- data %>% 
      run_clr_seurat()
    if(dataset_processing_variables[3] == 1)
    {
      normalized_data <- remove_batch_effect(normalized_data, log_transform = FALSE)
    }
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/CLR_SEURAT.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="LEAVE_ONE_OUT_TMM")
  {
    normalized_data<- data %>% 
      run_leave_one_out_tmm() 
    if(dataset_processing_variables[3] == 1)
    {
      normalized_data <- remove_batch_effect(normalized_data, log_transform = FALSE)
    }
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/LeaveOneOutTMM.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="LEAVE_ONE_OUT_SUBSAMPLING")
  {
    normalized_data<- data %>% 
      run_leave_one_out_subsampling() 
    if(dataset_processing_variables[3] == 1)
    {
      normalized_data <- remove_batch_effect(normalized_data, log_transform = FALSE)
    }
    output_table<-paste0(here("bin/NORMALIZED_DATASETS/"), tools::file_path_sans_ext(dataset), "/LeaveOneOutSubsampling.tsv")
    write_tsv(normalized_data, file=output_table)
  }
  else if(method=="EXPORT_FILTERED_DATA")
  {
    output_table<-paste0(here("bin/FILTERED_DATASETS/"), tools::file_path_sans_ext(dataset), ".tsv")
    write_tsv(data, file=output_table)
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
