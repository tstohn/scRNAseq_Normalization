library(stringr)

get_dataset_path<-function()
{
  setting_lines <- readLines("settings.ini")
  matches <- str_match(setting_lines,"datasets=(.*)")
  path=na.omit(as.vector(matches))
  if(!assertthat::are_equal(length(path), 2))
  {
    print("More than one datasets line are given in settings.ini")
  }
  return(path[2])
}

get_dataset_specific_column_names<-function(path)
{

  setting_lines <- readLines(paste0(tools::file_path_sans_ext(path),".ini"))
  
  generic_columns<-c("sample_id","ab_id","ab_count","batch_id","cluster_id","ab_type")
  specific_columns<-c()
  count<-1
  for(col in generic_columns)
  {
    matches <- str_match(setting_lines,paste0(col,"=(.*)"))
    path=na.omit(as.vector(matches))
    if(!assertthat::are_equal(length(path), 2))
    {
      print(paste0("More than one datasets line are given in ini file of dataset"))
    }
    specific_columns[count] = path[2]
    count<-count+1
  }
  
  return(specific_columns)
}

get_prefilter_thresholds<-function(path)
{
  
  setting_lines <- readLines(paste0(tools::file_path_sans_ext(path),".ini"))
  
  generic_columns<-c("threshold_minTotalABCountForSample","threshold_minUMIAB", "removeBatchEffect")
  thresholds<-c()
  count<-1
  for(col in generic_columns)
  {
    matches <- str_match(setting_lines,paste0(col,"=(.*)"))
    path=na.omit(as.vector(matches))
    if(!assertthat::are_equal(length(path), 2))
    {
      print(paste0("More than one datasets line are given in ini file of dataset"))
    }
    thresholds[count] = path[2]
    count<-count+1
  }
  
  return(thresholds)
}
