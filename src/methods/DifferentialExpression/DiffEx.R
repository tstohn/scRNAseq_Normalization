#Differential Expression Script
library(ggplot2)
library(gplots)

options(tidyverse.quiet = TRUE)
library(tidyverse)
library(edgeR)
library(scales)
library(dplyr)

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

#PREPROCESSING FOR ALL OTHER SUB ELEMENTS

getDataOnAllBCs<-function()
{
  data_ab <- read_tsv("/Users/t.stohn/Desktop/KATHY/RUN7/ABBarcodeProcessing.tsv")
  
  sumData <- data_ab %>% 
    group_by(AB_BARCODE, TREATMENT) %>%
    summarise(AB_COUNT = sum(AB_COUNT)/n()) %>%
    ungroup() %>%
    unique()
  sumData <- unique(sumData)
}

getDataOnAllBCs90Perc<-function()
{
  data_ab <- read_tsv("/Users/t.stohn/Desktop/ABBarcodeProcessing90PercUmi.tsv")
  
  return(data_ab)
}

getOldData<-function()
{
  data_ab <- read_tsv("/Users/t.stohn/Desktop/KATHY/RUN7/ABBarcodeProcessing.tsv")
  
  return(data_ab)
}

getDataOnBC1<-function()
{
  data_ab <- read_tsv("/Users/t.stohn/Desktop/KATHY/RUN7/ABBarcodeProcessing.tsv")
  
  #keep only BC1 as index
  REGEXSTRING <- "\\d+\\.(\\d+)\\.\\d+\\.\\d+"
  dataColumn <- str_match(data_ab$SingleCell_BARCODE, REGEXSTRING)
  bcNew <- (dataColumn[,2])
  newData <- data_ab
  newData$newSCBarcode = bcNew
  newData <- newData %>% select(-SingleCell_BARCODE)
  sumData <- newData %>% 
    group_by(AB_BARCODE, newSCBarcode) %>%
    summarise(AB_COUNT = sum(AB_COUNT)/n(), TREATMENT = TREATMENT)
  sumData <- unique(sumData)
  
  return(sumData)
}

seperateTreatmentDuplicates<-function(data)
{
  data$newSCBarcode <- as.integer(data$newSCBarcode)
  
  #map to treatment duplicates
  bc_map1 <- c(1,2,3,4,5,6,13,14,15,16,17,18,25,26,27,28,29,30,37,38,39,40,41,42,49,50,51,52,53,54,61,62,63,64,65,66,73,74,75,76,77,78,85,86,87,88,89,90)
  treatment_list <- c("u_1","u_2","A_1","A_2","B_1","B_2","u_1","u_2","A_1","A_2","B_1","B_2","u_1","u_2","A_1","A_2","B_1","B_2","u_1","u_2","A_1","A_2","B_1","B_2","u_1","u_2","A_1","A_2","B_1","B_2","u_1","u_2","A_1","A_2","B_1","B_2","u_1","u_2","A_1","A_2","B_1","B_2","u_1","u_2","A_1","A_2","B_1","B_2")
  map1 <- data.frame(realID=treatment_list, id=c(0:47) )

  x <- data %>%
    left_join(map1, by=c("newSCBarcode" = "id")) %>%
    mutate(newTreamtent = realID) %>%
    select(-realID)
  
  
  return(x)
}

getDataOnBC1NormalizedByWellCount<-function()
{
  #data_ab <- read_tsv("/Users/t.stohn/Desktop/KATHY/RUN7/ABBarcodeProcessing.tsv")
  data_ab <- read_tsv("/Users/t.stohn/Desktop/KATHY/RUN7/ABBarcodeProcessing90PercUmi.tsv")
  #data_ab <- read_tsv("/Users/t.stohn/Desktop/ABBarcodeProcessingMoreThan10Umis.tsv")
  
  #list of treatment by barcode, not used, just for one self to know ....
  #u1 <- c(1,13,25,37,49,61,73,85)
  #u2 <- c(2,14,26,38,50,62,74,86)
  #az1 <- c(3,15,27,39,51,63,75,87)
  #az2 <- c(4,16,28,40,52,64,76,88)
  #bez1 <- c(5,17,29,41,53,65,77,89)
  #bez2 <- c(6,18,30,42,54,66,78,90)
    
  #keep only BC1 as index
  REGEXSTRING <- "\\d+\\.(\\d+)\\.\\d+\\.\\d+"
  dataColumn <- str_match(data_ab$SingleCell_BARCODE, REGEXSTRING)
  bcNew <- (dataColumn[,2])
  
  newData <- data_ab
  newData$newSCBarcode = bcNew
  newData <- newData %>% select(-SingleCell_BARCODE)
  newData <- seperateTreatmentDuplicates(newData)
  
  #remove some strange single cells
  #newData <- newData[newData$AB_COUNT > 4 & newData$AB_COUNT < 100,]

  sumData <- newData %>% 
    group_by(newSCBarcode) %>%
    mutate(totalBC1Occurences = sum(AB_COUNT)) %>%
    ungroup() %>%
    group_by(AB_BARCODE, newSCBarcode) %>%
    summarise(AB_COUNT = sum(AB_COUNT), TREATMENT = newTreamtent, totalBC1Occurences =totalBC1Occurences) %>%
    ungroup() %>%
    mutate(AbCountNormalized = AB_COUNT/totalBC1Occurences) %>%
    unique()
  sumData <- unique(sumData)
  
  return(sumData)
}

getDataOnBC1NormalizedByWellCount_2<-function()
{

  data_ab <- read_tsv("/Users/t.stohn/Desktop/ABBarcodeProcessing90PercUmi.tsv")
  
  data_ab <- data_ab %>%
    group_by(SingleCell_BARCODE) %>%
    mutate(ReadsPerSingleCell = sum(AB_COUNT))%>%
    ungroup()
  data_ab <- data_ab[data_ab$ReadsPerSingleCell > 1,]
  #keep only BC1 as index
  REGEXSTRING <- "\\d+\\.(\\d+)\\.\\d+\\.\\d+"
  dataColumn <- str_match(data_ab$SingleCell_BARCODE, REGEXSTRING)
  bcNew <- (dataColumn[,2])
  newData <- data_ab
  newData$newSCBarcode = bcNew
  newData <- newData %>% select(-SingleCell_BARCODE)
  
  #remove some strange single cells
  #newData <- newData[newData$AB_COUNT > 4 & newData$AB_COUNT < 100,]
  
  sumData <- newData %>% 
    group_by(newSCBarcode) %>%
    mutate(totalBC1Occurences = sum(AB_COUNT)) %>%
    ungroup() %>%
    group_by(AB_BARCODE, newSCBarcode) %>%
    summarise(AB_COUNT = sum(AB_COUNT), TREATMENT = TREATMENT, totalBC1Occurences =totalBC1Occurences) %>%
    ungroup() %>%
    mutate(AbCountNormalized = AB_COUNT/totalBC1Occurences) %>%
    unique()
  sumData <- unique(sumData)
  
  return(sumData)
}

getDataWithCountsForEachAB<-function()
{
  data_RUN7 <- read_tsv("/Users/t.stohn/Desktop/KATHY/RUN7/ABBarcodeProcessing.tsv")
  #to get only high quality data
  abCount <- data_RUN7 %>%
    group_by(SingleCell_BARCODE) %>%
    mutate(occurences = n())
  abCount <- abCount[abCount$occurences == 30,]
  abCountGood <- abCount %>%
    group_by(SingleCell_BARCODE) %>%
    mutate(sum = sum(AB_COUNT)) %>%
    ungroup()
  abCountGood <-abCountGood[abCountGood$sum > 200 & abCountGood$sum < 1000,] #define here what you want
  abCountGood <- abCountGood %>% select(-sum, -occurences)
}
    

#HEAT MAP OF ANALYSIS
# make matrix of Barcode * treatment
##########################################################
# ANALYSIS ON CELLS DEFINED BY BC1
##########################################################
  
  data_ab <- getDataOnBC1()
  
  hist(data_ab$AB_COUNT, breaks = 10000, xlim = c(0,50))
  
  #make a matrix for median of AB VS. TREATMENT
  x <- data_ab %>%
    group_by(AB_BARCODE, TREATMENT) %>%
    mutate(median = median(AB_COUNT)) %>%
    mutate(sum = sum(AB_COUNT)) %>%
    ungroup() %>%
    select(AB_BARCODE, TREATMENT, sum) %>%
    unique() %>%
    pivot_wider(names_from = TREATMENT, values_from = sum) %>%
    column_to_rownames("AB_BARCODE")
  x <- x[c("untreated", "AZD6244", "BEZ235")]
  
  heatmapInput<-as.matrix(x[ order(row.names(x)), ])
  col_panel<-colorpanel(100, low="blue", mid="white", high="yellow")
  heatmap.2(heatmapInput, scale='none', col = cividis::cividis_pal(), trace = "none", dendrogram = "none", cexRow = 0.8, cexCol = 0.9, Rowv=FALSE,
            Colv=FALSE)
  
  #barplot of untreated
  untreated <- x %>%
    select(untreated) %>%
    rownames_to_column("Antibody target")
  untreated <- untreated[order(untreated$untreated, decreasing = TRUE),]
  ggplot(data = untreated, aes(x = reorder(`Antibody target`, -untreated), y = untreated)) + geom_bar(stat='identity') +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + ylab("Count") + ggtitle("Average Count of AB reads for untreated cells") +
    xlab("Antibody target")
  
  
#HEAT MAP OF ANALYSIS
  # make matrix of Barcode * treatment
##########################################################
# ANALYSIS ON CELLS SUMMED OVER ALL BCS; NO NORMALIZATION TOTAL COUNTS
##########################################################
  
  
  #SUM OVER ALL BARCODES
    newData <- getDataOnAllBCs()
    hist(sumData$AB_COUNT, breaks = 10000, xlim = c(0,50))
  
  #make a matrix for median
  x <- newData %>%
    pivot_wider(names_from = TREATMENT, values_from = AB_COUNT) %>%
    column_to_rownames("AB_BARCODE")
  x <- x[c("untreated", "AZD6244", "BEZ235")]
  
  heatmapInput<-as.matrix(x[ order(row.names(x)), ])
  col_panel<-colorpanel(100, low="blue", mid="white", high="yellow")
  heatmap.2(heatmapInput, scale='none', col = cividis::cividis_pal(), trace = "none", dendrogram = "none", cexRow = 0.8, cexCol = 0.9, Rowv=FALSE,
            Colv=FALSE)
  
  untreated <- x %>%
    select(untreated) %>%
    rownames_to_column("Antibody target")
  untreated <- untreated[order(untreated$untreated, decreasing = TRUE),]
  ggplot(data = untreated, aes(x = reorder(`Antibody target`, -untreated), y = untreated)) + geom_bar(stat='identity') +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + ylab("Count") + ggtitle("Average Count of AB reads for untreated cells") +
    xlab("Antibody target")

##########################################################
  # ANALYSIS ON CELLS WITH COUNTS FOR ALL 30 ABS
##########################################################

  data_ab <- getDataWithCountsForEachAB()
  
  #keep only BC1 as index
  sumData <- data_ab %>% 
    group_by(AB_BARCODE, SingleCell_BARCODE) %>%
    summarise(AB_COUNT = sum(AB_COUNT)/n(), TREATMENT = TREATMENT)
  sumData <- unique(sumData)
  hist(sumData$AB_COUNT, breaks = 10000, xlim = c(0,50))
  #make a matrix for median
  qualityData <- sumData %>%
    group_by(AB_BARCODE, TREATMENT) %>%
    mutate(median = median(AB_COUNT)) %>%
    ungroup() %>%
    select(AB_BARCODE, TREATMENT, median) %>%
    unique() %>%
    pivot_wider(names_from = TREATMENT, values_from = median) %>%
    column_to_rownames("AB_BARCODE")
  
  #heatmap
    x <- qualityData[c("untreated", "AZD6244", "BEZ235")]
    heatmapInput<-as.matrix(x)
    col_panel<-colorpanel(1000, low="black", high="yellow")
    heatmap.2(heatmapInput, scale='none', col = cividis::cividis_pal(), trace = "none", dendrogram = "none", cexRow = 0.8, cexCol = 0.9, Rowv=FALSE,
              Colv=FALSE)
  
  #heat after tmm normalization
    abCountNorm <- abCountGood %>%
      mutate(sample_id = SingleCell_BARCODE) %>%
      mutate(ab_id = AB_BARCODE) %>%
      mutate(treatment_id = TREATMENT) %>%
      mutate(ab_count = AB_COUNT) %>%
      select(sample_id, ab_id, treatment_id, ab_count)
    abCountNorm <- abCountNorm %>%
      run_tmm()
    summarizedData <- abCountNorm %>%
      group_by(ab_id, treatment_id) %>%
      summarise(average = mean(ab_count_normalized)) %>%
      pivot_wider(names_from = treatment_id, values_from = average) %>%
      column_to_rownames("ab_id")
    summarizedData <- summarizedData[c("untreated", "AZD6244", "BEZ235")]
    heatmapInputNorm<-as.matrix(summarizedData)
    col_panel<-colorpanel(1000, low="black", high="yellow")
    heatmap.2(heatmapInputNorm, scale='none', col = cividis::cividis_pal(), trace = "none", dendrogram = "none", cexRow = 0.8, cexCol = 0.9, Rowv=FALSE,
              Colv=FALSE)
  
    #barplot of untreated
    untreated <- x %>%
      select(untreated) %>%
      rownames_to_column("Antibody target")
    untreated <- untreated[order(untreated$untreated, decreasing = TRUE),]
    ggplot(data = untreated, aes(x = reorder(`Antibody target`, -untreated), y = untreated)) + geom_bar(stat='identity') +
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1)) + ylab("Count") + ggtitle("Average Count of AB reads for untreated cells") +
      xlab("Antibody target")

##########################################################
# BOXPLOT AND ANOVA ON BC1 WITH GREEN FRAME IF SIGNIFICANT
##########################################################
    
    data_ab <- getDataOnBC1NormalizedByWellCount()
    
    #make a matrix for median of AB VS. TREATMENT
    x <- data_ab %>%
      group_by(AB_BARCODE, TREATMENT) %>%
      mutate(median = median(AbCountNormalized)) %>%
      ungroup() %>%
      select(AB_BARCODE, TREATMENT, median) %>%
      unique() %>%
      pivot_wider(names_from = TREATMENT, values_from = median) %>%
      column_to_rownames("AB_BARCODE")
    x <- x[c("u_1","u_2", "A_1","A_2", "B_1","B_2")]
    
    heatmapInput<-as.matrix(x[ order(row.names(x)), ])
    col_panel<-colorpanel(100, low="blue", mid="white", high="yellow")
    heatmap.2(heatmapInput, scale='none', col = cividis::cividis_pal(), trace = "none", dendrogram = "none", cexRow = 0.8, cexCol = 0.9, Rowv=FALSE,
              Colv=FALSE)
    
    ABS <- unique(data_ab$AB_BARCODE)
    plotsList <- list()
    for(i in 1:length(ABS))
    {
      #get data
      ab <- ABS[[i]]
      dataTmp <- data_ab[data_ab$AB_BARCODE == ab,]
      dataTmp <- dataTmp %>%
        mutate(TREATMENT = fct_relevel(TREATMENT, 
                                       "u_1",
                                       "u_2",
                                       "A_1",
                                       "A_2",
                                       "B_1",
                                       "B_2"))
      #make Anova
      res.aov <- aov(AbCountNormalized ~ TREATMENT, data = dataTmp)
      # Summary of the analysis
      summary(res.aov)
      val <- round(summary(res.aov)[[1]][["Pr(>F)"]][1], digits = 3)
      
      themeNormal <- theme(panel.border = element_rect(color = "black",
                                                 fill = NA,
                                                 size = 0.5))
      themeHighlight <- theme(panel.border = element_rect(color = "green",
                                                       fill = NA,
                                                       size = 1))
      #store plot
      if(val <= 0.05)
      {
        themeVariable <- themeHighlight
      }
      else
      {
        themeVariable <- themeNormal
      }
        p <- ggplot(data = dataTmp, aes(x = TREATMENT, y = AbCountNormalized)) +geom_boxplot() + geom_jitter(width = 0.1) + xlab(paste0(ABS[i], " | p-val: ", val)) +
          themeVariable + ylab("AB count")
        plotsList[[i]] <- p
    }
    library("gridExtra")
    library(grid)
    grid_plot <- grid.arrange(grobs = plotsList, ncol = 5, nrow = 6, 
                              top = textGrob("Boxplots of AB count per treatment: \neach datapoint is an AB count for a specific well, \nnormalized by the total number of counts over all ABs in this well",
                              gp=gpar(fontsize=20,font=3)))
    
    ABS <- c("pS6","pAKT","pERK1/2","p4EBP1")
    plotsList <- list()
    for(i in 1:length(ABS))
    {
      #get data
      ab <- ABS[[i]]
      dataTmp <- data_ab[data_ab$AB_BARCODE == ab,]
      dataTmp <- dataTmp %>%
        mutate(TREATMENT = fct_relevel(TREATMENT, 
                                       "u_1",
                                       "u_2",
                                       "A_1",
                                       "A_2",
                                       "B_1",
                                       "B_2"))
      #make Anova
      res.aov <- aov(AbCountNormalized ~ TREATMENT, data = dataTmp)
      # Summary of the analysis
      summary(res.aov)
      val <- round(summary(res.aov)[[1]][["Pr(>F)"]][1], digits = 3)
      
      themeNormal <- theme(panel.border = element_rect(color = "black",
                                                       fill = NA,
                                                       size = 0.5))
      themeHighlight <- theme(panel.border = element_rect(color = "green",
                                                          fill = NA,
                                                          size = 1))
      #store plot
      if(val <= 0.05)
      {
        themeVariable <- themeHighlight
      }
      else
      {
        themeVariable <- themeNormal
      }
      p <- ggplot(data = dataTmp, aes(x = TREATMENT, y = AbCountNormalized)) +geom_boxplot() + geom_jitter(width = 0.1) + xlab(paste0(ABS[i], " | p-val: ", val)) +
        themeVariable + ylab("AB count")
      plotsList[[i]] <- p
    }
    library("gridExtra")
    library(grid)
    grid_plot <- grid.arrange(grobs = plotsList, ncol = 2, nrow = 2, 
                              top = textGrob("AB count based on BC1",
                                             gp=gpar(fontsize=20,font=3)))
    
    #LOOK AT ONE SPECIFIC EXAMPLE
      dataTmp <- data_ab[data_ab$AB_BARCODE == "pAKT",]
      ggplot(data = dataTmp, aes(x = TREATMENT, y = AbCountNormalized)) +geom_boxplot() + geom_jitter(width = 0.1)
      # Compute the analysis of variance
      res.aov <- aov(AbCountNormalized ~ TREATMENT, data = dataTmp)
      # Summary of the analysis
      summary(res.aov)
      val <- summary(res.aov)[[1]][["Pr(>F)"]][1]
  
##########################################################
# ONLY 90 perc UMIs kept - no more insight
##########################################################

  data <- getDataOnAllBCs90Perc()
  oldData <- getOldData()
  
  data <- data %>%
    group_by(SingleCell_BARCODE) %>%
    summarize(ReadsPerSingleCell = sum(AB_COUNT))%>%
    ungroup()
  
  oldData <- oldData %>%
    group_by(SingleCell_BARCODE) %>%
    summarize(ReadsPerSingleCell = sum(AB_COUNT))%>%
    ungroup()
  
  hist(data$ReadsPerSingleCell, breaks = 8000, xlim = c(0,10))
  hist(oldData$ReadsPerSingleCell, breaks = 80000, xlim = c(0,10))
  
  x<- data[data$ReadsPerSingleCell>1,]
  
##########################################################
# ANALYSIS ONLY ON CELLS WITH MANY READS
##########################################################
  #keep only cells with > 200 protein counts
  xcells <- data_RUN7 %>%
    group_by(SingleCell_BARCODE) %>%
    mutate(sum = sum(AB_COUNT)) %>%
    mutate(diffAbSum = n()) %>%
    select(SingleCell_BARCODE, sum, diffAbSum) %>%
    ungroup() %>%
    unique()
  best_10_thous_ells <- xcells[xcells$sum>100 & xcells$diffAbSum>10,]
  cellData <- data_RUN7[data_RUN7$SingleCell_BARCODE %in% best_10_thous_ells$SingleCell_BARCODE,]
  coutTopCells <- cellData %>%
    group_by(SingleCell_BARCODE) %>%
    mutate(sum = sum(AB_COUNT)) %>%
    select(SingleCell_BARCODE, sum) %>%
    ungroup() %>%
    unique()
  hist(coutTopCells$sum, breaks = 100000, xlim=c(0,1000))
  hist(cellData[cellData$AB_BARCODE=="GAPDH",]$AB_COUNT, breaks = 10000, xlim=c(0,200))
  
  cellData <- cellData %>% 
    group_by(SingleCell_BARCODE) %>%
    mutate(sum = sum(AB_COUNT)) %>%
    ungroup() %>%
    mutate(norm_count = AB_COUNT/sum)
  
  #make a matrix for median
  interestingProteins <- cellData[cellData$AB_BARCODE == "pS6" | cellData$AB_BARCODE == "pAKT" | cellData$AB_BARCODE == "pERK1/2" |cellData$AB_BARCODE == "p4EBP1",]
  qualityData <- interestingProteins %>%
    group_by(AB_BARCODE, TREATMENT) %>%
    mutate(median = median(norm_count)) %>%
    ungroup() %>%
    select(AB_BARCODE, TREATMENT, median) %>%
    unique() %>%
    pivot_wider(names_from = TREATMENT, values_from = median) %>%
    column_to_rownames("AB_BARCODE")
  #heatmap
  x <- qualityData[c("untreated", "AZD6244", "BEZ235")]
  heatmapInput<-as.matrix(x)
  col_panel<-colorpanel(1000, low="black", high="yellow")
  heatmap.2(heatmapInput, scale='none', col = cividis::cividis_pal(), trace = "none", dendrogram = "none", cexRow = 0.8, cexCol = 0.9, Rowv=FALSE,
            Colv=FALSE)
  
  #boxplot for a protein
  cellData$AB_BARCODE <- gsub('\\s+', '', cellData$AB_BARCODE)
  dataTmp <- cellData[cellData$AB_BARCODE == "pS6",]
  dataTmp<-dataTmp[dataTmp$AB_COUNT<50,]
  ggplot(data = dataTmp, aes(x = TREATMENT, y = AB_COUNT)) +geom_boxplot() + geom_jitter(width = 0.05, alpha = 0.1)

  cellData$AB_BARCODE <- gsub('\\s+', '', cellData$AB_BARCODE)
  ggData <- cellData[cellData$AB_BARCODE=="p4EBP1",]
  ggplot(ggData, aes(x = norm_count, fill = TREATMENT)) + 
    geom_histogram(bins = 70, position = "identity", alpha = 0.8) + 
    theme_minimal() + ggtitle("Guide Read Count Distribution") + xlab("Number of reads per single cell") +
    ylab("Frequency") + xlim(0,0.3)
  
  
  
  
  
