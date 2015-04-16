# Read in the counts.csv file and log2 normalise against GAPDH
# Will only work for mir342 data sets


get_all_raw_counts <- function (counts_file_list) {
  library('nesoni')
  library('reshape')
  library('ggplot2')
  
  both_counts  <- data.frame()
  
  for (counts_file in counts_file_list){
    all_counts  <- read.grouped.table(counts_file)
    counts <- all_counts$Count
    if (ncol(counts) ==10){
      counts_gene_expession <- counts[c(1,3:8),5:10]
    }
    else{
      counts_gene_expession <- counts[c(1,3:8),]
    }
    
    #total_norm_counts<- counts/colSums(counts)
    print(counts)
    
    
    melted_all_counts <- melt(as.matrix(counts_gene_expession ))
    # print(melted_all_counts)
    colnames(melted_all_counts) <- c('gene', 'sample', 'value')
    #normalise for total 
    #  gapdh <- (melted_all_counts[melted_all_counts$gene =='GAPDH',])
    # print(qplot(data = gapdh, x= sample, y= value )+ geom_point() +theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    print(ggplot(data= melted_all_counts, aes(x= gene ,y= value))+ geom_point()+ggtitle('PAT-Seq Raw Counts Testing Gene Expression')+ facet_wrap(~ sample, ncol =5) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    
  }
}
  

make_normalised_counts_frame <- function(counts_list, tail = 'no_tail', gene_expression = T, raw = F){

  
  both_normalised_counts <- data.frame()
  
  for (counts_file in counts_list){
    all_counts  <- read.grouped.table(counts_file)
    counts <- vector()
    if (tail == 'count')
    {
      counts <- all_counts$Tail_count
    }
    else if (tail == 'length'){
      counts <- all_counts$Tail
    }
    else{
      counts <- all_counts$Count
    }
    if (ncol(counts) ==10){
      counts <- counts[,5:10]
    }
    else{
      counts <- counts
    }
    
    
    if (raw != T){
      
      
      norm_count <- data.frame(lapply(counts, function(x) log2((x +0.5)/(x[42]+0.5))))
      rownames(norm_count) <- rownames(counts)
      counts<- norm_count
    }
    
    if (gene_expression == T){
      counts <- counts[1:8,]
    }
    
    
    melted_all_counts_replace <- melt(as.matrix(counts))
    colnames(melted_all_counts_replace) <- c('gene', 'sample', 'value')
    melted_all_counts <- melt(as.matrix(counts))
    colnames(melted_all_counts) <- c('gene', 'sample', 'value')
    melted_all_counts$sample <- melted_all_counts_replace$sample
    melted_all_counts <- cbind(melted_all_counts, rep (paste(counts_file), nrow(melted_all_counts)))
    colnames(melted_all_counts) <- c('gene', 'sample', 'value', 'experiment')
    
    both_normalised_counts  <- rbind(both_normalised_counts, melted_all_counts)
  }
  
  splitframes <- split (both_normalised_counts, both_normalised_counts$experiment)
  near_perfect_frame <- cbind(splitframes[[1]], splitframes[[2]])
  near_perfect_frame <- near_perfect_frame[,-5]
  near_perfect_frame[,4] <- rep('PAT-Seq', nrow(near_perfect_frame))
  near_perfect_frame[,7] <- rep('Re-PAT', nrow(near_perfect_frame))
  return (near_perfect_frame)
  
}

plot_counts <- function(counts_frame, name){
  print(ggplot(data= counts_frame, aes(x= value ,y= value.1))+ geom_point()+
          labs(x = 'PAT-Seq', y = 'RePat')+
          
          ggtitle(paste(name)))
  
}

compare_counts <- function (counts_file_list) {
  library('nesoni')
  library('reshape')
  library('ggplot2')
  
  both_normalised_counts  <- data.frame()
  
  for (counts_file in counts_file_list){
    all_counts  <- read.grouped.table(counts_file)
    counts <- all_counts$Count
    if (ncol(counts) ==10){
      counts <- counts[,5:10]
    }
    else{
      counts <- counts
    }
    
    normalised_frame <- data.frame(row.names = rownames (counts))
    for (col in counts){
      normalised_frame <- cbind(normalised_frame,log2(((col+0.5)/((col[42]+0.5)))) )
    } 
    melted_all_counts_replace <- melt(as.matrix(counts))
    colnames(melted_all_counts_replace) <- c('gene', 'sample', 'value')
    melted_all_counts <- melt(as.matrix(normalised_frame))
    colnames(melted_all_counts) <- c('gene', 'sample', 'value')
    melted_all_counts$sample <- melted_all_counts_replace$sample
    melted_all_counts <- cbind(melted_all_counts, rep (paste(counts_file), nrow(melted_all_counts)))
    colnames(melted_all_counts) <- c('gene', 'sample', 'value', 'experiment')
    
    
    #normalise for total 
    #  gapdh <- (melted_all_counts[melted_all_counts$gene =='GAPDH',])
    # print(qplot(data = gapdh, x= sample, y= value )+ geom_point() +theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    print(ggplot(data= melted_all_counts, aes(x= gene ,y= value))+ geom_point()+ggtitle('PAT-Seq Normalised to GAPDH Counts')+ facet_wrap(~ sample, ncol =3) + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
    both_normalised_counts  <- rbind(both_normalised_counts, melted_all_counts)
  }
  print (both_normalised_counts) 
  print(ggplot(data=both_normalised_counts, aes(x= value[1:144], y = value[289:432] ) ) + geom_point() +ggtitle('Normalised Counts:\n Empty vs Empty Counts PAT-Seq v RePAT')+labs(y= 'RePAT', x = 'PAT-Seq'))
  print(ggplot(data=both_normalised_counts, aes(x= value[145:288], y = value[433:576] ) ) + geom_point() +ggtitle('Normalised Counts:\n  Mir342 vs Mir342 PAT-Seq v RePAT')+labs(y= 'RePAT', x = 'PAT-Seq'))
  print(ggplot(data=both_normalised_counts, aes(x= value[1:144], y = value[145:288] ) ) + geom_point() +ggtitle('Normalised Counts:\n  Empty vs mir342 PAT-Seq')+labs(x= 'Empty', y = 'mir342'))
  print(ggplot(data=both_normalised_counts, aes(x= value[289:432], y = value[433:576] ) ) + geom_point() +ggtitle('Normalised Counts:\n Empty vs mir342 Re-PAT')+labs(x= 'Empty', y = 'mir342'))
  
  # Get fold change test - control 
  fold_change_patseq <-  both_normalised_counts$value[145:288] - both_normalised_counts$value[1:144]
  fold_change_patseq <- cbind(fold_change_patseq,  rep('patseq',length(fold_change_patseq) ))
  colnames(fold_change_patseq) <- c('experiment', 'FC')
  
  # Get fold change test - control 
  fold_change_repat <-both_normalised_counts$value[433:576]- both_normalised_counts$value[289:432] 
  fold_change_repat <- cbind(fold_change_repat, rep('repat',length(fold_change_repat) ))
  
  
  colnames(fold_change_repat) <- c('experiment', 'FC')
  
  
  both_fc <- as.data.frame(rbind(fold_change_patseq, fold_change_repat))
  
  print(ggplot (data = both_fc, aes(x= experiment, y = FC) ) +geom_boxplot())
  
}

counts_boxplots <- function(counts_frame, name){
  print(ggplot(data= counts_frame, aes(factor(sample), value))+ geom_boxplot()+
          labs(x = 'Sample', y = 'Normalised Counts')+
          ggtitle('PAT-Seq'))
  
  print(ggplot(data= counts_frame, aes(factor(sample.1) ,value.1))+ geom_boxplot()+
          labs(x = 'Sample', y = 'Normalised Counts')+
          ggtitle('Re-PAT'))
}
multivars <- function(counts_frame, experiment = 'PAT-Seq', control = T) {
  
  ordered_by_gene <- counts_frame[with(counts_frame, order(gene)),]
  split_by_gene <- (split(ordered_by_gene, ordered_by_gene$gene))
  
  control_variances <- vector() 
  test_variances <- vector() 
  total_variances<- vector()
  
  control_sum <- vector() 
  test_sum <- vector() 
  total_sum<- vector()
  
  for (df in split_by_gene){
    if (experiment == 'PAT-Seq'){
      
      control_variances <-c(control_variances, var(df$value[1:3]))
      test_variances <-c(test_variances, var(df$value[4:6]))
      total_variances<-c(total_variances, var(df$value))
      
      control_sum <-c(control_sum, sum(df$value[1:3]))
      test_sum <-c(test_sum, sum(df$value[4:6]))
      total_sum<-c(total_sum, sum(df$value))
    }
    else{
      control_variances <-c(control_variances, var(df$value.1[1:3]))
      test_variances <-c(test_variances, var(df$value.1[4:6]))
      total_variances<-c(total_variances, var(df$value.1))
      
      control_sum <-c(control_sum, sum(df$value.1[1:3]))
      test_sum <-c(test_sum, sum(df$value.1[4:6]))
      total_sum<-c(total_sum, sum(df$value.1))
    }
  }
  
  names(control_variances) <- counts_frame$gene[1:length(control_variances)]
  names(test_variances) <- counts_frame$gene[1:length(control_variances)]
  names(total_variances) <- counts_frame$gene[1:length(control_variances)]
  
  
  control_df <- as.data.frame(cbind(control_variances, control_sum))
  test_df <- as.data.frame(cbind(test_variances, test_sum))
  total_df <- as.data.frame(cbind(total_variances,total_sum))
  
  print (ggplot(control_df, aes(x= control_sum, y= control_variances)) + geom_point()+ggtitle (paste(experiment,'Empty Variance vs Count Following Log2 Normalisation to GAPDH') ))
  print (ggplot(test_df, aes(x= test_sum, y= test_variances)) + geom_point() +ggtitle (paste(experiment,'mir342 Variance vs Count Following Log2 Normalisation to GAPDH')) )
  print (ggplot(total_df, aes(x= total_sum, y= total_variances)) + geom_point() +ggtitle (paste(experiment,'All Variance vs Count Following Log2 Normalisation to GAPDH')) )
  
  if (control == T){
    return(control_variances)
  }
  else{
    return (test_variances)
  }
}

plot_vars <- function (var1, var2, title) {
  both <- as.data.frame (cbind(var1, var2))
  print(both)
  print (ggplot(both, aes(x= log2(var1), y= log2(var2))) + geom_point() +ggtitle (paste(title)) +labs(x= 'PAT-Seq', y = 'Re-PAT'))
}
box_plot_FC <- function (counts_frame){
  pat_seq <- counts_frame [,c(1,2:4)]
  re_pat <-  counts_frame [,c(1,5:7)]
  colnames (re_pat)<- colnames (pat_seq)
  plot_frame <- rbind(pat_seq,re_pat)
  
  
  by_sample <- split(plot_frame , plot_frame$sample) 
  
  
  pat_seq_control <- cbind(by_sample$EMPTYrep1$value, by_sample$EMPTYrep2$value, by_sample$EMPTYrep3$value)
  pat_seq_mir342 <-  cbind(by_sample$miR242rep1$value, by_sample$miR242rep2$value, by_sample$miR242rep13$value)
  
  pat_seq_control_means <- rowMeans(pat_seq_control)
  pat_seq_mir342_means <- rowMeans(pat_seq_mir342)
  
  repat_control <- cbind(by_sample$MDA.EV.rep1$value, by_sample$MDA.EV.rep2$value, by_sample$MDA.EV.rep3$value)
  repat_mir342 <-  cbind(by_sample$MDA.mir342.rep1$value, by_sample$MDA.mir342.rep2$value, by_sample$MDA.mir342.rep3$value)
  
  repat_control_means <- rowMeans(repat_control)
  repat_mir342_means <- rowMeans(repat_mir342)
  
  mean_FC_pat <- pat_seq_mir342_means - pat_seq_control_means
  mean_FC_repat <- repat_mir342_means - repat_control_means
  
  mean_FC_pat <- cbind(mean_FC_pat, rep('PAT-Seq', length(mean_FC_pat)) )
  mean_FC_repat <- cbind(mean_FC_repat, rep('Re-PAT', length(mean_FC_repat)) )
  fold_changes <- as.data.frame (rbind( mean_FC_pat, mean_FC_repat))
  colnames(fold_changes)<- c('Mean_Fold_Change', 'Experiment')
  fold_changes$Mean_Fold_Change <- as.numeric(fold_changes$Mean_Fold_Change)
  
  print(ggplot(data= fold_changes, aes(x= Experiment, y = Mean_Fold_Change))+geom_boxplot()+ labs (y = 'Log2 Normalised Counts\n Fold Changes')+
          ggtitle ('Normalised Fold Changes Boxplot'))
  print(fold_changes)
  
}


compare_experiments <- function(counts_list){
  
}