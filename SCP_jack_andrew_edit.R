

### This R script takes a list of two bam files, a gff file and a gene name or peak number for a gene of interest and returns
### a survival curve of poly A tail distribution over the samples (up to six at a time) based on peaks called in the gff file


SCP <- function(bam_list, gff, name, housekeeping_gene=FALSE, normalisation_factor = 1,  number_of_replicates = 3, legend = FALSE){
  library(Rsamtools)
  

  gff_file <- read.delim(gff, header=FALSE, comment.char="")
  
  #Prints out the order files were processed (the order of the input list)
  print (c('Order of bam file processing', bam_list))
  
  # Pull the poly A reads from each bam file, calls poly_A_puller
  processed_bam_files <- lapply(bam_list , poly_A_puller, gff_file, name, peak)
  print(str(processed_bam_files))
  
  conditions_list <- tapply(processed_bam_files, (seq_along(processed_bam_files)-1) %/% number_of_replicates,list)
  new_list <- list()
  for (reps in conditions_list){
    new_list <- c (new_list, list(unlist(reps, recursive = F)))
    }
  
  new_list <- lapply(new_list, sort)
  print(str(conditions_list))
  print(str(new_list))

  if (housekeeping_gene != FALSE){
    norm_factor <-number_of_reads1/number_of_reads2 
    cat ('normalisation factor(condition 1/ condition 2) =', norm_factor)
  }
  
  norm_reads_1 <-numeric()
  norm_reads_2 <-numeric()
  if (normalisation_factor != 1){
    number_of_reads2 <- number_of_reads2* normalisation_factor
    norm_reads_1 <- toString(c('normalised reads', paste(bam_list[[1]]),number_of_reads1))
    norm_reads_2 <- toString(c('normalised reads', paste(bam_list[[2]]),number_of_reads2))
    cat(norm_reads_1,'\n')
    cat(norm_reads_2)
  }
  
  # Make a reverse cumulative distribuition plot of the poly A reads
  elis <- lapply (new_list, ecdf)
  the_rest <- elis[c(2:length(bam_list))]
  ry <- lapply(new_list,max)
  r <- range(0,150)
  
  # Make graph have no box
  par(bty="l")
  # Statement to determine title of plot, plots first curve. 
  
  curve((1-elis[[1]](x))*100, from=r[1], to=r[2], col=1 , xlim=r, ylab= 'Percent Population (%)', xlab = 'Poly A length', main= paste(name),axes=FALSE)   
  
  # Add a few more axis options
  axis(1, pos=0)
  axis(2, pos=0, at= c(0,25,50,75,100), tick = 25)
  # Plots all the remaining curves 
  peak_plot <- function(ecdf_value, collist,r=r ){
    
    curve((1-ecdf_value(x))*100,  from=r[1], to=r[2], col= collist, xlim=r, ylab= 'Percent Population (%)', xlab = 'Poly A length', add=TRUE)
    
    
  }
  
  
  # Calls the peak plot function with colour parameters
  collist <- rainbow(15)
  mapply (peak_plot, the_rest, collist, MoreArgs= list(r=r))
  
  if(length(bam_list) >2){
    group_mean1<- mean(unlist(conditions_list[[1]]))
    groupsd1 <- sd(unlist(conditions_list[[1]]))
    list1 <- unlist(conditions_list[[1]])
    list2 <- unlist(conditions_list[[2]])
    se <- function(x) sqrt(var(x)/length(x))
    se1  <- se(list1)
    se2 <- se(list2)
    group_mean2<- mean(unlist(conditions_list[[2]]))
    groupsd2 <- sd(unlist(conditions_list[[2]]))
    
    segments(x0 =group_mean1 , y0 = 0, x1= group_mean1, y1 = 100,lty=3,col='red')
    segments(x0 =group_mean2 , y0 = 0, x1= group_mean2, y1 = 100,lty=3,col='blue')
    segments(x0 =group_mean1 , y0 = 100, x1= group_mean2, y1 = 100)
    #legend('topright', bty = "n", legend = c ('D=', paste(round(group_mean1-group_mean2,2)), '±', paste(round( se1+se2,2))))
    file <- 'PAT_seq.txt'
    
    if(name!= FALSE){
      write_matrix <- matrix(c(paste(name),paste(name),'N2', 'Gld-2', as.numeric(group_mean1), as.numeric(group_mean2)),nrow=2,ncol=3) 
    }
    else{
      write_matrix <- matrix(c(paste(peak),paste(peak),'N2', 'Gld-2', as.numeric(group_mean1), as.numeric(group_mean2)),nrow=2,ncol=3)   
    }
    
    
    write.table(write_matrix,file,append=TRUE,row.names=FALSE, col.names=FALSE,quote=FALSE)       
    
    
  }
  
  
  else {
    means <- lapply(processed_bam_files, mean)
    group_mean1<- as.numeric(means[1])
    group_mean2<- as.numeric(means[2])
    segments(x0 =group_mean1 , y0 = 100, x1= group_mean2, y1 = 100)
    
    segments(x0 =group_mean1 , y0 = 0, x1= group_mean1, y1 = 100,lty=3,col='red')
    segments(x0 =group_mean2 , y0 = 0, x1= group_mean2, y1 = 100,lty=3,col='blue')
    
    se <- function(x) sqrt(var(x)/length(x))
    se1  <- se(processed_bam_files[[1]])
    se2 <- se(processed_bam_files[[2]])
    
    legend('topright', bty = "n", legend = c ('D=', paste(round(group_mean1-group_mean2,2)), '±', paste(round( se1+se2,2))))
    
    
    file <- 'REPAT.txt'
    if(name!= FALSE){
      write_matrix <- matrix(c(paste(name),paste(name),'Group1', 'Group2', as.numeric(group_mean1), as.numeric(group_mean2)),nrow=2,ncol=3) 
    }
    else{
      write_matrix <- matrix(c(paste(peak),paste(peak),'Group1', 'Group2', as.numeric(group_mean1), as.numeric(group_mean2)),nrow=2,ncol=3)   
    }
    

    
    write.table(write_matrix,file,append=TRUE,row.names=FALSE, col.names=FALSE, quote=FALSE) 
    
    
    
  }
  
  
  if(legend == TRUE){
    legend("topright",bty = "n", legend = paste(bam_list), fill = rainbow (16) ,text.width=40)
  }
  
}
# Gets the poly A tails from a bam file and returns them as a list

gff_gene_finder <- function(gff, name) {
  ### Outpus the peaks matching the input gene name
  
  gff_peaks <-gff
  
  # Search through for gene of interest.
  index1 <- with(gff_peaks, grepl (paste('=',name,';',sep=""), gff_peaks[,9],fixed=TRUE))
  output <-gff_peaks[index1, ]
 
  return (output)
}
  

poly_A_puller<- function(bam_file, gff, name, peak){
  ### Pulls all reads identifying with a gene/peak of interst 
  
  # Run peak finder or name finder depending on input 
  goi<- gff_gene_finder(gff, name)
  
  # Sets the column names of the gff rows data frame  
  colnames(goi) <- c('chr', 'program', 'type', 'peak start', 'peak end','DNS','orientation', 'DNS2','information')
  #peakno <- goi[,'information'] 
  #print(peakno)
  

  # Remove duplicate peaks
  goi <- subset(goi, !duplicated(goi[,'peak start'])) 

  
  # Split the gene of interest by orientation
  split_goi <- split(goi, goi[,'orientation'],drop = TRUE)
  
  # Pull the poly A reads for strands in each oreintation from the gff file
  plus_reads <- split_goi[['+']]
  minus_reads <- split_goi[['-']]
  
  
  list_of_peaks <-  list()
  # For plus and minus reads, assign reads to peaks 
  all_poly_a_tails_minus <- numeric()
  last_line <- data.frame(matrix(1, nrow=1,ncol=10))
  if (length(minus_reads)>=1){
    for (line in 1:nrow(minus_reads)){
      peak <- minus_reads [line,]
      # If the start of a line is less than the end of a previous peak, move the start of the peak to the end of the previous peak. 
      if(peak[,4] <= last_line[,5] & peak[,5] >= last_line[,5]){
        peak[,4] <-  last_line[,5]       
      }
      
      param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(peak [,'chr'],IRanges(
        peak[,'peak start'], peak[,'peak end'] +5 )))
      
      result1 <- scanBam ( bam_file , param = param, isMinusStrand = TRUE)
      
      # The list of poly-A tails from the bam file
      no_of_as1 <- result1[[1]][[5]][[1]]
      
      reads_in_peak_neg <- length (no_of_as1)
      list_of_peaks <-  c(list_of_peaks, length(no_of_as1))
  
      
      # Reset last line to current peak
      
      
      # Add succesive peaks together
      all_poly_a_tails_minus <-c(all_poly_a_tails_minus,no_of_as1)
      last_line <- peak
      
    }
    
  }

  # The poly_A flag values 
  all_poly_a_tails_plus <- numeric()
  # Give last line a default value 
  last_line <- data.frame(matrix(1, nrow=1,ncol=10)) 
  # Last line is set to the last peak. If it overlaps another peak it will make its end the new peak's start. 
  if (length(plus_reads)>=1){
    for (line in 1:nrow(plus_reads)){
      
      peak <- plus_reads [line,]
      # If the start of a line in less than the end of a previous peak, move the start of the peak to the end of the previous peak. 
      if(peak[,4] <= last_line[,5] & peak[,5] >= last_line[,5]){
        peak[,4] <-  last_line[,5]       
      }
      
      
      param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(peak [,'chr'],IRanges(
        peak[,'peak start'] -305, peak[,'peak end']+305 )))
      
      result2 <- scanBam(bam_file , param=param, isMinusStrand = FALSE)
      
      # Calculates the length of the cigar and adds this to the pos of the positive bam reads to give us the 3' end
      result2[[1]][[3]]<-result2[[1]][[3]]+result2[[1]][[4]]-1
      
      result2 <- result2[[1]]
      
      if (length(result2[[5]][[1]])!= 0){
        
        df <- cbind(unlist(result2[[3]]), unlist(result2[[5]][[1]]))
        
   
        # Second set of precessig to account for 5'read ends at the ead of a peak
        if (nrow(df)>0){
          my.data.frame <-  subset(df, df[,1] >= (peak[,'peak start'] -5) & df[,1] <= (peak[,'peak end']+5))
          
          no_of_as2 <- my.data.frame[,2]
          list_of_peaks <-  c(list_of_peaks, length(no_of_as2))
          # Add succesive peaks together
          all_poly_a_tails_plus <-c(all_poly_a_tails_plus, no_of_as2)
          
          reads_in_peak_pos <- length (no_of_as2)
          
         
          
          
          
          last_line <- peak
          
        }       
        
      }}    
    
  } 
  # Combine the pulled reads back together
  all_poly_a_tails<- sort(c(all_poly_a_tails_plus, all_poly_a_tails_minus))
  
  return (all_poly_a_tails)
 }