gff_chromosome_subsetter <- function(counts_file, gff_file, chr) {
        
        
        peak_counts <- read.table(counts_file, sep= ",",header=TRUE, comment.char="",fill=TRUE, quote="")
        gff_peaks <-read.delim(gff_file, header=FALSE, comment.char="")
     
        
        
        
        found_chr <- subset(gff_peaks , gff_peaks [1] == chr)
        
        found_chr <-as.data.frame(found_chr)
     
        col_add<- gsub(".*id=","",found_chr[,9])
        
        coladd_2 <- gsub(";mean.*","",col_add)
        new_gff <- cbind(found_chr, coladd_2)
   
 
        extracted_peaks <- subset(peak_counts, sapply(peak_counts[,1], as.character) %in% sapply (new_gff[,'coladd_2'], as.character))
      
    
        return (extracted_peaks)
        
        # Get colmeans
        
        
}
        
        
      
        
