test_puller <- function(sam_file, gff_file, gene_of_interest) {
        
        
        #Outpus the peaks matching a certain gene. NM for gene of interest must be in quotes (' ')
        
        
        gff_peaks <- read.table (gff_file , sep = "\t", header = FALSE)
        
        numbers <- seq(1:length (gff_peaks [ ,1]))
        #Number the gff peaks
        numbered_gff <- cbind(gff_peaks, numbers)
        
        # Search through for gene of interest.
        
        numbered_gff[grep(gene_of_interest, numbered_gff[,9]), "gene"] <- gene_of_interest
        
        found_gene <- subset(numbered_gff, numbered_gff["gene"] == gene_of_interest)
        
        output_gff <-as.data.frame(found_gene)
        
        
                
        #read in your .sam file
        
        #sam_file <- read.table (sam_file, sep = "", header = FALSE, fill = TRUE)
        
                
        #remove all reads that don't have a ploy A flag
        
        a_flag_sam_file <- subset(sam_file, sam_file [,21]  != "")
        
        #bind an extra col of zeroes to sam file for later
        zero_sam_file<- cbind (a_flag_sam_file, c(0))
        
        #pull out reads
       
        last_peak <- data.frame(matrix(0, ncol =10, nrow = 1))
        
        for (I in 1:length(zero_sam_file[ ,1])) {
                
                read <- zero_sam_file[I,4]
                
                if (last_peak[,4] <= read &  last_peak[,5] >= read){
                        zero_sam_file[I,24] <- last_peak[10]
                }
                                
                else{
                logv = (output_gff[,4] <= read &  output_gff[,5] >= read)
                if (length(output_gff [logv,10]) ==1){
                        peak_no <- as.numeric(output_gff [logv,10])
                        last_peak <- output_gff [logv,] 
                        zero_sam_file[I,24] <- peak_no
                }  
        }
        }
        #Take from the output only lines assigned a peak. 
        reads <- subset (output, output[ ,24]!= 0)
        
        #Pull the AN flag numbers
        
        out <- as.numeric(gsub("AN:i:", "", reads [,21]))
        
        
        return  (reads)
        
        
}