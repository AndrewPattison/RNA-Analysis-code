DPP <- function(sam_input, gff_file, gene_of_interest) {
        
        
        ###Outputs the reads that fall anywhere inside peaks matching a given gene. Gene selected by NM_ for gene of interest, which must be in quotes (' ')###
        
        #Reads in the gff file
        gff_peaks <- gff_file  #read.table (gff_file , sep = "\t", header = FALSE)
        
        #Generates a sequence of numbers to number the gff peaks
        
        numbers <- seq(1:length (gff_peaks [ ,1]))
        
        #Numbers the gff peaks
        numbered_gff <- cbind(gff_peaks, numbers)
        
        # Searches through gff file for specified gene of interest
        
        numbered_gff[grep(gene_of_interest, numbered_gff[,9]), "gene"] <- gene_of_interest
        
        found_gene <- subset(numbered_gff, numbered_gff["gene"] == gene_of_interest)
        
        #Outputs only peaks matching given gene
        output_gff <-as.data.frame(found_gene)
        
        
                
        #read in your .sam file
        
        sam <-sam_input    #read.table (sam_file, sep = "", header = FALSE, fill = TRUE)
        
        #Remove all reads that don't contain a poly A flag
        a_flag_sam_file <- subset( sam, sam [,21] != "")
                
        
        #bind an extra col of zeroes to sam file for later
        zero_sam_file<- cbind (a_flag_sam_file, c(0))
        
        #pull out reads
       
        last_peak <- data.frame(matrix(0, ncol =10, nrow = 1))
        
        for (I in 1:length(zero_sam_file[ ,1])) {
                
                read <- zero_sam_file[I,4]
                
                if (last_peak[,4] -100 <= read &  last_peak[,5] >= read){
                        zero_sam_file[I,24] <- last_peak[10]
                }
                                
                else{
                logv = (output_gff[,4] -100 <= read &  output_gff[,5] >= read)
                if (length(output_gff [logv,10]) ==1){
                        peak_no <- as.numeric(output_gff [logv,10])
                        last_peak <- output_gff [logv,] 
                        zero_sam_file[I,24] <- peak_no
                }  
        }
        }
        # A file with annotations in c24 where reads have been matched to peaks
        output <- zero_sam_file
        
        #Take from the output only lines assigned a peak 
        reads <- subset (output, output[ ,24]!= 0)
        
        #Pull the AN flag numbers
        
        out <- suppressWarnings(as.numeric(gsub("AN:i:", "", reads [,21])))
        out <- out[!is.na(out)]
        
        #Returns an ordered list of poly A residues for specified gene 
        return  (sort(out))
        
        
}