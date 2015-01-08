
apply_DPP <- function(sam_list, gff_file, gene_of_interest){

gff_file <-  read.table (gff_file , sep = "\t", header = FALSE)


sam_files <- lapply(sam_list ,DPP, gff_file, gene_of_interest)


l1 <- sam_files[[1]]
l2 <- sam_files[[2]]


# Make a reverse cumulative distribuition plot (*100 to give percentage(I think, could very easily be wrong))
r<- range(l1,l2)
en <- ecdf(l1)
eg <- ecdf(l2)

ry <- max(length(l1),length(l2))

curve((1-en(x))*100, from=r[1], to=r[2], col="red", xlim=r, ylab= 'Percentage longer', xlab = 'Poly A length', main= paste(gene_of_interest))
curve((1-eg(x))*100, from=r[1], to=r[2], col="blue", add=TRUE)

legend("topright", legend = c(paste(sam_list[[1]]), paste(sam_list[[2]])), fill = c("red", "blue"))


}
DPP <- function(sam_input, gff_file, gene_of_interest) {
        
        
        ###Outputs the reads that fall anywhere inside peaks matching a given gene. Gene selected by NM_ for gene of interest, which must be in quotes (' ')###
        
        #Reads in the gff file
        gff_peaks <- gff_file  
        
        #Generates a sequence of numbers to number the gff peaks
        
        numbers <- 1:nrow(gff_peaks)
        
        #Numbers the gff peaks
        numbered_gff <- cbind(gff_peaks, numbers)
        colnames(numbered_gff) <- c('chr', 'program', 'type', 'read start', 'read end','DNS','orientation', 'DNS2','information','numbers')
        
        # Searches through gff file for specified gene of interest
        
        numbered_gff[grep(gene_of_interest, numbered_gff[,'information']), "gene"] <- gene_of_interest
        
        
        #Outputs only peaks matching given gene
        
        output_gff <- subset(numbered_gff, numbered_gff["gene"] == gene_of_interest)
        
       
        #read in your .sam file
        
        sam <- read.table (sam_input, sep = "", header = FALSE, fill = TRUE)
        
        #Remove all reads that don't contain a poly A flag
        a_flag_sam_file <- subset( sam, sam [,21] != "")
                
        
        #bind an extra col of zeroes to sam file for later
        zero_sam_file<- cbind (a_flag_sam_file, c(0))
        colnames (zero_sam_file) <- c('qname', 'flag', 'chr', 'position', 'mapping quality', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'quality', 'AS', 'XN', 'XM', 'XO', 'XG', 'NM', 'MD', 'YT', 'RG', 'AN', 'AD','AA', 'zeroes')
        
        #pull out reads
        #Last peak is to make code run faster if reads map sequentially 
        last_peak <- data.frame(matrix(0, ncol =10, nrow = 1))
        
        for (I in 1:nrow(zero_sam_file)) {
                
                read <- zero_sam_file[I,]
                ori <- bitwAND (read[I,2],16)
                if (ori ==0){
                        
                
                if (last_peak[,4] -100 <= read[,'position'] & last_peak[,5] >= read[,'position'] ){
                        zero_sam_file[I,24] <- last_peak[10]
                }
                                
                else{
                        
                logv = (output_gff[,'read start'] -100 <= read[,'position']  &  output_gff[,'read end'] >= read[,'position'] )
                
                if (logv == TRUE){
                        peak_no <- as.numeric(output_gff [logv,'numbers'])
                        last_peak <- output_gff [logv,] 
                       
                        zero_sam_file[I,'zeroes'] <- peak_no
                } } 
        }
        else {
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        }
        }
        # A file with annotations in c24 where reads have been matched to peaks
        output <- zero_sam_file
        
        #Take from the output only lines assigned a peak 
        reads <- subset (output, output[ ,'zeroes']!= 0)
        
        #Pull the AN flag numbers
        
        out <- suppressWarnings(as.numeric(gsub("AN:i:", "", reads [,'AN'])))
        out <- out[!is.na(out)]
        
        #Returns an ordered list of poly A residues for specified gene 
        return  (sort(out))
        
        
}