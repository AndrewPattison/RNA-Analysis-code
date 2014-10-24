peak_difference_finder <- function (bam_list, gff_file, gene_of_interest_nm, difference = difference){
        
        #Outpus the peaks matching a certain gene. NM for gene of interest must be in quotes (' ')
        
        
        gff_peaks <- read.table (gff_file , sep = "\t", header = FALSE)
        
        numbers <- 1:length (gff_peaks [ ,1])
        #Number the gff peaks
        numbered_gff <- cbind(gff_peaks, numbers)
        
        # Search through for gene of interest.
        
        numbered_gff[grep(gene_of_interest, numbered_gff[,9]), "gene"] <- gene_of_interest
        
        found_gene <- subset(numbered_gff, numbered_gff["gene"] == gene_of_interest)
        
        goi <-as.data.frame(found_gene))


        # Assign a number of reads to each peak and multiply by difference factor for second experiment.

colnames(goi) <- c('chr', 'program', 'type', 'peak start', 'peak end','DNS','orientation', 'DNS2','information')


# Split the gene of interest by orientation
split_goi <- split(goi, goi[,'orientation'],drop = TRUE)

# Pull the poly A reads for strands in each oreintation from the gff file
plus_reads <- split_goi[['+']]
minus_reads <- split_goi[['-']]


minus_read_count <- list()
all_poly_a_tails_minus<- numeric()
if (length(minus_reads)>=1){
        for (line in 1:nrow(minus_reads)){
                
                param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(minus_reads [,'chr'],IRanges(
                        minus_reads[line,'peak start'], minus_reads[line,'peak end'] +5 )))
                
                result1 <- scanBam ( bam_file , param=param, isMinusStrand = TRUE)
                no_of_as1 <- result1[[1]][[5]][[1]]
                
                reads_in_peak_neg <- length (no_of_as1)
                minus_read_count <- c(minus_read_count, reads_in_peak_neg)
                #add succesive peaks together
                all_poly_a_tails_minus<-c(all_poly_a_tails_minus,no_of_as1)
        }
}

plus_read_count <- list()
all_poly_a_tails_plus <- numeric()
if (ncol(plus_reads)>=1){
        for (line in 1:nrow(plus_reads)){
                
                param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(plus_reads [,'chr'],IRanges(
                        plus_reads[line,'peak start'] -200, plus_reads[line,'peak end']+5 )))
                
                result2 <- scanBam(bam_file , param=param, isMinusStrand = FALSE)
                
                result2[[1]][[3]]<-result2[[1]][[3]]+result2[[1]][[4]]-1
                df <- as.data.frame(result2)
                
                my.data.frame <- subset(df , df[,3] >= plus_reads[line,'peak start'] -5| df[,3] <= plus_reads[line,'peak end']+5)
                no_of_as2 <- my.data.frame[,5]
                #add succesive peaks together
                all_poly_a_tails_plus <-c(all_poly_a_tails_plus, no_of_as2)
                reads_in_peak_pos <- length (no_of_as2)
                plus_read_count <- c(plus_read_count,reads_in_peak_pos)
                
                
        }




}