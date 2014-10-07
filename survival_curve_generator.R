
### This R script takes a list of two bam files, a gff file and an NM number for a gene of interest and returns
### a survival curve of poly A tail distribution over the two samples based on peaks called in the gff file


survival_curve_generator <- function(bam_list, gff, gene_of_interest_NM){
        library(Rsamtools)
        
        #  Read the gff file
        
        gff_file <- read.table (gff , sep = "\t", header = FALSE)
        
        
        # Pull the poly A reads from each bam file
        processed_bam_files <- lapply(bam_list ,poly_A_puller, gff_file, gene_of_interest_NM)
       
        
        # Separate the pulled reads for plotting
        l1 <- processed_bam_files[[1]]
        l2 <- processed_bam_files[[2]]
        
        
        # Make a reverse cumulative distribuition plot of the poly A reads
        r<- range(l1,l2)
        en <- ecdf(l1)
        eg <- ecdf(l2)
        
        ry <- max(length(l1),length(l2))
        
        curve((1-en(x))*100, from=r[1], to=r[2], col="red", xlim=r, ylab= 'Percentage longer', xlab = 'Poly A length', main= paste(gene_of_interest_NM))
        curve((1-eg(x))*100, from=r[1], to=r[2], col="blue", add=TRUE)
        
        legend("topright", legend = c(paste(bam_list[[1]]), paste(bam_list[[2]])), fill = c("red", "blue"))
        
        
}
poly_A_puller<- function(bam_file,gff_file, gene_of_interst_NM){
        
        
        gff_gene_finder <- function(gff, gene_of_interest_NM) {
                
                
                # Outpus the peaks matching a certain gene. NM for gene of interest must be in quotes (' ')
                
                gene_of_interest <- gene_of_interest_NM
                gff_peaks <- gff_file
                
                numbers <- 1:length (gff_peaks [ ,1])
                # Number the gff peaks
                numbered_gff <- cbind(gff_peaks, numbers)
                
                # Search through for gene of interest.
                
                numbered_gff[grep(gene_of_interest, numbered_gff[,9]), "gene"] <- gene_of_interest
                
                found_gene <- subset(numbered_gff, numbered_gff["gene"] == gene_of_interest)
                
                output <-as.data.frame(found_gene)
                
                
                return (output)
                
        }     
        # Pull your gene of interest from the gff file. 

        goi <- gff_gene_finder(gff_file, gene_of_interst_NM)

        colnames(goi) <- c('chr', 'program', 'type', 'peak start', 'peak end','DNS','orientation', 'DNS2','information')


        # Split the gene of interest by orientation
        split_goi <- split(goi, goi[,'orientation'],drop = TRUE)

        # Pull the poly A reads for strands in each oreintation
        plus_reads <- split_goi[['+']]
        minus_reads <- split_goi[['-']]

        all_poly_a_tails_minus<- numeric()
        if (length(minus_reads)>1){
        for (line in 1:nrow(minus_reads)){
        
        param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(minus_reads [,'chr'],IRanges(
                minus_reads[,'peak start'], minus_reads[,'peak end'] )))
        
        result <- scanBam( bam_file , param=param, isMinusStrand = TRUE)
        no_of_as1 <- result[[1]][[5]][[1]]
        
        all_poly_a_tails_minus<-c(all_poly_a_tails_minus,no_of_as1)
        }
        }
        all_poly_a_tails_plus <- numeric()
        if (length(plus_reads)>1){
        for (line in 1:nrow(plus_reads)){
        
        param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(plus_reads [,'chr'],IRanges(
                plus_reads[,'peak start']-100, plus_reads[,'peak end'] )))
        
        result <- scanBam( bam_file , param=param, isMinusStrand = FALSE)
        no_of_as2 <- result[[1]][[5]][[1]]
        
        all_poly_a_tails_plus <-c(all_poly_a_tails_plus, no_of_as2)
        
        
}
}
# Combine the pulled reads back together
all_poly_a_tails<- sort(c(all_poly_a_tails_plus, all_poly_a_tails_minus))

return (all_poly_a_tails)

}