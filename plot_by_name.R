
### This R script takes a list of two bam files, a gff file and an NM number for a gene of interest and returns
### a survival curve of poly A tail distribution over the two samples based on peaks called in the gff file


survival_curve_generator_by_name <- function(bam_list, gff, name=FALSE, housekeeping_gene=FALSE,peak =FALSE){
        library(Rsamtools)   
        
        #  Read the gff file
        
        gff_file <- read.table (gff, sep = "\t", header = FALSE)
        
        
        print (c('Order of bam file processing', bam_list))
        
        
        # Pull the poly A reads from each bam file
        processed_bam_files <- lapply(bam_list , poly_A_puller, gff_file, name, peak)
        
        
        # Separate the pulled reads for plotting
        l1 <- processed_bam_files[[1]]
        l2 <- processed_bam_files[[2]]
        
        if (housekeeping_gene == TRUE) {
                difference <<- length(l1)/length(l2)
                print (toString(c(paste(bam_list[[2]]),'multiplied by',difference,'is equivalent to', paste(bam_list[[1]]))))
        }
        vodka_test <- ks.test(l1,l2, alternative="two.sided")
        print (vodka_test)
        
        
        # Make a reverse cumulative distribuition plot of the poly A reads
        r<- range(l1,l2)
        en <- ecdf(l1)
        eg <- ecdf(l2)
        
        
        ry <- max(length(l1),length(l2))
        if(peak!=FALSE){
        curve((1-en(x))*100, from=r[1], to=r[2], col="red", xlim=r, ylab= 'Percentage longer', xlab = 'Poly A length', main= c('peak',paste(peak)))}
        else{
                curve((1-en(x))*100,  from=r[1], to=r[2], col="red", xlim=r, ylab= 'Percentage longer', xlab = 'Poly A length', main= paste(name))}
        curve((1-eg(x))*100, from=r[1], to=r[2], col="blue", add=TRUE)
        
        legend("topright", legend = c(paste(bam_list[[1]]), paste(bam_list[[2]])), fill = c("red", "blue"),text.width=40)
        
        
}
poly_A_puller<- function(bam_file, gff, name, peak){
        
        gff_peak_finder <- function(gff_file, peak) {
                
                
                #Outpus the peaks matching a certain gene. NM for gene of interest must be in quotes (' ')
                
                
                gff_peaks <-gff
                
                numbers <- seq(1:length (gff_peaks [ ,1]))
                #Number the gff peaks
                numbered_gff <- cbind(gff_peaks, numbers)
                
                # Search through for gene of interest.
                
                
                
                found_gene <- subset(numbered_gff, numbered_gff[,10] == peak)
                
                output <-as.data.frame(found_gene)
                
                
                return (output)
                
        }
        
        
        gff_name_finder <- function(gff, name) {
                
                
                #Outpus the peaks matching a certain gene. NM for gene of interest must be in quotes (' ')
                
                
                gff_peaks <-gff
                
                numbers <- seq(1:length (gff_peaks [ ,1]))
                #Number the gff peaks
                numbered_gff <- cbind(gff_peaks, numbers)
                
                # Search through for gene of interest.
                
                numbered_gff[grep(name, numbered_gff[,9], ignore.case=TRUE), "name"] <- name
                
                found_gene <- subset(numbered_gff, numbered_gff["name"] == name)
                
                output <-as.data.frame(found_gene)
                
                
                return (output)
                
        }
                
             
        # Pull your gene of interest from the gff file. 
        
               
        if (peak != FALSE){
                goi<- gff_peak_finder(gff, peak)
                colnames(goi) <- c('chr', 'program', 'type', 'peak start', 'peak end','DNS','orientation', 'DNS2','information','numbers')
        }
        else{
                goi <- gff_name_finder(gff, name) 
                colnames(goi) <- c('chr', 'program', 'type', 'peak start', 'peak end','DNS','orientation', 'DNS2','information','numbers','name')
        }
        
        
        cat('peaks, used: \n')
        print(goi)
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
                        print (c('number of reverse strand reads per peak', line))
                        print (reads_in_peak_neg)
                        
                        #add succesive peaks together
                        all_poly_a_tails_minus<-c(all_poly_a_tails_minus,no_of_as1)
                        
                }
        }
        # The number of poly A values
        plus_read_count <- list()
        # The poly_A flag values 
        all_poly_a_tails_plus <- numeric()
        if (length(plus_reads)>=1){
                for (line in 1:nrow(plus_reads)){
                        
                        param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(plus_reads [,'chr'],IRanges(
                                plus_reads[line,'peak start'] -305, plus_reads[line,'peak end']+305 )))
                        
                        result2 <- scanBam(bam_file , param=param, isMinusStrand = FALSE)
                        
                        result2[[1]][[3]]<-result2[[1]][[3]]+result2[[1]][[4]]-1
                        df <- as.data.frame(result2)
                        
                        #second set of precessig to account for 5'read ends at the ead of a peak
                        my.data.frame <- subset(df , df[,3] >= plus_reads[line,'peak start'] -5| df[,3] <= plus_reads[line,'peak end']+5)
                        no_of_as2 <- my.data.frame[,5]
                                               
                        #add succesive peaks together
                        all_poly_a_tails_plus <-c(all_poly_a_tails_plus, no_of_as2)
                        
                        reads_in_peak_pos <- length (no_of_as2)
                        plus_read_count <- c(plus_read_count,reads_in_peak_pos)
                        print (c('number of forward strand reads per peak',line)) 
                        print (reads_in_peak_pos)
                        
                }
                
        }
        # Combine the pulled reads back together
        all_poly_a_tails<- sort(c(all_poly_a_tails_plus, all_poly_a_tails_minus))
        cat('next condition\n')
        
        
        return (all_poly_a_tails)
        
}