

### This R script takes a list of two bam files, a gff file and a gene name or peak number for a gene of interest and returns
### a survival curve of poly A tail distribution over the samples (up to six at a time) based on peaks called in the gff file


SCP <- function(bam_list, gff, name=FALSE, housekeeping_gene=FALSE,peak =FALSE, normalisation_factor = FALSE){
        library(Rsamtools)
        
        # Creates a .txt file named after the searched gene name 
        if (name != FALSE){
                fileConn <- file(paste(name,'.txt',sep = ""))
                writeLines('', fileConn)
                close(fileConn)
        }
        #  Read the gff file
        
        gff_file <- read.table (gff, sep = "\t", header = FALSE)
        
        #Prints out the order files were processed (the order of the input list)
        print (c('Order of bam file processing', bam_list))
        
        # Pull the poly A reads from each bam file, calls poly_A_puller
        processed_bam_files <- lapply(bam_list , poly_A_puller, gff_file, name, peak)
        
        # Samples of files from which to obtain plotting range
        l1 <- processed_bam_files[[1]]
        l2 <- processed_bam_files[[2]]
        
        number_of_reads1 <- length(l1)
        number_of_reads2 <- length(l2)
        
        if (housekeeping_gene != FALSE){
                norm_factor <-number_of_reads1/number_of_reads2 
                cat ('normalisation factor(condition 2/ condition 1) =', norm_factor)
        }
        
        norm_reads_1 <-numeric()
        norm_reads_2 <-numeric()
        if (normalisation_factor != FALSE){
                number_of_reads2 <- number_of_reads2* normalisation_factor
                norm_reads_1 <- toString(c('normalised reads', paste(bam_list[[1]]),number_of_reads1))
                norm_reads_2 <- toString(c('normalised reads', paste(bam_list[[2]]),number_of_reads2))
                cat(norm_reads_1,'\n')
                cat(norm_reads_2)
        }
        
        vodka_test <- ks.test(l1,l2, alternative="two.sided")
        print(t.test(l1,l2))
        print (vodka_test)
        
        
        # Make a reverse cumulative distribuition plot of the poly A reads
        elis <- lapply (processed_bam_files, ecdf)
        the_rest <- elis[c(2:length(bam_list))]
        ry <- lapply(processed_bam_files,max)
        r <- range(l1,l2)
        
        # If statement to determine title of plot, plots first curve. 
        if(peak!=FALSE){
                curve((1-elis[[1]](x))*100, from=r[1], to=r[2], col="red", xlim=r, ylab= 'Percentage longer', xlab = 'Poly A length', main= c('peak',paste(peak)))
        }
        else{
                curve((1-elis[[1]](x))*100, from=r[1], to=r[2], col="red", xlim=r, ylab= 'Percentage longer', xlab = 'Poly A length', main= paste(name))   
                
        }
        
        # Plots all the remaining curves 
        peak_plot <- function(ecdf_value, collist,r=r ){
                
                curve((1-ecdf_value(x))*100,  from=r[1], to=r[2], col=collist, xlim=r, ylab= 'Percentage longer', xlab = 'Poly A length', add=TRUE)
                
                
        }
        
        # Calls the peak plot function with colour parameters
        collist <- list("blue","green","orange", "yellow","purple")
        collist <- collist[c(1:length(bam_list)-1)]
        mapply (peak_plot, the_rest, collist, MoreArgs= list(r=r))
        
        # Plot legend 
        if (normalisation_factor!= FALSE){
                legend("topright", legend = c(paste(norm_reads_1),paste(norm_reads_2)), fill = c("red","blue","green","orange", "yellow","purple"),text.width=120)
        }
        else{
                legend("topright", legend = paste(bam_list), fill = c("red","blue","green","orange", "yellow","purple"),text.width=55)
        }
        
}
# Gets the poly A tails from a bam file and returns them as a list
poly_A_puller<- function(bam_file, gff, name, peak){
        
        gff_peak_finder <- function(gff_file, peak) {
                       
                
                #Outpus the peaks matching a certain gene
                
                
                gff_peaks <-gff
                
                numbers <- seq(1:length (gff_peaks [ ,1]))
                
                # Number the gff peaks
                numbered_gff <- cbind(gff_peaks, numbers)
                
                # Search through for gene of interest.
                
                found_gene <- subset(numbered_gff, numbered_gff[,10] == peak)
                
                output <-as.data.frame(found_gene)
                
                
                return (output)
                
        }
        
        
        gff_name_finder <- function(gff, name) {
                
                
                # Outpus the peaks matching the input gene name
                
                
                gff_peaks <-gff
                
                numbers <- seq(1:length (gff_peaks [ ,1]))
                # Number the gff peaks
                numbered_gff <- cbind(gff_peaks, numbers)
                
                # Search through for gene of interest.
                
                numbered_gff[grep(name, numbered_gff[,9], ignore.case=TRUE), "name"] <- name
                
                found_gene <- subset(numbered_gff, numbered_gff["name"] == name)
                
                output <-as.data.frame(found_gene)
                
                return (output)
                
        }
        
        
        # Run peak finder or name finder depending on input 
        
        
        if (peak != FALSE){
                goi<- gff_peak_finder(gff, peak)
                colnames(goi) <- c('chr', 'program', 'type', 'peak start', 'peak end','DNS','orientation', 'DNS2','information','numbers')
        }
        else{
                goi <- gff_name_finder(gff, name) 
                colnames(goi) <- c('chr', 'program', 'type', 'peak start', 'peak end','DNS','orientation', 'DNS2','information','numbers','name')
        }
        
        # Print peaks used for reference
        cat('peaks used: \n')
        print(goi)
        
        # Split the gene of interest by orientation
        split_goi <- split(goi, goi[,'orientation'],drop = TRUE)
        
        # Pull the poly A reads for strands in each oreintation from the gff file
        plus_reads <- split_goi[['+']]
        minus_reads <- split_goi[['-']]
        file <- toString(paste(name,'.txt',sep = ""))
        write(c('\n',bam_file), file ,append=TRUE)
        
        # For plus and minus reads, assign reads to peaks 
        all_poly_a_tails_minus<- numeric()
        last_line <- as.data.frame(data.frame(matrix(1, nrow=1,ncol=10))) 
        if (length(minus_reads)>=1){
                for (line in 1:nrow(minus_reads)){
                        peak <- minus_reads [line,]
                        # If the start of a line in less than the end of a previous read, move the start to the end of the previous read. 
                        if(peak[,4] <= last_line[,5]){
                                peak[,4] <-  last_line[,5]       
                        }
                        
                        param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(peak [,'chr'],IRanges(
                                peak[,'peak start'], peak[,'peak end'] +5 )))
                        
                        result1 <- scanBam ( bam_file , param=param, isMinusStrand = TRUE)
                        no_of_as1 <- result1[[1]][[5]][[1]]
                        
                        reads_in_peak_neg <- length (no_of_as1)
                        
                        print (c('number of reverse strand reads per peak', line))
                        print (reads_in_peak_neg)
                        
                        # Add succesive peaks together
                        all_poly_a_tails_minus<-c(all_poly_a_tails_minus,no_of_as1)
                        
                        # Print to output file for later identification of 3' UTR switching
                        
                        p <- c('Reverse strand peak',line,'Reads:',reads_in_peak_neg,'\n')
                        write.table(p, file ,append=TRUE, row.names = FALSE, col.names = FALSE, eol = "\t", quote = FALSE)
                }
        }
        
        
        # The poly_A flag values 
        all_poly_a_tails_plus <- numeric()
        # Give last line a default value 
        last_line <- as.data.frame(data.frame(matrix(1, nrow=1,ncol=10))) 
        # Last line is set to the last peak. If it overlaps another peak it will make its end the new peak's start. 
        if (length(plus_reads)>=1){
                for (line in 1:nrow(plus_reads)){
                        
                        peak <- plus_reads [line,]
                        # If the start of a line in less than the end of a previous read, move the start to the end of the previous read. 
                        
                        if(peak[,4] <= last_line[,5]){
                                peak[,4] <-  last_line[,5]       
                        }
                        
                       
                        param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(peak [,'chr'],IRanges(
                                peak[,'peak start'] -305, peak[,'peak end']+305 )))
                        
                        result2 <- scanBam(bam_file , param=param, isMinusStrand = FALSE)
                        result2[[1]][[3]]<-result2[[1]][[3]]+result2[[1]][[4]]-1
                        df <- as.data.frame(result2)
                        
                        # Second set of precessig to account for 5'read ends at the ead of a peak
                        if (nrow(df)>0){
                        my.data.frame <-  subset(df, df[,3] >= (peak[,'peak start'] -5) & df[,3] <= (peak[,'peak end']+5))
                        
                        no_of_as2 <- my.data.frame[,5]
                        
                        # Add succesive peaks together
                        all_poly_a_tails_plus <-c(all_poly_a_tails_plus, no_of_as2)
                        
                        reads_in_peak_pos <- length (no_of_as2)
                       
                        print (c('number of forward strand reads per peak',line)) 
                        print (reads_in_peak_pos)
                        
                        # Print to output file for later identification of 3' UTR switching
                        p <- c('Forward strand peak',line,'Reads:',reads_in_peak_pos,'\n')
                        write.table(p, file ,append=TRUE, row.names = FALSE, col.names = FALSE, eol = "\t", quote = FALSE)
                        
                        
                        last_line <- peak
                        
                        
                }}
                
        }
        # Combine the pulled reads back together
        all_poly_a_tails<- sort(c(all_poly_a_tails_plus, all_poly_a_tails_minus))
        
        cat('next condition\n')
        
        
        return (all_poly_a_tails)
        
}