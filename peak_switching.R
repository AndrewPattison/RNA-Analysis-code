# Gets the poly A tails from a bam file and returns them as a list
poly_A_peak_puller<- function(peak,bam_file, gff){
        library(Rsamtools)        
        
        
        #Outpus the peaks matching a certain gene
        
        
        gff_peaks <-gff
        
        numbers <- seq(1:length (gff_peaks [ ,1]))
        
        # Number the gff peaks
        numbered_gff <- cbind(gff_peaks, numbers)
        
        # Search through for gene of interest.
        
        found_gene <- subset(numbered_gff, numbered_gff[,10] == peak)
        
        output <-as.data.frame(found_gene)
        
        
        
        
        
        # Run peak finder or name finder depending on input 
        
        
        
        goi<- output
        colnames(goi) <- c('chr', 'program', 'type', 'peak start', 'peak end','DNS','orientation', 'DNS2','information','numbers')
        
        
        # Split the gene of interest by orientation
        split_goi <- split(goi, goi[,'orientation'],drop = TRUE)
        
        # Pull the poly A reads for strands in each oreintation from the gff file
        plus_reads <- split_goi[['+']]
        minus_reads <- split_goi[['-']]
        
        
        list_of_peaks <-  list()
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
                        list_of_peaks <-  c(list_of_peaks, length(no_of_as1))
                        
                        
                        
                        # Add succesive peaks together
                        all_poly_a_tails_minus<-c(all_poly_a_tails_minus,no_of_as1)
                        
                        # Print to output file for later identification of 3' UTR switching
                        
                        p <- c(line,reads_in_peak_neg)
                        #write.table(p, file ,append=TRUE, row.names = FALSE, col.names = FALSE, eol = "\t", quote = FALSE)
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
                                list_of_peaks <-  c(list_of_peaks, length(no_of_as2))
                                # Add succesive peaks together
                                all_poly_a_tails_plus <-c(all_poly_a_tails_plus, no_of_as2)
                                
                                reads_in_peak_pos <- length (no_of_as2)
                                
                                
                                
                                # Print to output file for later identification of 3' UTR switching
                                p <- c(line,reads_in_peak_pos)
                                #write.table(p, file ,append=TRUE, row.names = FALSE, col.names = FALSE, eol = "\t", quote = FALSE)
                                
                                
                                last_line <- peak
                                
                        }}
                
        }
        # Combine the pulled reads back together
        all_poly_a_tails<- sort(c(all_poly_a_tails_plus, all_poly_a_tails_minus))
        
        
        
        
        
        return (all_poly_a_tails)
        
}

peak_switching<- function(peak_list,bam_list, gff, number_of_replicates = 3, number_of_conditions =2){
        #This function currently only works for two conditions
        
        # Read in the gff file
        gff_file <- read.table (gff, sep = "\t", header = FALSE)
        
        # Set an empty read couts list for later use
        read_counts<- list()
        
        # Get differnt peaks for each BAM file and store them in read counts
        for (bam_file in bam_list){
                read_count_condition1 <- lapply (peak_list[[1]], poly_A_peak_puller, bam_file,gff_file)
                
                read_counts <- c(read_counts, read_count_condition1)
                
        }
        for (bam_file in bam_list){
                read_count_condition2 <- lapply (peak_list[[2]], poly_A_peak_puller, bam_file,gff_file)
                
                read_counts <- c(read_counts, read_count_condition2)
                
        }
        # Get the lengths of the poly A tail lists stored in read counts to get read counts
        length_list <- lapply (read_counts, length)
        # Put these values into a data frame
        length_frame <- do.call('rbind', length_list)
        # Add the list of file names at their corresponding counts 
        length_frame <- cbind(paste(bam_list), as.numeric(length_frame))
        
        # Order the data frame by name so that like conditions match up (may mave condition out of input order)
        length_frame <- length_frame[order(length_frame[,1], decreasing = FALSE),]
        
        
        # Get totals of peaks per conition
        totals <- colSums(matrix(as.numeric(length_frame[,2]), nrow=length(peak_list)))
        
        
        # Put totals in cols next to the vales that were added to reach them
        totals <- rep(totals, each = length(peak_list)) 
        
        # Calculate percentage peaks used
        length_frame <- cbind(length_frame, totals)
        
        length_frame[,2] <- as.numeric(length_frame[,2]) / (as.numeric(length_frame[,3]))*100
        
        # Split up your data frames into peak groups
        length_frame1 <- length_frame[,-3] [seq(1, nrow(length_frame), 2),]
        length_frame2 <- length_frame[,-3]  [seq(2, nrow(length_frame), 2),]
        
        stdev_lists1 <- split(matrix(as.numeric(length_frame1[,2])), ceiling(seq_along(length_frame1[,2])/number_of_replicates))
        print(stdev_lists1)
        stdev_lists1 <-lapply(stdev_lists1,sd)
        stdev_lists2 <- split(matrix(as.numeric(length_frame2[,2])), ceiling(seq_along(length_frame2[,2])/number_of_replicates))
        
        stdev_lists2 <-lapply(stdev_lists2,sd)
        print(stdev_lists2)
        
        
        sdevs<- as.numeric(c(stdev_lists1, stdev_lists2))
        sdevs <-cbind(paste(c(rep(bam_list[1],rep=2),rep(bam_list[(length(bam_list)/number_of_conditions)+1],rep=2))),sdevs)
        print(sdevs)
        sdevs <-sdevs [order(sdevs[,1], decreasing = TRUE),]
        print(sdevs)
        # Get the avearge proportion of these groups
        length_frame1 <- colMeans(matrix(as.numeric(length_frame1[,2]), nrow=length(bam_list)/number_of_conditions))
        length_frame2 <- colMeans(matrix(as.numeric(length_frame2[,2]), nrow=length(bam_list)/number_of_conditions))
        
        
        length_frame1 <- cbind(length_frame1)
        length_frame2 <- cbind(length_frame2)
        
        # Put them back together
        length_frame<- rbind(length_frame1,length_frame2)
        
        # Add names back in
        length_frame<-cbind(paste(c(rep(bam_list[1],rep=2),rep(bam_list[(length(bam_list)/number_of_conditions)+1],rep=2))),length_frame)
        
        # Order by sample
        length_frame <- length_frame[order(length_frame[,1], decreasing = TRUE),]
        
        avgs <- as.numeric(length_frame[,2]) 
        points<- c(1.5,4.5,7.5,10.5)
        
        
        barplot(main= 'peak switching', beside=TRUE, xlim=c(0,15), ylim = c(0,100), as.numeric(length_frame[,2]), xaxt= 'n' ,col=c('green','purple'),space =c(0,0), ylab = 'percentage peak used', axes=FALSE, xlab= '', width = 3:3, names= length_frame[,1])
        
        
        legend("topright", legend = paste('peak',peak_list), fill = c("green","purple"),text.width=5, bty = 'n')
        
        
        axis(2, pos=0, at= c(0,25,50,75,100), tick = 25)
        
        
        
        axis(1,hadj= 0.85,cex.axis=0.2,las=3, at=c(1.5,4.5,7.5,10.5), labels = length_frame[,1] )
        arrows(x0 = points, y0 = avgs, x1= points, y1 = avgs + as.numeric(sdevs[,2]), length=0.05, angle=90, code=2)
        
}
