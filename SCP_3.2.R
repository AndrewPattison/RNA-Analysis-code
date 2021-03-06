
### This R script takes a list of two bam files, a gff file and a gene name or 
### peak number for a gene of interest and returns a survival curve of poly A 
### tail distribution over the samples (up to six at a time) based on peaks 
### called in the gff file
### !to use!:
### SCP(gff, name, unequal_groups = FALSE, number_of_replicates = 3, combine = TRUE, two_curve = FALSE, save = FALSE, select = FALSE, plot_mean = TRUE, plot_legend = TRUE)
### gff <- gff file
### name <- name of gene of interest
### number of replicates defaults to 2
### if you have an uneven number of replicates between groups, set unequal_groups to TRUE and you will be asked to select groups.
### combine defaults to True which combines replicates into a single curve
### two curve defaults to False, when True produces two curve which divides 
### data to halves
### if select = T, bam files automatically extracted from working directory
### if "WD", wild type only, if "TD", non-wild type only, if "M", manually 
### select at console
### save = T, save by default name, save = F, dont save, or define save to the 
### desired save name
### plot_mean, plot_legend, logical parameter

#------------------------------------------------------------------------------
# function group of SCP
# user interacting function, main function, produces plot


create_bam_list <- function(select, bam_listWAR){
    all_files <- list.files(getwd(), pattern = "bam")
    bad_bam <- list.files(getwd(), pattern = "bai")
    bam_list <- setdiff(all_files, bad_bam)
    select_list <- list()
    return (bam_listWAR)
    if (select == FALSE){
        return(bam_list)
    }
    else if (select == "WT"){
        select_list <- bam_list[grep("WT", bam_list)]
    }
    else if (select == "TD"){
        select_list <- setdiff(bam_list, bam_list[grep("WT", bam_list)])
    }
    else if (select == "M"){
        select_list <- select.list(bam_list, multiple = T)
    }
    else{
      return (print('please select either "M" "WT" or "TD"'))
    }
    return(select_list)
}


make_plot_list <- function(processed_bam_files, number_of_replicates, unequal_groups){
    conditions_list <- list()
    if (unequal_groups == TRUE) {
      options <- seq(1,length(processed_bam_files))
      conditions_list <- select.list(options, multiple = TRUE, title= "test")
      
    }
    else{
      conditions_list <- tapply(processed_bam_files, (seq_along(processed_bam_files)-1) %/% number_of_replicates, list)
    }
      new_list <- list()
    for (reps in conditions_list){
        new_list <- c (new_list, list(unlist(reps, recursive = F)))
    }
    new_list <- lapply(new_list, sort)
    #print(conditions_list)
    return(new_list)
}


plot_curves <- function(new_list, colour_list, name){
    # Make a reverse cumulative distribuition plot of the poly A reads
    elis <- lapply (new_list, ecdf)
    ry <- lapply(new_list,max)
    r <- range(0,150)
    
    # Make graph have no box
    par(bty="l")
    par(mar=c(5.1,4.1,4.1,8.1), xpd =T)
    
    # Statement to determine title of plot, plots first curve. 
    # Plots all the remaining curves
    curve((1-elis[[1]](x))*100, from=r[1], to=r[2], col="white" , xlim=r, ylab= '', xlab = '', main= paste(name),axes=FALSE)   
    # Add a few more axis options
    axis(1, pos=0, at= c(0,50,100,150), tick = 50)
    axis(2, pos=0, at= c(0,25,50,75,100), tick = 25)
    # Plots all the remaining curves
    peak_plot <- function(elis, colour_list,r=r ){
        curve((1-elis(x))*100,  from=r[1], to=r[2], col=colour_list, xlim=r, ylab= '', xlab = '', add=TRUE)   
    }
    # Calls the peak plot function with colour parameters
    mapply (peak_plot, elis, colour_list, MoreArgs= list(r=r))
    return()
}


plot_mean_legend <- function(bam_list, number_of_replicates, processed_bam_files, legend_list, colour_list, plot_mean, plot_legend){
    mean_division <- tapply(processed_bam_files, (seq_along(processed_bam_files)-1) %/% (length(bam_list)/2), list)
    print(str(mean_division))
    if(length(bam_list) > 1 & number_of_replicates != length(bam_list) & plot_mean == TRUE){
        group_mean1<- mean(unlist(mean_division[[1]]))
        groupsd1 <- sd(unlist(mean_division[[1]]))
        list1 <- unlist(mean_division[[1]])
        list2 <- unlist(mean_division[[2]])
        se <- function(x) sqrt(var(x)/length(x))
        se1  <- se(list1)
        se2 <- se(list2)
        group_mean2<- mean(unlist(mean_division[[2]]))
        groupsd2 <- sd(unlist(mean_division[[2]]))
        segments(x0 =group_mean1 , y0 = 0, x1= group_mean1, y1 = 100,lty=3,col='black')
        segments(x0 =group_mean2 , y0 = 0, x1= group_mean2, y1 = 100,lty=3,col='red')
        segments(x0 =group_mean1 , y0 = 100, x1= group_mean2, y1 = 100)
        legend(bty = "n", x= min(group_mean1, group_mean2)-5,y = 110, legend = paste('D=', paste(round(group_mean2-group_mean1,2)), '±', paste(round( se1+se2,2))))   
    }
    if (plot_legend == TRUE){
        legend("topright",bty = "n",inset= c(-0.1,-0.02), legend = legend_list, fill = colour_list ,text.width=40)
    }
    
    # Place axis labels in correct postion
    text(75,-14, labels = 'Poly A Tail Length')
    text(-14,50, labels = 'Percent Population (%)', srt =90)
    return()
}

set_colour <- function(bam_list, number_of_replicates){
    if(length(bam_list)/number_of_replicates <= 8){
        default <- palette("default")
        colour_list <- default[1:(length(bam_list)/number_of_replicates)]     
    }
    else{
        colour_list <- rainbow(length(bam_list)/number_of_replicates)
    }
    return(colour_list)
}


plot_save <- function(name, save){
    if (save == TRUE){
        dir.create(file.path("./sCURVES/"), showWarnings= FALSE)
        setEPS(width= 10)
        file <- paste('./sCURVES/', name,'.eps', sep = "")
        postscript(file)
    }
    else if(save == FALSE){
    }
    else {
        dir.create(file.path("./sCURVES/"), showWarnings= FALSE)
        setEPS(width= 10)
        file <- paste('./sCURVES/', save,'.eps', sep = "")
        postscript(file)
    }
}


SCP<- function(bam_listWAR, gff, name, unequal_groups = FALSE, number_of_replicates = 3, combine = TRUE, two_curve = FALSE, save = FALSE, select = FALSE, plot_mean = TRUE, plot_legend = TRUE){
    library(Rsamtools)
    # Write plot output to file
    plot_save(name, save)
    bam_list <- create_bam_list(select, bam_listWAR)
    if(length(bam_list) %% number_of_replicates != 0){
        print("Bam files does not match number of replicates, the program will continue, however please check if the data is correct.")
    }
    if(combine == FALSE){
        number_of_replicates = 1
    }
    if(two_curve == TRUE){
        number_of_replicates = length(bam_list)/2
    }
    
    legend_list <- bam_list[seq(1, length(bam_list), number_of_replicates)]
    colour_list <- set_colour(bam_list, number_of_replicates)
    gff_file <- read.delim(gff, header=FALSE, comment.char="")
    #Prints out the order files were processed (the order of the input list)
    print (c('Order of bam file processing', bam_list))
    # Pull the poly A reads from each bam file, calls poly_A_puller
    processed_bam_files <- lapply(bam_list, poly_A_puller, gff_file, name)
    plot_list <- make_plot_list(processed_bam_files, number_of_replicates, unequal_groups)
    plot_curves(plot_list, colour_list, name)
    plot_mean_legend(bam_list, number_of_replicates, processed_bam_files, legend_list, colour_list, plot_mean, plot_legend)
    if (save != FALSE){
        dev.off()
    }
}


#-------------------------------------------------------------------------------
# function group of poly_A_puller
# data interacting function, extract poly A length infomation from bam and gff

gff_gene_finder <- function(gff, name){
    ### Outpus the peaks matching the input gene name   
    gff_peaks <- gff
    
    index1 <- with(gff_peaks, grepl (ignore.case = T,paste('=',name,';',sep=""), gff_peaks[,9]))
    output <-gff_peaks[index1, ] 
    
    if (nrow(output)==0){
      index1 <- with(gff_peaks, grepl (ignore.case = T,paste('=',name,'/',sep=""), gff_peaks[,9]))
    }
    output <-gff_peaks[index1, ] 
    if (nrow(output)==0){
        index1 <- with(gff_peaks, grepl (ignore.case = T,paste('=',name,'$',sep=""), gff_peaks[,9]))
    }
    output <-gff_peaks[index1, ]
    print(output)
    return(output)
}


minus_pull <- function(minus_reads, bam_file){
    list_of_peaks <- list()
    all_poly_a_tails_minus <- numeric()
    last_line <- data.frame(matrix(1, nrow=1,ncol=10))
    if (length(minus_reads) >= 1){
        for (line in 1:nrow(minus_reads)){
            peak <- minus_reads [line,]
            # If the start of a line is less than the end of a previous peak, move the start of the peak to the end of the previous peak. 
            if(peak[,4] <= last_line[,5] & peak[,5] >= last_line[,5]){
                peak[,4] <-  last_line[,5]       
            }
            param <- ScanBamParam(what=c('qname','pos','qwidth','strand'),tag=c('AN'), which=GRanges(peak [,'chr'],IRanges(
                peak[,'peak start'], peak[,'peak end'] +5 )))
            result1 <- scanBam (bam_file , param = param, isMinusStrand = TRUE)
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
    return(all_poly_a_tails_minus)
}


plus_pull <- function(plus_reads, bam_file){
    list_of_peaks <- list()
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
            }    
        }
    } 
    return(all_poly_a_tails_plus)
}


split_geneofinterest <- function(gff, name){
    goi<- gff_gene_finder(gff, name)
    # Sets the column names of the gff rows data frame  
    colnames(goi) <- c('chr', 'program', 'type', 'peak start', 'peak end','DNS','orientation', 'DNS2','information')
    # Remove duplicate peaks
    goi <- subset(goi, !duplicated(goi[,'peak start'])) 
    # Split the gene of interest by orientation
    split_goi <- split(goi, goi[,'orientation'],drop = TRUE)
    return(split_goi)
}


poly_A_puller <- function(bam_file, gff, name){
    ### Pulls all reads identifying with a gene/peak of interst 
    split_goi <- split_geneofinterest(gff, name)
    # Pull the poly A reads for strands in each oreintation from the gff file
    plus_reads <- split_goi[['+']]
    minus_reads <- split_goi[['-']]
    list_of_peaks <- list()
    # For plus and minus reads, assign reads to peaks
    all_poly_a_tails_minus <- minus_pull(minus_reads, bam_file)
    all_poly_a_tails_plus <- plus_pull(plus_reads, bam_file)
    # Combine the pulled reads back together
    all_poly_a_tails <- sort(c(all_poly_a_tails_plus, all_poly_a_tails_minus))
    return(all_poly_a_tails)
}
