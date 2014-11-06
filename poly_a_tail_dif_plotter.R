poly_a_tail_diff_plotter<- function(condition){

PAT_seq <- read.table('PAT_seq.txt')
REPAT_seq <- read.table('REPAT.txt')

both <- cbind (PAT_seq, REPAT_seq)
print(both)
# Split the data frame from N2 into Gld-2
split <- split(both, both[,2])
print(split)
Gld_2<- split[[1]]
N2<- split[[2]]

VALS <- cbind(Gld_2[,3],Gld_2[,6], N2[,3], N2[,6])
colnames(VALS)<- c('Gld-2 PAT-seq', 'Gld-2 RE-pat','N2 PAT-seq', 'N2 RE-pat')
xaxis <- 1: nrow(VALS)

labs <- sort(c ('Pos-1', 'Mex-1', 'Neg-1', 'Act-3', 'Nhl-2', 'Gld-3 peak 1', 'Gld-3 peak 2', 'Glp-1', 'Apx-1', 'Zif-1', 'Mom-2', 
           'Nos-2', 'Cye-1', 'Tpxl-1', 'Idi-1', 'Cpg-1', 'Ergo-1', 'Nra-2', 'Zk813.3','Mex-3', 'Pie-1', 'Spn-4', 'Lip-1', 'Daz-1', 'Rme-2',
           'Puf-5', 'Xnd-1', 'Fem-3', 'Gpd-4', 'Oma-2', 'Egg-1', 'Pup-2'))

# Correct for an error in plottting order
labs <- c(labs [1:2], labs[4], labs [3], labs [5:32])

print(labs)

if (condition== 'Gld-2'){
        mat = as.matrix(rbind(VALS[,1], VALS[,2]))
        colnames(mat) <- labs
        print(mat)
        par(cex.axis=0.7)
        par(mgp = c(3, 0.3, 0))
        barplot(beside =TRUE, mat , main= paste(condition), xlab= 'Gene',ylab= 'Mean Ploy-A Tail Length', col=c('red','blue') , yaxt='n', ylim= c(0,50), las=2, space =c(0.5,4))
        legend("topright",bty = "n", legend = c('PAT-seq', 'RE-pat'), fill = c("red","blue"))
        par(cex.axis=1, mgp = c(3, 1, 0))
        axis(2, at = c(0,10,20, 30, 40,50), pos= c(0,50))
        
}
else{ 
        mat = as.matrix(rbind(VALS[,3], VALS[,4]))
        colnames(mat) <- labs
        print(mat)
        par(cex.axis=0.7)
        par(mgp = c(3, 0.3, 0))
        barplot(beside =TRUE, mat , main= paste(condition), xlab= 'Gene',ylab= 'Mean Ploy-A Tail Length', col=c('red','blue') , yaxt='n', ylim= c(0,50), las=2, space =c(0.5,4))
        #legend("topright",bty = "n", legend = c('PAT-seq', 'RE-pat'), fill = c("red","blue"))
        par(cex.axis=1, mgp = c(3, 1, 0))
        axis(2, at = c(0,10,20, 30, 40,50), pos= c(0,50))
}

}