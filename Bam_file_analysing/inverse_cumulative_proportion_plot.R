#After sorting output from DPP.R l1 = Gld2 l2 = N2

        
en <- ecdf(l2)
eg<- ecdf(l1)
x <- range(l2,l1)
ry <- max(length(l2),length(l1))

curve((1-en(x))*100, from=r[1], to=r[2], col="red", xlim=r, ylab= 'Percentage longer', xlab = 'Poly A length', main= 'Pos1')
curve((1-eg(x))*100, from=r[1], to=r[2], col="blue", add=TRUE)

legend("topright", 
       legend = c("N2, "Gld2"), 
       fill = c("red"", "blue"))