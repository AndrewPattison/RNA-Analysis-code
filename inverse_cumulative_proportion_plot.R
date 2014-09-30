#After sorting output from DPP.R l1 = Gld2 l2 = N2

l1 <- sort('Output condition 1')
l2 <- sort('Output condition 2')


# Make a reverse cumulative distribuition plot (*100 to give percentage(I think, could very easily be wrong))
r<- range(l1,l2)
en <- ecdf(l1)
eg <- ecdf(l2)

ry <- max(length(l1),length(l2))

curve((1-en(x))*100, from=r[1], to=r[2], col="red", xlim=r, ylab= 'Percentage longer', xlab = 'Poly A length', main= 'MOM-2')
curve((1-eg(x))*100, from=r[1], to=r[2], col="blue", add=TRUE)

legend("topright", legend = c("N2", "Gld2"), fill = c("red", "blue"))


