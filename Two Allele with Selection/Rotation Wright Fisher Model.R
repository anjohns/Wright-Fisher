#Andrew Johnson
#Wright-Fisher model: two alleles, selection

#set initial conditions---------------------------------------------------------------

#generation number
generations <- 1000
genNum <- c(1:generations)

#population number
popNum <- 10000

#allele frequencies
p <- .5
q <- (1 - p)

#allele fitness coefficient
sel <- 0.001

#create data frame of appropriate dimensions------------------------------------------

genTable <- data.frame(matrix(ncol = 2, nrow = genNum))
colnames(genTable)=c("p","q")

#binomial sample iterated over genNum generations-------------------------------------

for(i in genNum){
  
  #adding variable values from previous generation to the dataframe
  genTable[i,1]=p
  genTable[i,2]=q

  #calculate the fitness-weighted probability
  pFit <- ((p^2 + ( p * q * ( 1 - (.5 * sel)))) / 
             (p^2 + (2 * p * q * (1 - (.5*sel))) + (q^2 * (1-sel))))
  
  #binomial sample of two variants, p and q in this case
  newDist <- rbinom(1, popNum, pFit)
  p <- newDist/popNum
  q <- 1 - p
}

#plot the table to see the allelic progressions---------------------------------------
plot(genTable$p, col="goldenrod2", type = "l", xlab="Generations", ylab="", 
     ylim= c(0,1), main = " :) ")
lines(genTable$q, col="cadetblue3", type = "l")
par(bg= "white")
legend("right",legend = colnames(genTable),col=c("goldenrod2", "cadetblue3"),bg="white",lwd=3, cex = .8)
