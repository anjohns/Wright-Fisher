#Andrew Johnson
#Wright-Fisher model: n alleles, selection

#NEEDS-- MULTIPLE ALLELES
#        MULTINOMIAL DRAW 
#        MUTATION/ INTRODUCTION OF NEW ALLELES (THAT WILL NOT BE PART OF THE REPRODUCTIVE POOL)

#set initial conditions------------------------------------------------------------------------------------

#initial population number
popNum <- 1000

#set the number of alleles in the system
numAlleles <- 5
alleleVariants <- c(1:numAlleles)

#creates n number allele variables, all initally equal in frequency
alleleName <- c()
freqVector <- c()
for(i in alleleVariants){
  nam <- paste("p", i, sep = "")
  assign(nam, 1/numAlleles)
  alleleName[i] <- nam
  freqVector[i] <- 1/numAlleles
  #designate an allelic fitness coefficient-- these are equally assigned, but they don't have to be
  alFit <- paste(nam, "Sel", sep = "")
  assign(alFit, 1)
}

#generation counter, intitally set to zero
gens <- 0

#create data frame of appropriate dimensions--------------------------------------------------------------

genTable <- data.frame(matrix(nrow = 0, ncol = numAlleles))
colnames(genTable) <- alleleName

#binomial sample iterated over n generations until extinction or fix.-------------------------------------
while(!(1 %in% freqVector)){
  
  #adds one generation to the counter
  gens <- gens + 1
  
  #adding variable values from previous generation to the dataframe
  genTable <- rbind(genTable, freqVector)

  
  #calculate the fitness-weighted probability
  pFit <- ((p^2 + ( p * q * ( 1 - (.5 * sel)))) / 
             (p^2 + (2 * p * q * (1 - (.5*sel))) + (q^2 * (1-sel))))
  
  #binomial sample of two variants, p and q in this case
  newDist <- rbinom(1, popNum, pFit)
  p <- newDist/popNum
  q <- 1 - p
}

colnames(genTable) <- alleleName

#plot the results-----------------------------------------------------------------------------------------
plot(genTable$p, col="goldenrod2", type = "l", xlab="Generations", ylab="", 
     ylim= c(0,1), main = " :) ")
lines(genTable$q, col="cadetblue3", type = "l")
par(bg= "white")
legend("right",legend = colnames(genTable),col=c("goldenrod2", "cadetblue3"),bg="white",lwd=3, cex = .8)

#prints generations until fixation or extinction----------------------------------------------------------
pFix <- 0
pExt <- 0
if(p > q){
  pFix <- 1
  sprintf("after %g generations, p fixed at 1", gens)
}else{
  pExt <- 1
  sprintf("after %g generations, p went extinct", gens)
}

lapply(freqVector, )


