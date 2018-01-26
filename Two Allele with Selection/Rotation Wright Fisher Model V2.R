#Andrew Johnson
#Wright-Fisher model: two alleles, selection

#NEEDS-- MULTIPLE ALLELES
#        MULTINOMIAL DRAW 
#        MUTATION/ INTRODUCTION OF NEW ALLELES (THAT WILL NOT BE PART OF THE REPRODUCTIVE POOL)

#set initial conditions------------------------------------------------------------------------------------

#initial population number
popNum <- 1000

#allele frequencies
p <- .5
q <- (1 - p)

#allele fitness coefficient
sel <- 0.1

#generation counter, intitally set to zero
gens <- 0

#create data frame of appropriate dimensions--------------------------------------------------------------

genTable <- data.frame(matrix(ncol = 2, nrow = 1))
colnames(genTable)=c("p","q")

#binomial sample iterated over n generations until extinction or fix.-------------------------------------

while(p > 0 && p < 1){
  
  #adds one generation to the counter
  gens <- gens + 1
  
  #adding variable values from previous generation to the dataframe
  genTable <- rbind(genTable)
  genTable[gens,1]=p
  genTable[gens,2]=q
  
  #calculate the fitness-weighted probability
  pFit <- ((p^2 + ( p * q * ( 1 - (.5 * sel)))) / 
             (p^2 + (2 * p * q * (1 - (.5*sel))) + (q^2 * (1-sel))))
  
  #binomial sample of two variants, p and q in this case
  newDist <- rbinom(1, popNum, pFit)
  p <- newDist/popNum
  q <- 1 - p
}

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


