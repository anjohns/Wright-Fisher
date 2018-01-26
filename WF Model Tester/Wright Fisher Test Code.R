#MODEL TEST FOR FIXATION PROBABILITIES

totTest <- 1000
testNum <- c(1:totTest)

#allele fitness coefficient
sel <- .01

fixQuan <- 0

#make a table of each test generations----------------------------------------------------
#testTable <- data.frame(generations=seq(1,totTest))
testTable <- data.frame(matrix(ncol = 3, nrow = totTest))
colnames(testTable)=c("generations","pFix", "pExt")

for(j in testNum){
  pFix <- 0
  pExt <- 0
  #population number
  popNum <- 10000
  diploidPop <- 2 * popNum
  
  #allele frequencies
  p <- .5
  q <- (1 - p)
  
  #allele fitness coefficient
  sel <- 0.0001
  
  #generation counter, intitally set to zero
  gens <- 0
  
  while(p > 0 && p < 1){
    
    #adds one generation to the counter
    gens <- gens + 1
    
    #calculate the fitness-weighted probability
    pFit <- ((p^2 + ( p * q * ( 1 - (.5 * sel)))) / 
               (p^2 + (2 * p * q * (1 - (.5*sel))) + (q^2 * (1-sel))))
    
    #binomial sample of two variants, p and q in this case
    newDist <- rbinom(1, diploidPop, pFit)
    p <- newDist/diploidPop
    q <- 1 - p
  }
  if(p > q){
    pFix <- 1
  }else{
    pExt <- 1
  }
  testTable[j, 1]= gens
  testTable[j, 2]= pFix
  testTable[j, 3]= pExt
}

#probability of fixation or extinction----------------------------------------------------------------------
fixProb <- (1- exp(-2 * popNum * sel * .5)) / (1- exp(-2 * popNum * sel))

#plot the test results----------------------------------------------------------------------------------------
hist(testTable$generations, xlab = "Generations", col = "cadetblue3", 
     main = "Generations to fixation or extinction")
abline(v = mean(testTable$generations), lwd = 2, col = "goldenrod2")

pFixAvg <- mean(testTable$pFix)
pFixAvg
fixProb













