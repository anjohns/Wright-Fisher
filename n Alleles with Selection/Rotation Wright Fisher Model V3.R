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
frequencyVector <- c()
fitnessVector <- c()
for(i in alleleVariants){
  nam <- paste("p", i, sep = "")
  assign(nam, 1/numAlleles)
  alleleName[i] <- nam
  frequencyVector[i] <- 1/numAlleles
  #creates a selection coefficient vector. can be set as a random number. default at 1 (no selection) 
  fitnessVector[i] <- 1
}
censusVector <- frequencyVector * popNum

#generation counter, intitally set to zero
gens <- 0

#create data frame of appropriate dimensions--------------------------------------------------------------

genTable <- data.frame(matrix(nrow = 0, ncol = numAlleles))
colnames(genTable) <- alleleName

#binomial sample iterated over n generations until extinction or fix.-------------------------------------
while(!(1 %in% frequencyVector)){

  #adds one generation to the counter
  gens <- gens + 1
  
  #adding variable values from previous generation to the dataframe
  genTable <- rbind(genTable, frequencyVector)

  #calculate the fitness-weighted probability
  probabilityVector <- (fitnessVector * frequencyVector)/ (sum(fitnessVector * frequencyVector))
  
  #binomial sample of two variants, p and q in this case
  censusVector <- rmultinom(1, popNum, probabilityVector)
  censusVector <- as.vector(censusVector)
  frequencyVector <- censusVector/popNum

}

#temporary fix (rename columns)-- while loop rbind command erases the column titles: reason unknown
colnames(genTable) <- alleleName

#plot the results-----------------------------------------------------------------------------------------
plot(genTable$p1, type = "l", xlab="Generations", ylab="", 
     ylim= c(0,1), main = " :) ")
for(i in numAlleles){
  lines(genTable[i], type = "l")
}
lines(genTable$p2, col = "red", type = "l")
par(bg= "white")
legend("right",legend = colnames(genTable),col=c("goldenrod2", "cadetblue3"),bg="white",lwd=3, cex = .8)

#prints generations until fixation or extinction----------------------------------------------------------














