#Andrew Johnson
#Wright-Fisher model: n alleles, selection

#set initial conditions------------------------------------------------------------------------------------

#initial conditions
popNum <- 1000
numAlleles <- 6
mutationRate <- .00005


alleleVariants <- c(1:numAlleles)
theta <- mutationRate*popNum
alleleName <- c()
censusVector <- c()
frequencyVector <- c()
fitnessVector <- c()
#generation counter, intitally set to zero
gens <- 0
newNum <- 0

randomizerSet <- c(-1000:1000)

#splits the population close to evenly to establish initial allele distributions
alleleCensusOtherThanLast <- popNum %/% numAlleles
lastAlleleCensus <- popNum - (alleleCensusOtherThanLast * (numAlleles - 1))

#throws if allele numer is smaller than 1
if(numAlleles < 1){
  stop("There has to be at least one allele for this model to function")
}


#creates n number allele variables, all initally equal in frequency
for(i in alleleVariants){
  #creates a variable name for each allele
  nam <- paste("Allele", i, sep = "")
  
  #begins the model with all allele frequencies being close to equal
  if(i < numAlleles){
    assign(nam, (alleleCensusOtherThanLast/popNum))
    censusVector[i] <- alleleCensusOtherThanLast
    frequencyVector[i] <- alleleCensusOtherThanLast/popNum
  }else{
    assign(nam, (lastAlleleCensus/popNum))
    censusVector[i] <- lastAlleleCensus
    frequencyVector[i] <- lastAlleleCensus/popNum
  }
  
  #adds the variable name to the name vector (string variable)
  alleleName[i] <- nam
  
  #creates a selection coefficient vector. can be set as a random number. default at 1 (no selection) 
  fitnessVector[i] <- 1 + sample(randomizerSet, 1, replace = TRUE)/10000
}

#create data frame of appropriate dimensions--------------------------------------------------------------
genTable <- data.frame(matrix(nrow = 0, ncol = numAlleles))

#multinomial sample iterated over n generations until extinction or fix.----------------------------------

#THIS WHILE CONDITION IS NOT APPROPRIATE FOR MODELING MUTATION, AS THERE CAN BE RANDOM MUTATION FROM
#ALLELES THAT HAVE FIXED AT 1. NEEDS INTELLECTUAL REVISION 
while(!(1 %in% frequencyVector)){

  #adds one generation to the counter
  gens <- gens + 1
  
  #machinery for adding new mutant alleles to the model-------------
  if(newNum > 0){
    for(j in 1:newNum){
      index <- j + numAlleles
      #----makes new variables corresponding to new mutant alleles------------
      #creates a variable name for each allele
      nam <- paste("Allele", index, sep = "")
      
      #mutated allele begins as a single individual in the population
      assign(nam, 1/popNum)
      
      #adds the variable name to the name vector (string variable)
      alleleName[index] <- nam
      
      #----updates census, frequency, and selection vector------
      #updates censusVector to include new mutants
      censusVector[index] <- 1
      
      #adds allele frequency to the frequency vector
      frequencyVector[index] <- censusVector[index]/popNum
      
      #creates a selection coefficient vector. here set randomly
      #the mutants are given a higher fitness here for interesting results
      fitnessVector[index] <- 1 + sample(randomizerSet, 1, replace = TRUE)/10000
      
      #----updates number of variants to include the new mutant and adds a column for the new variable----
      #updates allele variants range
      alleleVariants[index] <- index
      
      #adds a column to genTable where 
      genTable <- cbind(genTable, frequencyVector[index])
      genTable[1:(gens - 1), index] <- NA
    }
    
    #adds the new alleles to the sum total of alleles in numAlleles
    numAlleles <- numAlleles + newNum
    
    #resets newNum for the next iteration
    newNum <- 0
  }
  
  #--machinery reshaping the population, frequency, and probability distributions--------
  
  #adding variable values from previous generation to the dataframe
  genTable <- rbind(genTable, frequencyVector)
  
  for(z in 1:length(numAlleles)){
    if(frequencyVector[z] < 1/popNum){
      frequencyVector[z] <- 0
    }
  }

  #calculate the fitness-weighted probability
  probabilityVector <- (fitnessVector * frequencyVector)/ (sum(fitnessVector * frequencyVector))
  
  #number of new alleles in the population due to mutuation
  #theta is the mutation rate and the number of alleles is a poisson draw with mean = theta
  newNum <- rpois(1, theta)
  
  #multinomial sample to update the next generation distribution of alleles
  censusVector <- rmultinom(1, (popNum - newNum), probabilityVector)
  censusVector <- as.vector(censusVector)
  
  #censusVector append numNew 1's to end of censusVector
  frequencyVector <- censusVector/popNum
  
  #if(any census number in censusvector is 0){
  #   set that census and freq index to NULL -- something like genTable[X, Y] <- NA
  #   This should be conserved for each subsequent generation
  #}

}

#names columns of genTable 
colnames(genTable) <- alleleName

#plot the results-----------------------------------------------------------------------------------------


#vector of colors-- fourteen colors listed here to handle up to fourteen alleles. should be expanded 
#if additional alleles need to be visualized
listOfColors <- c("cyan3", "darkgoldenrod2", "coral2", "deepskyblue3", "firebrick", "darkgreen", "lightpink2",
                  "seagreen", "tan3", "sienna3", "aquamarine2", "orangered3", "mediumseagreen",
                  "slateblue")

#creates a plot
plot(genTable$Allele1, type = "l", lwd = 1, xlab="Generations", ylab="", 
     ylim= c(0,1), xlim = NULL, col = listOfColors[1], main = " :) ")

#cycles through each column of genTable and adds it to the plot
for(k in alleleVariants){
  lines(genTable[k], col = listOfColors[k], type = "l", lwd = 1)
}

#plot background color
par(bg= "white")

#creates plot legend
legend("right",legend = colnames(genTable), col = listOfColors, bg="snow",lwd=3, cex = .8)