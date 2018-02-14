#Andrew Johnson
#Wright-Fisher model: negative frequency dependent selection (alleles at smaller freq. have higher fitnesses)

library(beepr)

#set initial conditions------------------------------------------------------------------------------------

#initial parameters
popNum <- 1000
numAlleles <- 3
mutationRate <- .00005
generations <- 10000

alleleVariants <- c(1:numAlleles)
theta <- mutationRate*popNum
alleleName <- c()
censusVector <- c()
censusVectorNA <- c()
frequencyVector <- c()
fitnessVector <- c()

#counter, intitally set to zero
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
  
  #creates a selection coefficient vector. in this model each successive allele will be more fit
  # than the last. default at 1  
  fitnessVector[i] <- 1 + 1/censusVector

}

#create data frame of appropriate dimensions--------------------------------------------------------------
genTable <- data.frame(matrix(nrow = 0, ncol = numAlleles))

#multinomial sample iterated over n generations-----------------------------------------------------------

for(e in 1:generations){
  
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
      fitnessVector[index] <- 1 + 1/censusVector[index]
      
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
  
  #prevents frequencies from being lower than 1/N
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
  censusVectorNA <- replace(censusVector, censusVector==0, NA)
  
  #censusVector append numNew 1's to end of censusVector
  frequencyVector <- censusVector/popNum
  
  fitnessVector <- 1 + 1/censusVectorNA
  fitnessVector <- replace(fitnessVector, is.na(fitnessVector), 0)
  
}

#names columns of genTable 
colnames(genTable) <- alleleName

#replace all 0's in gentable with null values
genTable[genTable == 0] <- NA

#plot the results-----------------------------------------------------------------------------------------

#creates a plot
matplot(y = genTable, type = 'l', lty = 1, xlab = "Generations", ylab = "Frequency", 
        main = "Negative Frequency Dependent Fitness")

#plot background color
par(bg= "white")


#Signals that the code is finished
beep()



