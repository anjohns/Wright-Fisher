#Andrew Johnson
#Wright-Fisher model: positive frequency dependent selection (alleles at smaller freq. have lower fitnesses)
#                     time decay built in

library(beepr)
library(data.table)

#parameterized curve like a parabola 
#run with theya higher
#compare the shape of the success curves between models


#set initial conditions------------------------------------------------------------------------------------

#initial parameters
popNum <- 1000
numAlleles <- 3
mutationRate <- .0005
generations <- 500
#set freqDepCoefficient to shape the positive freq fitness curve. Larger numbers will make the curve
#more linear, smaller fractions will make the curve rise rapidly to ~1 
freqDepCoefficient <- 1

alleleVariants <- c(1:numAlleles)
theta <- mutationRate*popNum
alleleName <- c()
censusVector <- c()
censusVectorNA <- c()
frequencyVector <- c()
fitnessVector <- c()
ageVector <- c()

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
  
  ageVector[i] <- 0
  
}

#assigns fitnesses based on positive frequency dependency
fitnessVector <- 1 + (frequencyVector / (frequencyVector + freqDepCoefficient))

#create data frame of appropriate dimensions--------------------------------------------------------------
genTable <- data.frame(matrix(nrow = generations, ncol = numAlleles))

#multinomial sample iterated over n generations-----------------------------------------------------------

for(e in 1:generations){
  
  #adds one generation to the counter
  gens <- gens + 1
  ageVector <- ageVector + 1
  
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
      
      ageVector[index] <- 1
      
      #----updates census, frequency, and selection vector------
      #updates censusVector to include new mutants
      censusVector[index] <- 1
      
      #adds allele frequency to the frequency vector
      frequencyVector[index] <- censusVector[index]/popNum
      
      #creates a selection coefficient vector. here set randomly
      #the mutants are given a higher fitness here for interesting results
      fitnessVector[index] <- 1 + (frequencyVector[index] / (frequencyVector[index] + freqDepCoefficient))
      
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
  genTable[e,] <- frequencyVector
  
  #prevents frequencies from being lower than 1/N
  replace(frequencyVector, (frequencyVector < (1/popNum)), 0)
  
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
  
  fitnessVector <- (1 + (frequencyVector / (frequencyVector + freqDepCoefficient))) / ageVector
  fitnessVector <- replace(fitnessVector, is.na(fitnessVector), 0)
  
}

#names columns of genTable 
colnames(genTable) <- alleleName

#replace all 0's in gentable with null values
genTable[genTable == 0] <- NA

#plot the results-----------------------------------------------------------------------------------------

#creates a plot
matplot(y = genTable, type = 'l', lty = 1, xlab = "Generations", ylab = "Frequency", 
        main = "Positive Frequency Dependence with Age Decay")

#plot background color
par(bg= "white")

matplot(y= genTable$Allele70, type = 'l', lty = 1, xlab = "Generations", ylab = "Frequency", 
        main = "Pos Freq Dep Fitness Age Decay, zoomed", xlim = c(115,140))

#Signals that the code is finished
beep()

