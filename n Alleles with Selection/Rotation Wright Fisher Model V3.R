#Andrew Johnson
#Wright-Fisher model: n alleles, selection

#NEEDS-- MUTATION/ INTRODUCTION OF NEW ALLELES (THAT WILL NOT BE PART OF THE REPRODUCTIVE POOL)

library(ggplot2)

#set initial conditions------------------------------------------------------------------------------------

#initial population number
popNum <- 1000

#set the number of alleles in the system
#This number should be >= 1
numAlleles <- 5
alleleVariants <- c(1:numAlleles)

#throws if allele numer is smaller than 1
if(numAlleles < 1){
  stop("There has to be at least one allele for this model to function")
}else{

#creates n number allele variables, all initally equal in frequency
alleleName <- c()
frequencyVector <- c()
fitnessVector <- c()
for(i in alleleVariants){
  #creates a variable name for each allele
  nam <- paste("Allele", i, sep = "")
  #begins the model with all allele frequencies being equal
  assign(nam, 1/numAlleles)
  #adds the variable name to the name vector (string variable)
  alleleName[i] <- nam
  #adds allele frequency to the frequency vector
  frequencyVector[i] <- 1/numAlleles
  #creates a selection coefficient vector. can be set as a random number. default at 1 (no selection) 
  fitnessVector[i] <- 1 + i/10000
}
censusVector <- frequencyVector * popNum

#generation counter, intitally set to zero
gens <- 0

#create data frame of appropriate dimensions--------------------------------------------------------------
genTable <- data.frame(matrix(nrow = 0, ncol = numAlleles))

#multinomial sample iterated over n generations until extinction or fix.----------------------------------
while(!(1 %in% frequencyVector)){

  #adds one generation to the counter
  gens <- gens + 1
  
  #adding variable values from previous generation to the dataframe
  genTable <- rbind(genTable, frequencyVector)

  #calculate the fitness-weighted probability
  probabilityVector <- (fitnessVector * frequencyVector)/ (sum(fitnessVector * frequencyVector))
  
  #multinomial sample to update the next generation distribution of alleles
  censusVector <- rmultinom(1, popNum, probabilityVector)
  censusVector <- as.vector(censusVector)
  frequencyVector <- censusVector/popNum

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
plot(genTable$Allele1, type = "l", xlab="Generations", ylab="", 
     ylim= c(0,1), xlim = NULL, col = listOfColors[1], main = " :) ")

#cycles through each column of genTable and adds it to the plot
for(i in alleleVariants){
  lines(genTable[i], col = listOfColors[i], type = "l")
}

#plot background color
par(bg= "white")

#creates plot legend
legend("right",legend = colnames(genTable), col = listOfColors, bg="snow",lwd=3, cex = .8)



#end-----------------------------------------------------------------------------------------------------
#this bracket references the allele number catch
}











