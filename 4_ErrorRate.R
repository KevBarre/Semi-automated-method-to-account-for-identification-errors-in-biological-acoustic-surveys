# This script allows to compute real false positive rate for each species  and for each maximum error risk 
# tolerance using the error risk modelling performed in the script 1_ErrorRiskModelling.R.
#############################################################################################
# By Kévin Barré, Julie Pauwels and Yves Bas
#############################################################################################

rm(list = ls())

library(data.table) # version 1.11.4
library(boot) # version 1.3-20

setwd("./") # Your work directory

# Output of the automatic identification software
datatot <- fread("./data_total_speciesgroup.csv")

# Reading the file saved in the script 1 containing confidence indexes corresponding to each 
# maximum error risk tolerance (50%, 40%, 30%, 20%, 10%) 
ErrorRiskThresholds <- read.csv("./ErrorRiskThresholds.csv")

# Merging both files to get thresholds for each species in datatot
data = as.data.frame(merge(datatot, ErrorRiskThresholds[,c(1:8)],by.x="simp_species",by.y="Species"))
# Adding a column for the calcultation on raw data
ER0 = seq(0,0,length.out=nrow(datatot))
data = cbind(data, ER0)

tabErr = c()
Thresholds = c("ER0","ER50", "ER60", "ER70", "ER80", "ER90")

# Calculation of real false positive rates for each species at each thresholds (i.e. each maximum
# error risk tolerance)
for (i in unique(as.factor(data$simp_species))){
  
  # Subset on the species
  dataSp=subset(data, data$simp_species==i)
  ErrorRateTotal = c()
  
  for (j in Thresholds) {
    # Subset on the threshold
    dataT=subset(dataSp,dataSp$confidence_index>dataSp[,j])
    ErrorRate = 0
    Err = 0
    # Calculation of error probability corresponding to each confience index from the equation modelling
    # the error risk in the sript 1 on each rows (Err column) of the subset data
    Err = 1-inv.logit(dataT$Int+dataT$Estimate*dataT$confidence_index)
    # Averaging of all real rates
    ErrorRate = sum(Err)/nrow(dataT)
    ErrorRateTotal = c(ErrorRateTotal, round(ErrorRate,4))
    rm(ErrorRate, Err)
  }
  tabErr = rbind(tabErr, ErrorRateTotal)
}

# Adding Myotis nattereri from the "specie" column ("simp_species" only contains the "Myotis spp" level)
dataSp=subset(data, data$species=="Myonat")
ErrorRateTotal = c()

for (j in Thresholds) {
  
  dataT=subset(dataSp,dataSp$confidence_index>dataSp[,j])
  ErrorRate = 0
  Err = 0
  
  Err = 1-inv.logit(dataT$Int+dataT$Estimate*dataT$confidence_index)
  ErrorRate = sum(Err)/nrow(dataT)
  ErrorRateTotal = c(ErrorRateTotal, round(ErrorRate,4))
  rm(ErrorRate, Err)
}
tabErr = rbind(tabErr, ErrorRateTotal)

colnames(tabErr) = Thresholds
row.names(tabErr) = c(unique(as.character(data$simp_species)), "Myonat")
# Saving
write.csv(tabErr, "./ErrorRate_datatotal.csv", row.names=T)

