# This script allows to model the error risk in bat calls automatic identification made with 
# any software, by modelling the success probability of an identification as a function of the 
# confidence index provided by the software.Using the models, we can then determine the confidence
# index threshold necessary to have the desired maximum error risk in the dataset (from 50% to 10%).
# This allows to make an objective data selection for datasets too large to be fully checked. 
# Then we checked if the automatic identification success and false negatives are biased by the 
# environnemental context of the recording stations.
#############################################################################################
# By Kévin Barré, Julie Pauwels and Yves Bas
#############################################################################################

rm(list=ls())

library(boot) # version 1.3-20
library(lme4) # version 1.1-18-1
library(data.table) # version 1.11.4

setwd("./") # Your work directory

# csv file containing all the bat passes manually checked. It must contain at least those columns:
# - the site ID ("ID")
# - the species identified by the software ("simp_species")
# - the confidence index of the identification ("confidence_index")
# - the success/failure of automatic identification (0=error and 1=success) ("success")
data1 <- fread("./data_check.csv")

# csv file containing environmental variables to test the consistency of confidence indexes in relation to
# ennvironmental conditions
data2 <- fread("./environmental_variables.csv")

# Merging both files using the user ID column
datacheck <- as.data.frame(merge(data1, data2, by.x ="ID", by.y="ID"))

# Vector creation to store the models intercept (Int) and Estimate, the confidence index thresholds for each 
# maximum error risk tolerance tested (ER 50%, 40%, 30%, 20% and 10%),the false negative rates (FN) and
# false positive rates (FP),the estimate (EstimateVariable) and p-value (PVariable)for each environmental variable test, 
# the number of manual checks per species (Ncheck) and the proportion identification success per species (Success)
Int = vector() 
Estimate = vector() 
ER50 = vector() 
ER60 = vector() 
ER70 = vector()
ER80 = vector()
ER90 = vector() 
FNraw = vector()
FN50 = vector() 
FN60 = vector()
FN70 = vector() 
FN80 = vector() 
FN90 = vector() 
FPraw = vector()
FP50 = vector()
FP60 = vector() 
FP70 = vector() 
FP80 = vector() 
FP90 = vector() 
Pforest = vector()
Purban = vector()
Pwetland = vector()
Phedgerowlength = vector()
Pedge = vector()
EstimateForest = vector()
EstimateUrban = vector()
EstimateWetland = vector()
EstimateHedgerowlength = vector()
EstimateEdge = vector()
Ncheck = vector()
Success = vector()

# A logistic model is build for each species: success probability of automatic species identification as 
# a function of the confidence index provided by the software
for (i in unique(datacheck$simp_species))
{
  # Species selection
  data = datacheck[which(datacheck$simp_species == i),] 
  # Logistic modelling
  m <- glm(success ~ confidence_index, family=binomial(link = logit), data = data) 
  # Prediction scaling
  x <- seq(0, 1, 0.001) 
  # Logistic predictions
  y = inv.logit(m$coefficients[1] + m$coefficients[2]*x) 
  # Store of the intercept
  Int = c(Int, m$coefficients[1]) 
  # Store of the estimate
  Estimate = c(Estimate, m$coefficients[2])
  # Store of confidence indexes threshold corresponding to 
  # each maximum error risk tolerance (50%, 40%, 30%, 20% and 10%)
  ER50 = c(ER50, min(subset(x, y>0.5)))
  ER60 = c(ER60, min(subset(x, y>0.6)))
  ER70 = c(ER70, min(subset(x, y>0.7)))
  ER80 = c(ER80, min(subset(x, y>0.8)))
  ER90 = c(ER90, min(subset(x, y>0.9)))
  # Calculation of the false negative rate for each maximum error risk tolerance
  FNraw = c(FNraw, sum(subset(data$success, data$confidence_index<min(subset(x, y>0))))/sum(data$success))
  FN50 = c(FN50, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.5))))/sum(data$success))
  FN60 = c(FN60, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.6))))/sum(data$success))
  FN70 = c(FN70, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.7))))/sum(data$success))
  FN80 = c(FN80, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.8))))/sum(data$success))
  FN90 = c(FN90, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.9))))/sum(data$success))
  # Calculation of the false positive rate for each maximum error risk tolerance
  FPraw = c(FPraw, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0)))))
  FP50 = c(FP50, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.5))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.5)))))
  FP60 = c(FP60, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.6))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.6)))))
  FP70 = c(FP70, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.7))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.7)))))
  FP80 = c(FP80, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.8))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.8)))))
  FP90 = c(FP90, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.9))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.9)))))
  # Store the number of manual checks and identification success proportion
  Ncheck = c(Ncheck, length(data$success))
  Success = c(Success, sum(subset(data$success, data$success != 0))*100/length(data$success))
}

## Adding Myotis nattereri from the "species" column ("simp_species" only contains the "Myotis spp" level)
i = "Myonat"
# Species selection
data = datacheck[which(datacheck$species == i),] 
# Logistic modelling
m <- glm(success ~ confidence_index, family=binomial(link = logit), data = data) 
# Prediction scaling
x <- seq(0, 1, 0.001) 
# Logistic predictions
y = inv.logit(m$coefficients[1] + m$coefficients[2]*x) 
# Storage of the intercepts
Int = c(Int, m$coefficients[1]) 
# Storage of the estimates
Estimate = c(Estimate, m$coefficients[2])
# Storage of confidence indexes from the automatic identification corresponding to 
# each maximum error risk tolerance (50%, 40%, 30%, 20% and 10%)
ER50 = c(ER50, min(subset(x, y>0.5)))
ER60 = c(ER60, min(subset(x, y>0.6)))
ER70 = c(ER70, min(subset(x, y>0.7)))
ER80 = c(ER80, min(subset(x, y>0.8)))
ER90 = c(ER90, min(subset(x, y>0.9)))
# Calculation of the false negative rate generated by each maximum error risk tolerance
FNraw = c(FNraw, sum(subset(data$success, data$confidence_index<min(subset(x, y>0))))/sum(data$success))
FN50 = c(FN50, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.5))))/sum(data$success))
FN60 = c(FN60, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.6))))/sum(data$success))
FN70 = c(FN70, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.7))))/sum(data$success))
FN80 = c(FN80, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.8))))/sum(data$success))
FN90 = c(FN90, sum(subset(data$success, data$confidence_index<min(subset(x, y>0.9))))/sum(data$success))
# Calculation of the false positive rate for each maximum error risk tolerance
FPraw = c(FPraw, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0)))))
FP50 = c(FP50, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.5))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.5)))))
FP60 = c(FP60, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.6))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.6)))))
FP70 = c(FP70, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.7))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.7)))))
FP80 = c(FP80, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.8))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.8)))))
FP90 = c(FP90, sum(subset((1-data$success), data$confidence_index>=min(subset(x, y>0.9))))/length(subset(data$confidence_index, data$confidence_index>=min(subset(x, y>0.9)))))
# Storage of the number of manual checks
Ncheck = c(Ncheck, length(data$success))
Success = c(Success, sum(subset(data$success, data$success != 0))*100/length(data$success))


## Test of the dependancy of automatic identification success probability to the environmental context

for (i in unique(datacheck$simp_species))
{
  # Species selection
  data2 = datacheck[which(datacheck$simp_species == i),]
  # Run models on species for which there are successes and errors
  if ((sum(data2$success) != 0) & (sum(data2$success) != length(data2$success)))
  {
  # Models for each environmental variables
  mforest <- glmer(success ~ scale(distforest) + (1|ID), family=binomial(link = logit), data = data2)
  murban <- glmer(success ~ scale(disturban) + (1|ID), family=binomial(link = logit), data = data2)
  mwetland <- glmer(success ~ scale(distwetland) + (1|ID), family=binomial(link = logit), data = data2)
  mhedgerow <- glmer(success ~ scale(hedgerowlength) + (1|ID), family=binomial(link = logit), data = data2)
  medge <- glmer(success ~ edgetype + (1|ID), family=binomial(link = logit), data = data2)
  # p-values and estimates storage
  Pforest = c(Pforest, round(coef(summary(mforest))[2,'Pr(>|z|)'], 3))
  Purban = c(Purban, round(coef(summary(murban))[2,'Pr(>|z|)'], 3))
  Pwetland = c(Pwetland, round(coef(summary(mwetland))[2,'Pr(>|z|)'], 3))
  Phedgerowlength = c(Phedgerowlength, round(coef(summary(mhedgerow))[2,'Pr(>|z|)'], 3))
  Pedge = c(Pedge, round(coef(summary(medge))[2,'Pr(>|z|)'], 3))
  EstimateForest = c(EstimateForest, round(summary(mforest)$coefficients[2], 3))
  EstimateUrban = c(EstimateUrban, round(summary(murban)$coefficients[2], 3))
  EstimateWetland = c(EstimateWetland, round(summary(mwetland)$coefficients[2], 3))
  EstimateHedgerowlength = c(EstimateHedgerowlength, round(summary(mhedgerow)$coefficients[2], 3))
  EstimateEdge = c(EstimateEdge, round(summary(medge)$coefficients[2], 3))
  } 
  # If there are only identification successes or only identification errors, models can't be run and NAs are implemented 
   else {
    Pforest = c(Pforest, "NA")
    Purban = c(Purban, "NA")
    Pwetland = c(Pwetland, "NA")
    Phedgerowlength = c(Phedgerowlength, "NA")
    Pedge = c(Pedge, "NA")
    EstimateForest = c(EstimateForest, "NA")
    EstimateUrban = c(EstimateUrban, "NA")
    EstimateWetland = c(EstimateWetland, "NA")
    EstimateHedgerowlength = c(EstimateHedgerowlength, "NA")
    EstimateEdge = c(EstimateEdge, "NA")
  }
}

## Adding Myotis nattereri from the "species" column ("simp_species" only contains the "Myotis spp" level)
i = "Myonat"
data2 = datacheck[which(datacheck$species == i),]
# Models for each tested environmental variables
mforest <- glmer(success ~ scale(distforest) + (1|ID), family=binomial(link = logit), data = data2)
murban <- glmer(success ~ scale(disturban) + (1|ID), family=binomial(link = logit), data = data2)
mwetland <- glmer(success ~ scale(distwetland) + (1|ID), family=binomial(link = logit), data = data2)
mhedgerow <- glmer(success ~ scale(hedgerowlength) + (1|ID), family=binomial(link = logit), data = data2)
medge <- glmer(success ~ edgetype + (1|ID), family=binomial(link = logit), data = data2)
# p-values and estimates storage
Pforest = c(Pforest, round(coef(summary(mforest))[2,'Pr(>|z|)'], 3))
Purban = c(Purban, round(coef(summary(murban))[2,'Pr(>|z|)'], 3))
Pwetland = c(Pwetland, round(coef(summary(mwetland))[2,'Pr(>|z|)'], 3))
Phedgerowlength = c(Phedgerowlength, round(coef(summary(mhedgerow))[2,'Pr(>|z|)'], 3))
Pedge = c(Pedge, round(coef(summary(medge))[2,'Pr(>|z|)'], 3))
EstimateForest = c(EstimateForest, round(summary(mforest)$coefficients[2], 3))
EstimateUrban = c(EstimateUrban, round(summary(murban)$coefficients[2], 3))
EstimateWetland = c(EstimateWetland, round(summary(mwetland)$coefficients[2], 3))
EstimateHedgerowlength = c(EstimateHedgerowlength, round(summary(mhedgerow)$coefficients[2], 3))
EstimateEdge = c(EstimateEdge, round(summary(medge)$coefficients[2], 3))


# Placing results in a dataframe
ErrorRiskThresholds = as.data.frame(cbind(Species = c(simp_species = unique(datacheck$simp_species), "Myonat"), Int, Estimate
                                       , ER50, ER60, ER70, ER80, ER90, FNraw,
                                       FN50, FN60, FN70, FN80, FN90, FPraw,FP50, FP60, FP70, FP80, FP90, Ncheck, Success,
                                       EstimateForest, Pforest, EstimateUrban, Purban, EstimateWetland, 
                                       Pwetland, EstimateHedgerowlength, Phedgerowlength, EstimateEdge,Pedge))

# Saving the dataframe
write.csv(ErrorRiskThresholds, "./ErrorRiskThresholds.csv", row.names=F)

## Test of the dependency of false negatives to the environmental context
# Vector creations which will contain a binomial variable in relation with false negatives
# 1: false negative of the species
# 0: no false negative of the species
Barbar = seq(0,0,length.out=nrow(datacheck))
Eptser = seq(0,0,length.out=nrow(datacheck))
Myosp = seq(0,0,length.out=nrow(datacheck))
Nyclei = seq(0,0,length.out=nrow(datacheck))
Nycnoc = seq(0,0,length.out=nrow(datacheck))
Pipkuh = seq(0,0,length.out=nrow(datacheck))
Pipnat = seq(0,0,length.out=nrow(datacheck))
Pippip = seq(0,0,length.out=nrow(datacheck))
Plesp = seq(0,0,length.out=nrow(datacheck))
Rhihip = seq(0,0,length.out=nrow(datacheck))
# Vector creations to store esitmates and p-values
Pforest = vector()
Purban = vector()
Pwetland = vector()
Phedgerowlength = vector()
Pedge = vector()
EstimateForest = vector()
EstimateUrban = vector()
EstimateWetland = vector()
EstimateHedgerowlength = vector()
EstimateEdge = vector()
# Merging new columns with datacheck
data = cbind(datacheck, Barbar, Eptser, Myosp, Nyclei, Nycnoc, Pipkuh, Pipnat, Pippip, Plesp, Rhihip)
# Coding each species column :
# 1: false negative of the species
# 0: no false negative of the species
for (i in 17:26) {
  for (j in 1:nrow(data)) {
    data[j,i] <- ifelse((colnames(data)[i] == data$manual_check[j]) & (colnames(data)[i] != data$simp_species[j]),"1","0")
  }
}
# Subset to only select lines containing either a false negative or a true positive of the species
# Then modelling of the probability to have a false negative according to environmental variables
species = c("Barbar", "Eptser", "Myosp", "Nyclei", "Nycnoc", "Pipkuh", "Pipnat", "Pippip", "Plesp", "Rhihip")
for (i in species) {
  print(i) ; start = Sys.time()
  datas = subset(data, subset = (data$simp_species == i) | (data[i] == "1"))
  datas2 = subset(datas, datas$manual_check == i)
    # Run models on species for which there are successes and errors
  datas2[i] <- as.numeric(datas2[,i])
  if ((sum(datas2[i]) != 0) & (sum(datas2[i]) != length(datas2[i])))
  {
  # Models for each environmental variables
  mforest <- glmer(as.numeric(datas2[,i]) ~ scale(datas2$distforest) + (1|ID), family=binomial(link = logit), data = datas2)
  murban <- glmer(as.numeric(datas2[,i]) ~ scale(datas2$disturban) + (1|ID), family=binomial(link = logit), data = datas2)
  mwetland <- glmer(as.numeric(datas2[,i]) ~ scale(datas2$distwetland) + (1|ID), family=binomial(link = logit), data = datas2)
  mhedgerow <- glmer(as.numeric(datas2[,i]) ~ scale(datas2$hedgerowlength) + (1|ID), family=binomial(link = logit), data = datas2)
  medge <- glmer(as.numeric(datas2[,i]) ~ datas2$edgetype + (1|ID), family=binomial(link = logit), data = datas2)
  # p-values and estimates storage
  Pforest = c(Pforest, round(coef(summary(mforest))[2,'Pr(>|z|)'], 3))
  Purban = c(Purban, round(coef(summary(murban))[2,'Pr(>|z|)'], 3))
  Pwetland = c(Pwetland, round(coef(summary(mwetland))[2,'Pr(>|z|)'], 3))
  Phedgerowlength = c(Phedgerowlength, round(coef(summary(mhedgerow))[2,'Pr(>|z|)'], 3))
  Pedge = c(Pedge, round(coef(summary(medge))[2,'Pr(>|z|)'], 3))
  EstimateForest = c(EstimateForest, round(summary(mforest)$coefficients[2], 3))
  EstimateUrban = c(EstimateUrban, round(summary(murban)$coefficients[2], 3))
  EstimateWetland = c(EstimateWetland, round(summary(mwetland)$coefficients[2], 3))
  EstimateHedgerowlength = c(EstimateHedgerowlength, round(summary(mhedgerow)$coefficients[2], 3))
  EstimateEdge = c(EstimateEdge, round(summary(medge)$coefficients[2], 3))
  } 
    # If there are only identification successes or only identification errors, models can't be run and NAs are implemented 
    else {
      Pforest = c(Pforest, "NA")
      Purban = c(Purban, "NA")
      Pwetland = c(Pwetland, "NA")
      Phedgerowlength = c(Phedgerowlength, "NA")
      Pedge = c(Pedge, "NA")
      EstimateForest = c(EstimateForest, "NA")
      EstimateUrban = c(EstimateUrban, "NA")
      EstimateWetland = c(EstimateWetland, "NA")
      EstimateHedgerowlength = c(EstimateHedgerowlength, "NA")
      EstimateEdge = c(EstimateEdge, "NA")
  }
}

# Storing results in a dataframe
FalseNegativesCheck = as.data.frame(cbind(Species = c(species = unique(species)), 
                                          EstimateForest, Pforest, EstimateUrban, Purban, EstimateWetland, 
                                          Pwetland, EstimateHedgerowlength, Phedgerowlength, EstimateEdge,Pedge))

# Saving the dataframe
write.csv(FalseNegativesCheck, "./FalseNegativesCheck.csv", row.names=F)
