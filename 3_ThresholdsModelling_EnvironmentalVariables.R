# This script allows to check the consistency of bat species response to environmental
# variables according to the choosen threshold to sort data and to reduce the error rate 
# (R scripts 1 & 2). It is a safe approach consisting in the repetition of analyses at  
# different false positive tolerances. When results are consistent between thresholds 
# (e.g. between 50% and 10% of maximum error risk), this allows to study species with low  
# occurrences for which a restrictive threshold (e.g. 10%) removes too much data and do not 
# allow for robust statistical analysis.
#############################################################################################
# By Kévin Barré, Julie Pauwels and Yves Bas
#############################################################################################

rm(list = ls())

library(data.table) # vesrion 1.11.4
library(lme4) # version 1.1-18-1
library(beepr) # version 1.3

setwd("./") # Your work directory

# Reading the file saved in the script 2: number of bat passes per species per site and per threshold
Table <- fread("./CrossTableFiltered.csv")

# csv file containing environmental variables
Envt <- fread("./environmental_variables.csv")

data = as.data.frame(merge(Table,Envt,by.x ="V1",by.y="ID")) # V1 = auto-generated 1st column name containing ID in dataset "Table"

# Scaling of environmental variables
sdistF<-scale(data$distforest)
sdistW<-scale(data$distwetland)
sdistU<-scale(data$disturban)
shedg<-scale(data$hedgerowlength)

# Test of possible correlation issues between environemental variables
type<-as.factor(data$edgetype)
kruskal.test(sdistF, edgetype)
kruskal.test(sdistW, edgetype) 
kruskal.test(sdistU, edgetype) 
kruskal.test(shedg, edgetype) 
cor(sdistF,sdistW)
cor(sdistF,sdistU)
cor(sdistW,sdistU)
cor(shedg,sdistU)
cor(shedg,sdistW)
cor(shedg,sdistF)

# Loading VIF (variance inflation factor) function for GLMM to check for collinearity issues
# For each variable the VIF value should be < 2 
vif.mer <- function (fit) {
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

# Modelling the response of bat species (number of bat passes per night) to environmental
# variables for each false positive tolerance in the automatic identification
tab = c()
for (i in 2:73) { # Columns of data containing bat activity per species per threshold
  print(names(data)[i]) ; start = Sys.time()
  f = paste0(names(data)[i], "~ edgetype + sdistF + sdistW + sdistU + shedg + (1|date)")
  try = tryCatch({ # In order to not stop the process if a model do not converge
    if (strsplit(names(data)[i], "_")[[1]][1] %in% c("Nycnoc", "Nyclei", "Rhifer")){
      m<-glmer(as.formula(f), family = "poisson", data = data) # for species with poisson error distribution
    } else {m<-glmer.nb(as.formula(f), data=data)} # for species with negative binomial error distribution
  }, error = function(e) {
    print(paste0("Error in model ", names(data)[i]))
  })
  if (is.character(try) == F) {
    # store variables estimates
    s = as.data.frame(summary(m)$coefficients)
    s = cbind(variable = rownames(s),s)
    # store VIF values
    row.names(s) = NULL
    s$vif = c(NA, as.vector(vif.mer(m)))
  } else {
    # if the model did not run, stores NA
    s = as.data.frame(matrix(ncol = 6, nrow = 6, data = NA))
    colnames(s) = colnames(tab)
  }
  tab = rbind(tab, s)
  print(Sys.time() - start)
}

beep(2)

rnames = c()
for (i in 2:73) rnames = c(rnames, rep(names(data)[i], nrow(s)))
tab = cbind(species_thresholds = rnames, tab)
# Saving
write.csv(tab, "./Results_VariablesModelling.csv", row.names=F)


# Calculation of occurrences and number of  passes per species and per threshold
sdata = as.data.frame(data[,c(2:73)]) # Subset on species/thresholds columns
OccurrenceTotal = c()
NbTotal = c()
for (i in 1:ncol(sdata))
{
  print(i)
  datasp = sdata[,i] 
  Occurrence = length(subset(datasp, datasp > 0))/length(datasp)
  OccurrenceTotal = c(OccurrenceTotal, round(Occurrence, 3))
  rm(Occurrence)
  Nb = sum(datasp)
  NbTotal = c(NbTotal, Nb)
  rm(Nb)
}
# Dataframe
OccNb = cbind(OccurrenceTotal, NbTotal)
OccNb = cbind(espece_seuil = colnames(sdata), OccNb)
# Saving
write.csv(OccNb, "./OccurrencesAbundances2.csv", row.names=F)
