## Semi-automated method to account for identification errors in biological acoustic surveys

To respond to the increasing use of biodiversity Passive Acoustic Monitorng (PAM), methods 
for detecting sound events, extracting numerous features, and automatically identifying species 
have been developed. However, automatic identification can generate large rate of errors which 
could affect response of bats to environmental variables and pressures. These R scripts propose a 
cautious method to account for identification errors in acoustic surveys without fully check records, 
using recordings of bat passes.

We proposed from manual checks (~ 1 900 bat passes) of a representative sample of the outputs (~ 212 000) 
of Tadarida automatic identification software to model the identification success probability of species 
as a function of the provided confidence index. This first step step (Script 1) allows to predict the confidence
index threshold necessary to have the desired maximum error risk in the dataset (from 50% to 10%).
This allows to make an objective data selection for datasets too large to be fully checked. 
We then investigated the effect of setting different maximum error rate tolerance under which data should
be discarded (Script 2), by repeating a large-scale analysis of bat activity response to habitat variables, 
and checking for consistency in the results (Script 3). The last script (Script 4) allows to compute real false 
positive rate for each species  and for each maximum error risk tolerance using the error risk modelling performed 
in the script 1_ErrorRiskModelling.R.

## Scripts
1_ErrorRiskModelling

2_RawdataSorting

3_ThresholdsModelling_EnvironmentalVariables

4_ErrorRate

## Corresponding data
data_check.csv # manual checks subset

data_total.csv # raw dataset from Tadarida software outputs

environmental_variables.csv # environmental data for each sites included in data_total.csv
