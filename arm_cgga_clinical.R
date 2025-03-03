
# ARM for clinical brain tumor data

# import necessary libraries
library(viridis)
library(arules)
library(TSP)
library(data.table)
library(ggplot2)
library(dplyr)
library(devtools)
library(purrr)
library(tidyr)
library(tcltk)
library(arulesViz)

# read in data as transaction data
data <- read.transactions("cgga_clinical_for_arm.csv", #cleaned and transformed csv file
                          rm.duplicates = FALSE, # do not remove duplicates
                          format = "basket", 
                          sep=",")



# implement apriori to get rules
rules = arules::apriori(data, parameter = list(support=.05, 
                                                 confidence=.9, minlen=2))
#check that  data is imported in correct format
inspect(data) 

# plot of which items are most frequent
itemFrequencyPlot(data, topN=20, type="absolute")

# get top 15 rules,sorted by support
sorted_rules_support = sort(rules, by="support", decreasing=TRUE)
inspect(sorted_rules_support[1:15])
summary(sorted_rules_support)

# get top 15 rules, sorted by confidence
sorted_rules_confidence = sort(rules, by="confidence", decreasing=TRUE)
inspect(sorted_rules_confidence[1:15])
summary(sorted_rules_confidence)

# obtain top 15 rules, sorted by lift
sorted_rules_lift = sort(rules, by="lift", decreasing=TRUE)
inspect(sorted_rules_lift[1:15])
summary(sorted_rules_lift)

# select for rules with specific RHS
deceased_rules = apriori(data=data,parameter = list(supp=.005, conf=.01, minlen=2),
                         appearance = list(default="lhs", rhs="deceased"),
                         control=list(verbose=FALSE))
sorted_deceased_rules = sort(deceased_rules, decreasing=TRUE, by="confidence")
inspect(sorted_deceased_rules[1:30])

# select rules with specific LHS
recurrent_rules <- apriori(data=data,parameter = list(supp=.005, conf=.01, minlen=2),
                           appearance = list(default="rhs", lhs="recurrent tumor"),
                           control=list(verbose=FALSE))
recurrent_rules <- sort(recurrent_rules, decreasing=TRUE, by="support")
inspect(recurrent_rules[1:15])

#visualize how support, conf, and lift relate for top 15 rules, sorted by support
sub_sorted_rules <- head(sort(sorted_rules, by="support"),10)
plot(sub_sorted_rules)

#visualize relationship between rules interactively
plot(sub_sorted_rules, method="graph", engine="interactive")