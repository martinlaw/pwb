##### How to approach #####
#
# Run this file
#
# Run simon_sim.R and examine the output (plot, and tables if desired)
#
# Run GS_sim.R and examine the output (plot, and tables if desired)
#
# To see to see workings (if desired), look at simon_example.R for an example
# outside of R functions, and pwbSimon.R and pwbGS.R for actual functions.


install.packages("librarian")
librarian::shelf(ggplot2, clinfun, xtable, DHARMa, martinlaw/curtailment, mjg211/singlearm)

# functions for finding precision weighted bias and displaying results:
source("R/pwbSimon.R")
source("R/pwbGS.R")
source("R/showTable.R")

#### Guide to new files ####
#
# simon_example.R
# For a single Simon design, obtains bias for a single true response rate.
# Script is self-contained to more easily examine the code.
#
# pwbSimon.R
# Function. For a single Simon design, obtains bias for a single true response rate.
#
# simon_sim.R
# Uses the function pwbSimon() to obtain bias for a range of true response rates.
# Finds estimates for UMVUE.
# Plots the bias.
# Uses the function showTable() to show results (bias and SE) for single response rates
#
#
#
# pwbGS.R
# Function. For a single group sequential design, obtains bias for a single true response rate.
#
# GS_sim.R
# Uses the function pwbGS() to obtain bias for a range of true response rates.
# Plots the bias.
# Uses the function showTable() to show results (bias and SE) for single response rates.
#


