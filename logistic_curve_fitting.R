# 25 Sept 2017
# code to fit a logistic curve for our dose/OD curve 
# author <christa.caggiano@ucsf.edu> 

# needs dose response curve package for use
# install.packages("drc")
library("drc")

# sets working directory to appropriate repo 
setwd("/Users/Christa.Caggiano/Documents/UCSF_year1/PUBS/") 
df = read.csv("doubling_time.csv")

dose = df$Dose # x is dose 
OD = df$OD..actual. # OD is response
OD = y/.13 # fa

mL <- drm(OD ~ dose, data = df, fct = L.3(), type = "continuous") # fits a 3 parameter continuous curve
mL # parameters for logisitic equation 

plot(mL, type="all", col=454, main="logistic curve for dose/OD data", panel.first=grid()) # plots curve 

