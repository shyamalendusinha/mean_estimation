rm(list = ls())

setwd("~/mean_estimation")

data = read.table('~/mean_estimation/real_data_examples/Brown_batting_data.txt', header=TRUE, sep=",", quote="")

data$Season1.AB = data$AB.4. + data$AB.5. + data$AB.6.
data$Season1.H = data$H.4. + data$H.5. + data$H.6.
data$Season2.AB = data$AB.7. + data$AB.8. + data$AB.9.10.
data$Season2.H = data$H.7. + data$H.8. + data$H.9.10.
data$Season1.X = asin(sqrt((data$Season1.H+0.25)/(data$Season1.AB+0.5)))
data$Season2.X = asin(sqrt((data$Season2.H+0.25)/(data$Season2.AB+0.5)))
data$Season1.S2 = 1/(4*data$Season1.AB)
data$Season2.S2 = 1/(4*data$Season2.AB)

datafull = data[which(data$Season1.AB > 10),]

datapicher = data1[which(data1$Pitcher.==1),]
datanonpicher = data1[which(data1$Pitcher.==0),]

source("~/mean_estimation/all_R_functions/DPMMfunction.R")

set.seed(1234)

MCMCoutput_ng = DPMM.ng.sknown(datafull$Season1.X, 1/datafull$Season1.S2)

datafull$mu.est = MCMCoutput_ng$mu.est
datafull1 = datafull[which(datafull$Season2.AB > 10),]

TSEP.DPMM_full = (sum((datafull1$mu.est-datafull1$Season2.X)^2-datafull1$Season2.S2))/
  (sum((datafull1$Season1.X-datafull1$Season2.X)^2-datafull1$Season2.S2))


MCMCoutput_ng_picher = DPMM.ng.sknown(datapicher$Season1.X, 1/datapicher$Season1.S2)

datapicher$mu.est = MCMCoutput_ng_picher$mu.est
datapicher1 = datapicher[which(datapicher$Season2.AB > 10),]

TSEP.DPMM_picher = (sum((datapicher1$mu.est-datapicher1$Season2.X)^2-datapicher1$Season2.S2))/
  (sum((datapicher1$Season1.X-datapicher1$Season2.X)^2-datapicher1$Season2.S2))


MCMCoutput_ng_nonpicher = DPMM.ng.sknown(datanonpicher$Season1.X, 1/datanonpicher$Season1.S2)

datanonpicher$mu.est = MCMCoutput_ng_nonpicher$mu.est
datanonpicher1 = datanonpicher[which(datanonpicher$Season2.AB > 10),]

TSEP.DPMM_nonpicher = (sum((datanonpicher1$mu.est-datanonpicher1$Season2.X)^2-datanonpicher1$Season2.S2))/
  (sum((datanonpicher1$Season1.X-datanonpicher1$Season2.X)^2-datanonpicher1$Season2.S2))


print(paste("TSEP.DPMM_full",TSEP.DPMM_full))
print(paste("TSEP.DPMM_picher",TSEP.DPMM_picher))
print(paste("TSEP.DPMM_nonpicher",TSEP.DPMM_nonpicher))
