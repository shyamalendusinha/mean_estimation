
# Image 1

rm(list = ls())

setwd("~/mean_estimation")

library(VGAM)
library(invgamma)
library(rmutil)



# Image 2
library(MCMCpack)
library(invgamma)
library(truncnorm)
library(VGAM)
library(rootSolve)
library(quadprog)
rinvgamma = invgamma::rinvgamma
dinvgamma = invgamma::dinvgamma

load("~/mean_estimation/real_data_examples/dataprostate.RData")

controldata = prostatedata[,1:50]
treatmentdata = prostatedata[,51:102]

data = controldata

set.seed(1234)
data1 = data[,sample(dim(data)[2],3)]
#6 31 30 selected with that seed

xbar.all.vec = mud = apply(data,1,mean)
S2.all.vec = sigma2d = apply(data,1,var)

q = dim(data1)[1]
n = dim(data1)[2]
ids = rep(1:q,n)
x.vec = c(data1)
xbar.vec = tapply(x.vec,ids,mean)
S2.vec = tapply(x.vec,ids,var)
ni.vec = as.vector(table(ids))
S2.vec.xbar = S2.vec/ni.vec

n1 = 100
n2 = 100
mu.grid = seq(min(xbar.all.vec)-IQR(xbar.all.vec), max(xbar.all.vec)+IQR(xbar.all.vec),
              length.out=2*n1+1)[seq(2,2*n1,2)]
sigma2.grid = seq(max(min(S2.all.vec)-IQR(S2.all.vec),0.001), max(S2.all.vec)+IQR(S2.all.vec),
                  length.out=2*n2+1)[seq(2,2*n2,2)]

prec.grid = 1/sigma2.grid

source("~/mean_estimation/all_R_functions/DPMMfunction.R")
 
MCMCoutput = DPMM.nig(x.vec, ids, gamma.pi = 0.1, k = 10,
                       Burnin = 5000, Simsize = 5000, mu.grid = mu.grid,
                       sigma2.grid = sigma2.grid,
                       hyperparameters = c('Default'))

MCMCoutput_ng = DPMM.ng(x.vec, ids, gamma.pi = 0.1, k = 10,
                          Burnin = 5000, Simsize = 5000,
                         mu.grid = mu.grid, prec.grid = prec.grid,
                          hyperparameters = c('Default'))
 
 
save(list = ls(), file='prostatedataimage.RData')

load('prostatedataimage.RData')

pdf('prostateplot1.pdf',height = 10, width = 15)

par(mfrow=c(2,2))
mar.default = c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 

# plot(xbar.vec,S2.vec,cex=0.3,col='green',xlab='sample mean', ylab='sample var')
# points(xbar.all.vec,S2.all.vec,cex=0.3,col='red')
# legend()

plot(xbar.vec,S2.vec,cex=0.1, xlim = c(min(xbar.vec),max(xbar.vec)), ylim = c(min(S2.vec),
                                                                             max(S2.vec)),
     xlab = expression(bar(X)[.j]), ylab = expression(S[.j]^2) 
#, main = "sample mean vs. variance based on 3 columns"
)
plot(xbar.all.vec,S2.all.vec,cex=0.1, xlim = c(min(xbar.vec),max(xbar.vec)), ylim = c(min(S2.vec),
                                                                                     max(S2.vec)),
     xlab = expression(mu[j]), ylab = expression(sigma[j]^2)
     #, main = "sample mean vs. variance based on 50 columns"
     )

plot(MCMCoutput$mu.est,MCMCoutput$sigma2.est,cex=0.1, xlim = 
       c(min(xbar.vec),max(xbar.vec)), ylim = c(min(S2.vec),max(S2.vec)),
     xlab = expression(hat(mu)[j]), ylab = expression(hat(sigma[j])^2)
# , main = "estimated mean vs. variance"
)

par(mfrow=c(1,1))

dev.off()



pdf('prostateplot2.pdf',height = 10, width = 15)

par(mfrow=c(2,2))
mar.default = c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 

plot(mud,xbar.vec,cex =0.1,xlab=expression(mu[j]), ylab=expression(bar(X)[.j]),
     xlim=c(min(xbar.vec),max(xbar.vec)),ylim=c(min(xbar.vec),max(xbar.vec))
     # , main = "sample mean based on 50 columns vs. sample mean based on 3 columns"
     )
abline(a=0,b=1)

plot(mud,MCMCoutput$mu.est,cex =0.1,xlab=expression(mu[j]), ylab=expression(hat(mu)[j]),
     xlim=c(min(xbar.vec),max(xbar.vec)),ylim=c(min(xbar.vec),max(xbar.vec))
     # , main = "sample mean based on 50 columns vs. estimated mean based on 3 columns"
     )
abline(a=0,b=1)

plot(sigma2d,S2.vec,cex =0.1,xlab=expression(sigma[j]^2), ylab=expression(S[.j]^2),
     xlim=c(min(S2.vec),max(S2.vec)),ylim=c(min(S2.vec),max(S2.vec))
     # , main = "sample variance based on 50 columns vs. sample variance based on 3 columns"
     )
abline(a=0,b=1)

plot(sigma2d,MCMCoutput$sigma2.est,cex =0.1,xlab=expression(sigma[j]^2), ylab=expression(hat(sigma[j])^2),
     xlim=c(min(S2.vec),max(S2.vec)),ylim=c(min(S2.vec),max(S2.vec))
     #, main = "sample variance based on 50 columns vs. estimated variance based on 3 columns"
     )
abline(a=0,b=1)


par(mfrow=c(1,1))

dev.off()


d.mu.true = density(mud, n=length(mu.grid), from=mu.grid[1], to=mu.grid[length(mu.grid)])$y

d.mu.sample.est = density(xbar.vec, n=length(mu.grid), from=mu.grid[1], to=mu.grid[length(mu.grid)])$y

d.sigma2.true = density(c(sigma2d,-sigma2d), n=length(sigma2.grid), from=sigma2.grid[1], to=sigma2.grid[length(sigma2.grid)])$y

d.sigma2.sample.est = density(c(S2.vec,-S2.vec), n=length(sigma2.grid), from=sigma2.grid[1], to=sigma2.grid[length(sigma2.grid)])$y


pdf('prostateplot4.pdf',height = 5, width = 10)

par(mfrow=c(1,2))

plot(density(xbar.all.vec), lty=2,xlab = expression(mu[j]), 
     ylab = 'Density' #, main='True density of mean based on 50 columns'
     )
plot(density(c(S2.all.vec,-S2.all.vec)), xlim = c(0.1,3.5),lty=2,xlab = expression(sigma[j]^2), 
     ylab = 'Density' #, main='True density of variance based on 50 columns'
     )

par(mfrow=c(1,1))

dev.off()


pdf("examplesknown.pdf",height = 15, width = 10)


par(mfrow = c(4,2))
for(i in 1:8){
  filename = paste("diffq_sigma2known",i,".csv",sep = "")
  exam1 = read.csv(filename)
  serieslen = dim(exam1)[2]-3
  allother = vector()
  for (t in 1:serieslen){
    allother = c(allother,paste("value.",t,sep=""))
  }
  data1 = t(as.matrix(exam1[c(1,7,9,10,12,18),c("value","value.1","value.2","value.3",
                                              "value.4","value.5","value.6","value.7","value.8","value.9","value.10",
                                              "value.11","value.12")]))
  
  matplot(data1[,1], data1[,-1], lty=1:5, type = "b", pch=1:5, cex=1, col=1:5, xlab = 'q',
          ylab = 'MSE',main= paste('Example',i,sep=" "))
  legend("topright", legend = c("oracle.XKB","SURE.M.XKB","SURE.SG.XKB",
                                "GL.WMBZ","NIG mixture") , col=1:5, lty=1:5, pch=1:5, cex = 1)
}

par(mfrow = c(1,1))

dev.off()



pdf("examplesunknown.pdf",height = 15, width = 10)


par(mfrow = c(3,2))
for(i in 9:14){
  filename = paste("diffq",i,".csv",sep = "")
  exam1 = read.csv(filename)
  serieslen = dim(exam1)[2]-3
  allother = vector()
  for (t in 1:serieslen){
    allother = c(allother,paste("value.",t,sep=""))
  }
  data1 = t(as.matrix(exam1[c(1,7,9,10,12,18,23),c("value",allother)]))
  
  matplot(data1[,1], data1[,-1], lty=1:6, type = "b", pch=1:6, cex=1, col=1:6, xlab = 'q',
          ylab = 'MSE',main= paste('Example',i,sep=" "))
  legend("topright", legend = c("oracle.XKB","SURE.M.XKB","SURE.SG.XKB",
                                "GL.WMBZ","SURE.M.Double","NIG mixture") , col=1:6, lty=1:6, pch=1:6, cex = 1)
}

par(mfrow = c(1,1))

dev.off()


pdf("s-examplesunknown.pdf",height = 15, width = 10)


par(mfrow = c(3,2))
for(i in 9:14){
  filename = paste("diffq",i,".csv",sep = "")
  exam1 = read.csv(filename)
  serieslen = dim(exam1)[2]-3
  allother = vector()
  for (t in 1:serieslen){
    allother = c(allother,paste("value.",t,sep=""))
  }
  data1 = t(as.matrix(exam1[c(1,25,30),c("value",allother)]))
  
  matplot(data1[,1], data1[,-1], lty=1:5, type = "b", pch=1:5, cex=1, col=1:5, xlab = 'q',
          ylab = 'MSE',main= paste('Example',i,sep=" "))
  legend("topright", legend = c("SURE.M.Double",
                                "NIG mixture") , col=1:5, lty=1:5, pch=1:5, cex = 1)
}

par(mfrow = c(1,1))

dev.off()


boxdata = read.csv('diffq11-box-sort.csv')

pdf("sensitivity_gamma3.pdf",height = 10, width = 10)

par(mfrow = c(2,2))
boxplot(boxdata$prob.2_ng1~boxdata$q, 
        main = expression(paste(gamma,'=',0.1)),ylab=expression(pi[(9)]), xlab='q')
boxplot(boxdata$prob.2_ng2~boxdata$q,
        main = expression(paste(gamma,'=',10)), ylab=expression(pi[(9)]), xlab='q')
boxplot(boxdata$prob.2_ng3~boxdata$q,
        main = expression(paste(gamma,'=',50)), ylab=expression(pi[(9)]), xlab='q')
boxplot(boxdata$prob.2_ng4~boxdata$q,
        main = expression(paste(gamma,'=',100)), ylab=expression(pi[(9)]), xlab='q')
par(mfrow = c(1,1))
dev.off()

pdf("sensitivity_gamma4.pdf",height = 10, width = 10)

par(mfrow = c(2,2))
boxplot(boxdata$prob.3_ng1~boxdata$q, 
        main = expression(paste(gamma,'=',0.1)),ylab=expression(pi[(8)]), xlab='q')
boxplot(boxdata$prob.3_ng2~boxdata$q,
        main = expression(paste(gamma,'=',10)), ylab=expression(pi[(8)]), xlab='q')
boxplot(boxdata$prob.3_ng3~boxdata$q,
        main = expression(paste(gamma,'=',50)), ylab=expression(pi[(8)]), xlab='q')
boxplot(boxdata$prob.3_ng4~boxdata$q,
        main = expression(paste(gamma,'=',100)), ylab=expression(pi[(8)]), xlab='q')
par(mfrow = c(1,1))
dev.off()

pdf("sensitivity_gamma5.pdf",height = 10, width = 10)
par(mfrow = c(2,2))
boxplot(boxdata$avg.loss.mu.est.DPMM_ng1~boxdata$q,
        main = expression(paste(gamma,'=',0.1)), ylab='MSE', xlab='q')
boxplot(boxdata$avg.loss.mu.est.DPMM_ng2~boxdata$q,
        main = expression(paste(gamma,'=',10)), ylab='MSE', xlab='q')
boxplot(boxdata$avg.loss.mu.est.DPMM_ng3~boxdata$q,
        main = expression(paste(gamma,'=',50)), ylab='MSE', xlab='q')
boxplot(boxdata$avg.loss.mu.est.DPMM_ng4~boxdata$q,main = expression(paste(gamma,'=',100)),
        ylab='MSE', xlab='q')
par(mfrow = c(1,1))
dev.off()

pdf("sensitivity_gamma6.pdf",height = 10, width = 10)
par(mfrow = c(2,2))
boxplot(boxdata$avg.loss.sigma2d.est.DPMM_ng1~boxdata$q,
        main = expression(paste(gamma,'=',0.1)), ylab='MSE', xlab='q', ylim=c(0,3))
boxplot(boxdata$avg.loss.sigma2d.est.DPMM_ng2~boxdata$q,
        main = expression(paste(gamma,'=',10)), ylab='MSE', xlab='q', ylim=c(0,3))
boxplot(boxdata$avg.loss.sigma2d.est.DPMM_ng3~boxdata$q,
        main = expression(paste(gamma,'=',50)), ylab='MSE', xlab='q', ylim=c(0,3))
boxplot(boxdata$avg.loss.sigma2d.est.DPMM_ng4~boxdata$q,main = expression(paste(gamma,'=',100)),
        ylab='MSE', xlab='q', ylim=c(0,3))
par(mfrow = c(1,1))
dev.off()
