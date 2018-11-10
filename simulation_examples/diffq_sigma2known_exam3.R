
rm(list = ls())

library(foreach)
library(doParallel)


#********Data Storing**************************

outputnumber = '3'

setwd("~/mean_estimation")
outputfile = paste("diffq_sigma2known",outputnumber,".Rdata",sep="")
outputfilecsv = paste("diffq_sigma2known",outputnumber,".csv",sep="")

#************Data Simulation parameters*********************
Bq = 1000
q.upper.lim = 500
q.lower.lim = 20
q.interval = 40
B = Bq*((q.upper.lim-q.lower.lim)/q.interval+1)

q.vec = vector()
for (b in 1:B){
  q.vec[b] = ceiling(b/Bq)*q.interval+(q.lower.lim-q.interval)
}

#****Data simulation parameters******************
sigma2min = 0.1
sigma2max = 1

r.biv.density = function(q){
  sigma2d = runif(q, min=sigma2min, max=sigma2max)
  mud = sigma2d
  return(list(mud=mud, sigma2d=sigma2d))
}

source("~/mean_estimation/all_R_functions/DPMMfunction.R")
source("~/mean_estimation/all_R_functions/SUREmethods.R")
source("~/mean_estimation/all_R_functions/grouplinearfunction_all.R")

cl = min(40,B)
registerDoParallel(cl)

alloutput_diffq = foreach(b = 1:B,.combine='cbind',.inorder=FALSE) %dopar% {
  
  set.seed(1234+47*b)
  q = q.vec[b]
  n1 = 100
  n2 = 100
  random.sample = r.biv.density(q)
  mud = random.sample$mud
  sigma2d = random.sample$sigma2d
  errord = rnorm(q)
  x.vec = mud+sqrt(sigma2d)*errord
  
  mu.grid = seq(min(x.vec)-IQR(x.vec), max(x.vec)+IQR(x.vec), 
                length.out=2*n1+1)[seq(2,2*n1,2)]
  sigma2.grid = seq(max(min(sigma2d)-IQR(sigma2d),0.001), max(sigma2d)+IQR(sigma2d), 
                    length.out=2*n2+1)[seq(2,2*n2,2)]
  
  prec.grid = 1/sigma2.grid
  
  
  #naive density estimation
  
  avg.loss.mu.est.sample.mean = mean((mud-x.vec)^2)
  
  #********EBMLE***********
  
  mu.est.EBMLE = thetahat.EBMLE.XKB(x.vec, sigma2d)
  
  avg.loss.mu.est.EBMLE = mean((mud-mu.est.EBMLE)^2)
  
  #**********************
  
  #*******EBMOM**********
  
  mu.est.EBMOM = thetahat.EBMOM.XKB(x.vec, sigma2d)
  
  avg.loss.mu.est.EBMOM = mean((mud-mu.est.EBMOM)^2)
  
  #*****************************
  
  #*************JS estimate*****
  
  mu.est.JS = thetahat.JS.XKB(x.vec, sigma2d)
  
  avg.loss.mu.est.JS = mean((mud-mu.est.JS)^2)
  
 
  #****************************
  
  #********Oracle**************
  
  mu.est.oracle = thetahat.oracle.XKB(x.vec, sigma2d, mud)
  
  avg.loss.mu.est.oracle = mean((mud-mu.est.oracle)^2)
  

  #************************************
  
  #********large class oracle**************
  
  m.or = mean(mud)
  or.loss = vector()
  
  for (i in 1:q){
    if(mud[i] < x.vec[i] & mud[i] > m.or){
      or.loss[i] = 0
    }
    else if(mud[i] > x.vec[i] & mud[i] < m.or){
      or.loss[i] = 0
    } else{
      or.loss[i] = min((x.vec-mud)^2,(mud-m.or)^2)
    }
  }
  avg.loss.mu.est.true.oracle = mean(or.loss)
  
  
  #*****SURE.G method*****************
  
  mu.est.SURE.G = thetahat.SURE.G.XKB(x.vec, sigma2d)
  
  avg.loss.mu.est.SURE.G = mean((mud-mu.est.SURE.G)^2)
  
 
  #************************************
  
  #*****SURE.M method*****************
  
  mu.est.SURE.M = thetahat.SURE.M.XKB(x.vec, sigma2d)
  
  avg.loss.mu.est.SURE.M = mean((mud-mu.est.SURE.M)^2)
  
  
  #************************************
  
  #*****SURE.SG  method*****************
  
  mu.est.SURE.SG = thetahat.SURE.SG.XKB(x.vec, sigma2d)
  
  avg.loss.mu.est.SURE.SG = mean((mud-mu.est.SURE.SG)^2)
  
  #************************************
  
  #*****SURE.SM  method*****************
  
  mu.est.SURE.SM = thetahat.SURE.SM.XKB(x.vec, sigma2d)
  
  avg.loss.mu.est.SURE.SM = mean((mud-mu.est.SURE.SM)^2)
  

  #************************************
  
  #************************************************************  
  
  # group_linear: num bins = n^1/3
  mu.est.gl = grouplinear(x.vec, sigma2d)
  
  avg.loss.mu.est.gl = mean((mud-mu.est.gl)^2)
  
  
  #************************************************************
  
  # group_linear: sure(equal-bins)
  mu.est.gl.SURE = grouplinear.sure(x.vec, sigma2d, kmax=min(ceiling(q^(1/3)/0.8),q))
  
  avg.loss.mu.est.gl.SURE = mean((mud-mu.est.gl.SURE)^2)
  
  #************************************************************
  
  # dynamic_group_linear_all_division
  mu.est.gl.dynamic = GroupSure(x.vec, sigma2d)
  
  avg.loss.mu.est.gl.dynamic = mean((mud-mu.est.gl.dynamic)^2)
  
  
  # dynamic_group_linear_with_minimum_bin_size_constraint
  mu.est.gl.dynamicMin = GroupSureMin(x.vec, sigma2d, q^(2/3)*0.8)
  
  avg.loss.mu.est.gl.dynamicMin = mean((mud-mu.est.gl.dynamicMin)^2)
  
  mu.est.gl.dynamicMin2 = GroupSureMin(x.vec, sigma2d, q^(2/3))
  
  avg.loss.mu.est.gl.dynamicMin2 = mean((mud-mu.est.gl.dynamicMin2)^2)
  
  mu.est.gl.dynamicMin3 = GroupSureMin(x.vec, sigma2d, q^(2/3)*1.2)
  
  avg.loss.mu.est.gl.dynamicMin3 = mean((mud-mu.est.gl.dynamicMin3)^2)
  

  #**************************************************************
 
  #*********NG mixture method**********************************
  
  MCMCoutput_ng = DPMM.ng.sknown(x.vec, 1/sigma2d, gamma.pi=0.1, 
                                 k=10, Burnin=3000, Simsize=5000, mu.grid=mu.grid, prec.grid=prec.grid, 
                                 hyperparameters = c('Default'))
  
  avg.loss.mu.est.DPMM_ng = mean((MCMCoutput_ng$mu.est-mud)^2)
  

  #************************************************************
  
  if(b%%50==0){
    print(b)
  }
  
  #************************************************************
  
  
  #*************************************************************
  
  obj1 = 
    c(q,
      avg.loss.mu.est.true.oracle,
      avg.loss.mu.est.sample.mean,
      avg.loss.mu.est.EBMLE,
      avg.loss.mu.est.EBMOM,
      avg.loss.mu.est.JS,
      avg.loss.mu.est.oracle,
      avg.loss.mu.est.SURE.G,
      avg.loss.mu.est.SURE.M,
      avg.loss.mu.est.SURE.SG,
      avg.loss.mu.est.SURE.SM,
      avg.loss.mu.est.gl,
      avg.loss.mu.est.gl.SURE,
      avg.loss.mu.est.gl.dynamic,
      avg.loss.mu.est.gl.dynamicMin,
      avg.loss.mu.est.gl.dynamicMin2,
      avg.loss.mu.est.gl.dynamicMin3,
      avg.loss.mu.est.DPMM_ng)
  
  obj1
}

distinctq = unique(as.vector(alloutput_diffq[1,]))
paperoutput = array(0,c(dim(alloutput_diffq)[1],length(distinctq)))

for(kk in 1:length(distinctq)){
  paperoutput[,kk] = apply(alloutput_diffq[,which(as.vector(alloutput_diffq[1,])==distinctq[kk])],1,mean)
}

finalans = data.frame(c("q",
                        "avg.loss.mu.est.true.oracle",
                        "avg.loss.mu.est.sample.mean",
                        "avg.loss.mu.est.EBMLE",
                        "avg.loss.mu.est.EBMOM",
                        "avg.loss.mu.est.JS",
                        "avg.loss.mu.est.oracle",
                        "avg.loss.mu.est.SURE.G",
                        "avg.loss.mu.est.SURE.M",
                        "avg.loss.mu.est.SURE.SG",
                        "avg.loss.mu.est.SURE.SM",
                        "avg.loss.mu.est.gl",
                        "avg.loss.mu.est.gl.SURE",
                        "avg.loss.mu.est.gl.dynamic",
                        "avg.loss.mu.est.gl.dynamicMin",
                        "avg.loss.mu.est.gl.dynamicMin2",
                        "avg.loss.mu.est.gl.dynamicMin3",
                        "avg.loss.mu.est.DPMM_ng"),
                      paperoutput)

names(finalans) = c("measure", rep("value",((q.upper.lim-q.lower.lim)/q.interval+1)))

print(finalans)

write.csv(finalans, file=outputfilecsv)

save(list = ls(), file = outputfile)
