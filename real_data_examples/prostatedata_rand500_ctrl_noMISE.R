rm(list = ls())

B = 100

library(foreach)
library(doParallel)

setwd("~/Desktop/finalRdir/DPMM/v10")

outputfile = "prostate_500ctrl_noMISE.Rdata"
outputfilecsv = "prostate_500ctrl_noMISE.csv"

load("dataprostate.RData")

source("DPMMfunction.R")
source("SUREmethods.R")
source("grouplinearfunction_all.R")

controldata = prostatedata[,1:50]
treatmentdata = prostatedata[,51:102]

data = controldata

xbar.all.vec = apply(data,1,mean)
S2.all.vec = apply(data,1,var)

n1 = 100
n2 = 100

mu.grid = seq(min(xbar.all.vec)-IQR(xbar.all.vec), max(xbar.all.vec)+IQR(xbar.all.vec), 
              length.out=2*n1+1)[seq(2,2*n1,2)]
sigma2.grid = seq(max(min(S2.all.vec)-IQR(S2.all.vec),0.001), max(S2.all.vec)+IQR(S2.all.vec), 
                  length.out=2*n2+1)[seq(2,2*n2,2)]

prec.grid = 1/sigma2.grid


cl = min(20,B)
registerDoParallel(cl)

alloutput = foreach(b = 1:B,.combine='cbind',.inorder=FALSE) %dopar% {
  
  set.seed(1234+43*b)
  
  rows = sample(dim(data)[1],500)
  cols = sample(dim(data)[2],3)
  
  data1 = data[rows,cols]
  mud = xbar.all.vec[rows]
  sigma2d = S2.all.vec[rows]
  q = dim(data1)[1]
  n = dim(data1)[2]
  ids = rep(1:q,n)
  x.vec = c(data1)
  xbar.vec = tapply(x.vec,ids,mean)
  S2.vec = tapply(x.vec,ids,var)
  ni.vec = as.vector(table(ids))
  S2.vec.xbar = S2.vec/ni.vec

  #naive density estimation
  avg.loss.mu.est.sample.mean = mean((mud-xbar.vec)^2)
  avg.loss.sigma2d.est.sample.var = mean((sigma2d-S2.vec)^2) 

  
  #********EBMLE***********
  
  mu.est.EBMLE = thetahat.EBMLE.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.EBMLE = mean((mud-mu.est.EBMLE)^2)
  
  #**********************
  
  
  #*******EBMOM**********
  
  mu.est.EBMOM = thetahat.EBMOM.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.EBMOM = mean((mud-mu.est.EBMOM)^2)
  
  #*****************************
  
  #*************JS estimate*****
  
  mu.est.JS = thetahat.JS.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.JS = mean((mud-mu.est.JS)^2)
  
  #****************************
  
  #********Oracle**************
  
  mu.est.oracle = thetahat.oracle.XKB(xbar.vec, S2.vec.xbar, mud)
  
  avg.loss.mu.est.oracle = mean((mud-mu.est.oracle)^2)
  

  #************************************
  
  #********large class oracle**************
  
  m.or = mean(mud)
  or.loss = vector()
  
  for (i in 1:q){
    if(mud[i] < xbar.vec[i] & mud[i] > m.or){
      or.loss[i] = 0
    }
    else if(mud[i] > xbar.vec[i] & mud[i] < m.or){
      or.loss[i] = 0
    } else{
      or.loss[i] = min((xbar.vec-mud)^2,(mud-m.or)^2)
    }
  }
  avg.loss.mu.est.true.oracle = mean(or.loss)
  
  #*****SURE.G method*****************
  
  mu.est.SURE.G = thetahat.SURE.G.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.SURE.G = mean((mud-mu.est.SURE.G)^2)
  
 
  #************************************
  
  #*****SURE.M method*****************
  
  mu.est.SURE.M = thetahat.SURE.M.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.SURE.M = mean((mud-mu.est.SURE.M)^2)
  
  #************************************
  
  #*****SURE.SG  method*****************
  
  mu.est.SURE.SG = thetahat.SURE.SG.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.SURE.SG = mean((mud-mu.est.SURE.SG)^2)
  
  #************************************
  
  #*****SURE.SM  method*****************
  
  mu.est.SURE.SM = thetahat.SURE.SM.XKB(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.SURE.SM = mean((mud-mu.est.SURE.SM)^2)
  
  
  #*************************************************
  
  #************************************************************  
  
  # group_linear: num bins = n^1/3
  mu.est.gl = grouplinear(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.gl = mean((mud-mu.est.gl)^2)
  
  #************************************************************
  
  # group_linear: sure(equal-bins)
  mu.est.gl.SURE = grouplinear.sure(xbar.vec, S2.vec.xbar, kmax=min(ceiling(q^(1/3)/0.8),q))
  
  avg.loss.mu.est.gl.SURE = mean((mud-mu.est.gl.SURE)^2)
  
  #************************************************************
  
  # dynamic_group_linear_all_division
  mu.est.gl.dynamic = GroupSure(xbar.vec, S2.vec.xbar)
  
  avg.loss.mu.est.gl.dynamic = mean((mud-mu.est.gl.dynamic)^2)
  
  
  # dynamic_group_linear_with_minimum_bin_size_constraint
  mu.est.gl.dynamicMin = GroupSureMin(xbar.vec, S2.vec.xbar,q^(2/3)*0.8)
  
  avg.loss.mu.est.gl.dynamicMin = mean((mud-mu.est.gl.dynamicMin)^2)
  
  
  mu.est.gl.dynamicMin2 = GroupSureMin(xbar.vec, S2.vec.xbar,q^(2/3))
  
  avg.loss.mu.est.gl.dynamicMin2 = mean((mud-mu.est.gl.dynamicMin2)^2)
  

  mu.est.gl.dynamicMin3 = GroupSureMin(xbar.vec, S2.vec.xbar,q^(2/3)*1.2)
  
  avg.loss.mu.est.gl.dynamicMin3 = mean((mud-mu.est.gl.dynamicMin3)^2)
  
  #************************************************************
  
  
  #**********SURE.M Double*************************************
  
  mu.est.SURE.M.Double.object = estimates.SURE.M.JLPZ(xbar.vec,S2.vec.xbar,n)
  
  mu.est.SURE.M.Double = mu.est.SURE.M.Double.object$mu.est
  
  sigma2d.est.SURE.M.Double = mu.est.SURE.M.Double.object$sigma2.est
  
  avg.loss.mu.est.SURE.M.Double = mean((mud-mu.est.SURE.M.Double)^2)
  avg.loss.sigma2d.est.SURE.M.Double = mean((sigma2d-sigma2d.est.SURE.M.Double)^2)
  
  #*************************************************************
  
  #**********SURE.G Double**************************************
  
  
  mu.est.SURE.G.Double.object = estimates.SURE.G.JLPZ(xbar.vec,S2.vec.xbar,n)
  
  mu.est.SURE.G.Double = mu.est.SURE.G.Double.object$mu.est
  
  sigma2d.est.SURE.G.Double = mu.est.SURE.G.Double.object$sigma2.est
  
  avg.loss.mu.est.SURE.G.Double = mean((mud-mu.est.SURE.G.Double)^2)
  avg.loss.sigma2d.est.SURE.G.Double = mean((sigma2d-sigma2d.est.SURE.G.Double)^2)
  
  #*************************************************************
  
  #**********SURE.SM Double**************************************
  
  
  mu.est.SURE.SM.Double.object = estimates.SURE.SM.JLPZ(xbar.vec,S2.vec.xbar,n)
  
  mu.est.SURE.SM.Double = mu.est.SURE.SM.Double.object$mu.est
  
  sigma2d.est.SURE.SM.Double = mu.est.SURE.SM.Double.object$sigma2.est
  
  avg.loss.mu.est.SURE.SM.Double = mean((mud-mu.est.SURE.SM.Double)^2)
  avg.loss.sigma2d.est.SURE.SM.Double = mean((sigma2d-sigma2d.est.SURE.SM.Double)^2)
  
  #*************************************************************
  
  #**********SURE.SG Double**************************************
  
  mu.est.SURE.SG.Double.object = estimates.SURE.SG.JLPZ(xbar.vec,S2.vec.xbar,n)
  
  mu.est.SURE.SG.Double = mu.est.SURE.SG.Double.object$mu.est
  
  sigma2d.est.SURE.SG.Double = mu.est.SURE.SG.Double.object$sigma2.est
  
  avg.loss.mu.est.SURE.SG.Double = mean((mud-mu.est.SURE.SG.Double)^2)
  avg.loss.sigma2d.est.SURE.SG.Double = mean((sigma2d-sigma2d.est.SURE.SG.Double)^2)
  
  #*************************************************************
  
  #*********NIG mixture method**********************************
  
  MCMCoutput = DPMM.nig(x.vec, ids, gamma.pi=0.1, 
                        k=10, Burnin=5000, Simsize=5000, mu.grid=mu.grid, sigma2.grid=sigma2.grid, 
                        hyperparameters = c('Default'))
  
  avg.loss.mu.est.DPMM = mean((MCMCoutput$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM = mean((MCMCoutput$sigma2.est-sigma2d)^2)
  
  #*********NG mixture method**********************************
  
  
  MCMCoutput_ng = DPMM.ng(x.vec, ids, gamma.pi=0.1, 
                          k=10, Burnin=5000, Simsize=5000, mu.grid=mu.grid, prec.grid=prec.grid, 
                          hyperparameters = c('Default'))
  
  avg.loss.mu.est.DPMM_ng = mean((MCMCoutput_ng$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM_ng = mean((1/MCMCoutput_ng$prec.est-sigma2d)^2)


  #*************************************************************

  #*********NIG mixture method**********************************
  
  MCMCoutput1 = DPMM.nig(x.vec, ids, gamma.pi=0.1, 
                        k=1, Burnin=5000, Simsize=5000, mu.grid=mu.grid, sigma2.grid=sigma2.grid, 
                        hyperparameters = c('Default'))
  
  avg.loss.mu.est.DPMM1 = mean((MCMCoutput1$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM1 = mean((MCMCoutput1$sigma2.est-sigma2d)^2)
  
  #*********NG mixture method**********************************
  
  
  MCMCoutput_ng1 = DPMM.ng(x.vec, ids, gamma.pi=0.1, 
                          k=1, Burnin=5000, Simsize=5000, mu.grid=mu.grid, prec.grid=prec.grid, 
                          hyperparameters = c('Default'))
  
  avg.loss.mu.est.DPMM_ng1 = mean((MCMCoutput_ng1$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM_ng1 = mean((1/MCMCoutput_ng1$prec.est-sigma2d)^2)
  

  
  #*************************************************************
  
  if(b%%50==0){
    print(b)
  }
  
  #************************************************************
  
  
  #*************************************************************
  
  obj1 = 
    c(avg.loss.mu.est.true.oracle,
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
      avg.loss.mu.est.SURE.M.Double,
      avg.loss.mu.est.SURE.G.Double,
      avg.loss.mu.est.SURE.SG.Double,
      avg.loss.mu.est.SURE.SM.Double,
      avg.loss.mu.est.DPMM,
      avg.loss.mu.est.DPMM_ng,
      avg.loss.mu.est.DPMM1,
      avg.loss.mu.est.DPMM_ng1,
      
      avg.loss.sigma2d.est.sample.var,
      avg.loss.sigma2d.est.SURE.M.Double,
      avg.loss.sigma2d.est.SURE.G.Double,
      avg.loss.sigma2d.est.SURE.SG.Double,
      avg.loss.sigma2d.est.SURE.SM.Double,
      avg.loss.sigma2d.est.DPMM,
      avg.loss.sigma2d.est.DPMM_ng,
      avg.loss.sigma2d.est.DPMM1,
      avg.loss.sigma2d.est.DPMM_ng1
    )
  
  
  obj1
  
}

SUREcomparision = apply(alloutput,1,mean)

finalans = data.frame(c("avg.loss.mu.est.true.oracle",
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
                        "avg.loss.mu.est.SURE.M.Double",
                        "avg.loss.mu.est.SURE.G.Double",
                        "avg.loss.mu.est.SURE.SG.Double",
                        "avg.loss.mu.est.SURE.SM.Double",
                        "avg.loss.mu.est.DPMM",
                        "avg.loss.mu.est.DPMM_ng",
                        "avg.loss.mu.est.DPMM1",
                        "avg.loss.mu.est.DPMM_ng1",
                        
                        "avg.loss.sigma2d.est.sample.var",
                        "avg.loss.sigma2d.est.SURE.M.Double",
                        "avg.loss.sigma2d.est.SURE.G.Double",
                        "avg.loss.sigma2d.est.SURE.SG.Double",
                        "avg.loss.sigma2d.est.SURE.SM.Double",
                        "avg.loss.sigma2d.est.DPMM",
                        "avg.loss.sigma2d.est.DPMM_ng",
                        "avg.loss.sigma2d.est.DPMM1",
                        "avg.loss.sigma2d.est.DPMM_ng1"),
                      SUREcomparision)

names(finalans) = c("measure", "value")

print(finalans)

write.csv(finalans,file = outputfilecsv)

save(list = ls(), file = outputfile)

