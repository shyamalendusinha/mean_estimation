
rm(list = ls())

library(foreach)
library(doParallel)

library(invgamma)

#********Data Storing**************************

outputnumber = '11-box-sort'

setwd("~/mean_estimation/sensitivity_gamma_exam11")
outputfilecsv = paste("diffq",outputnumber,".csv",sep="")
outputfile = paste("diffq",outputnumber,".Rdata",sep="")

#************Data Simulation parameters*********************
Bq = 100
q.upper.lim = 500
q.lower.lim = 20
q.interval = 40
B = Bq*((q.upper.lim-q.lower.lim)/q.interval+1)


#****Data simulation parameters******************
kd = 2
meand = c(2,10)
lambdad = c(2,4)
sigmashape = c(5,3)
sigmarate = c(2,3)
pid = c(0.95,0.05)

r.biv.density = function(q){
  mud = sigma2d =vector()
  zd = sample(kd, q, prob=pid, replace=TRUE)
  for (i in 1:q){
    sigma2d[i] = rinvgamma(1, shape=sigmashape[zd[i]], rate=sigmarate[zd[i]])
    mud[i] = rnorm(1, mean=meand[zd[i]], sd=sqrt(sigma2d[i]/lambdad[zd[i]]))
  }
  return(list(mud=mud, sigma2d=sigma2d))
}

source("DPMMfunction_sensitivity.R")
# source("SUREmethods.R")
# source("grouplinearfunction_all.R")

cl = min(40,B)
registerDoParallel(cl)

alloutput_diffq = foreach(b = 1:B,.combine='cbind',.inorder=FALSE) %dopar% {
  
  set.seed(1234+47*b)
  
  n = 4
  q = ceiling(b/Bq)*q.interval+(q.lower.lim-q.interval)
  ni.vec = rep(n,q)
  n1 = 100
  n2 = 100
  
  random.sample = r.biv.density(q)
  mud = random.sample$mud
  sigma2d = random.sample$sigma2d
  errord = rnorm(n*q)
  x.vec = rep(mud,ni.vec)+sqrt(rep(sigma2d,ni.vec))*errord
  ids = rep(1:q,ni.vec)
  xbar.vec = tapply(x.vec,ids,mean)
  S2.vec = tapply(x.vec,ids,var)
  sigma2d.xbar = sigma2d/ni.vec
  S2.vec.xbar = S2.vec/ni.vec
  mu.grid = seq(min(xbar.vec)-IQR(xbar.vec), max(xbar.vec)+IQR(xbar.vec), 
                length.out=2*n1+1)[seq(2,2*n1,2)]
  sigma2.grid = seq(max(min(S2.vec)-IQR(S2.vec),0.001), max(S2.vec)+IQR(S2.vec), 
                    length.out=2*n2+1)[seq(2,2*n2,2)]
  
  
  # #naive density estimation
  # avg.loss.mu.est.sample.mean = mean((mud-xbar.vec)^2)
  # avg.loss.sigma2d.est.sample.var = mean((sigma2d-S2.vec)^2) 
  # 
  # 
  # #********EBMLE***********
  # 
  # mu.est.EBMLE = thetahat.EBMLE.XKB(xbar.vec, S2.vec.xbar)
  # 
  # avg.loss.mu.est.EBMLE = mean((mud-mu.est.EBMLE)^2)
  # 
  # #**********************
  # 
  # 
  # #*******EBMOM**********
  # 
  # mu.est.EBMOM = thetahat.EBMOM.XKB(xbar.vec, S2.vec.xbar)
  # 
  # avg.loss.mu.est.EBMOM = mean((mud-mu.est.EBMOM)^2)
  # 
  # #*****************************
  # 
  # #*************JS estimate*****
  # 
  # mu.est.JS = thetahat.JS.XKB(xbar.vec, S2.vec.xbar)
  # 
  # avg.loss.mu.est.JS = mean((mud-mu.est.JS)^2)
  # 
  # #****************************
  # 
  # #********Oracle**************
  # 
  # mu.est.oracle = thetahat.oracle.XKB(xbar.vec, S2.vec.xbar, mud)
  # 
  # avg.loss.mu.est.oracle = mean((mud-mu.est.oracle)^2)
  # 
  # 
  # #************************************
  # 
  # #********large class oracle**************
  # 
  # m.or = mean(mud)
  # or.loss = vector()
  # 
  # for (i in 1:q){
  #   if(mud[i] < xbar.vec[i] & mud[i] > m.or){
  #     or.loss[i] = 0
  #   }
  #   else if(mud[i] > xbar.vec[i] & mud[i] < m.or){
  #     or.loss[i] = 0
  #   } else{
  #     or.loss[i] = min((xbar.vec-mud)^2,(mud-m.or)^2)
  #   }
  # }
  # avg.loss.mu.est.true.oracle = mean(or.loss)
  # 
  # #*****SURE.G method*****************
  # 
  # mu.est.SURE.G = thetahat.SURE.G.XKB(xbar.vec, S2.vec.xbar)
  # 
  # avg.loss.mu.est.SURE.G = mean((mud-mu.est.SURE.G)^2)
  # 
  # #************************************
  # 
  # #*****SURE.M method*****************
  # 
  # mu.est.SURE.M = thetahat.SURE.M.XKB(xbar.vec, S2.vec.xbar)
  # 
  # avg.loss.mu.est.SURE.M = mean((mud-mu.est.SURE.M)^2)
  # 
  # #************************************
  # 
  # #*****SURE.SG  method*****************
  # 
  # mu.est.SURE.SG = thetahat.SURE.SG.XKB(xbar.vec, S2.vec.xbar)
  # 
  # avg.loss.mu.est.SURE.SG = mean((mud-mu.est.SURE.SG)^2)
  # 
  # #************************************
  # 
  # #*****SURE.SM  method*****************
  # 
  # mu.est.SURE.SM = thetahat.SURE.SM.XKB(xbar.vec, S2.vec.xbar)
  # 
  # avg.loss.mu.est.SURE.SM = mean((mud-mu.est.SURE.SM)^2)
  # 
  # 
  # #************************************
  # 
  # #************************************************************  
  # 
  # # group_linear: num bins = n^1/3
  # mu.est.gl = grouplinear(xbar.vec, S2.vec.xbar)
  # 
  # avg.loss.mu.est.gl = mean((mud-mu.est.gl)^2)
  # 
  # 
  # #************************************************************
  # 
  # # group_linear: sure(equal-bins)
  # mu.est.gl.SURE = grouplinear.sure(xbar.vec, S2.vec.xbar, kmax=min(ceiling(q^(1/3)/0.8),q))
  # 
  # avg.loss.mu.est.gl.SURE = mean((mud-mu.est.gl.SURE)^2)
  # 
  # #************************************************************
  # 
  # # dynamic_group_linear_all_division
  # mu.est.gl.dynamic = GroupSure(xbar.vec, S2.vec.xbar)
  # 
  # avg.loss.mu.est.gl.dynamic = mean((mud-mu.est.gl.dynamic)^2)
  # 
  # 
  # 
  # # dynamic_group_linear_with_minimum_bin_size_constraint
  # mu.est.gl.dynamicMin = GroupSureMin(xbar.vec, S2.vec.xbar,q^(2/3)*0.8)
  # 
  # avg.loss.mu.est.gl.dynamicMin = mean((mud-mu.est.gl.dynamicMin)^2)
  # 
  # 
  # mu.est.gl.dynamicMin2 = GroupSureMin(xbar.vec, S2.vec.xbar,q^(2/3))
  # 
  # avg.loss.mu.est.gl.dynamicMin2 = mean((mud-mu.est.gl.dynamicMin2)^2)
  # 
  # 
  # mu.est.gl.dynamicMin3 = GroupSureMin(xbar.vec, S2.vec.xbar,q^(2/3)*1.2)
  # 
  # avg.loss.mu.est.gl.dynamicMin3 = mean((mud-mu.est.gl.dynamicMin3)^2)
  # 
  # 
  # #************************************************************
  # 
  # 
  # #**********SURE.M Double**************************
  # 
  # mu.est.SURE.M.Double.object = estimates.SURE.M.JLPZ(xbar.vec,S2.vec.xbar,n)
  # 
  # mu.est.SURE.M.Double = mu.est.SURE.M.Double.object$mu.est
  # 
  # sigma2d.est.SURE.M.Double = mu.est.SURE.M.Double.object$sigma2.est
  # 
  # avg.loss.mu.est.SURE.M.Double = mean((mud-mu.est.SURE.M.Double)^2)
  # avg.loss.sigma2d.est.SURE.M.Double = mean((sigma2d-sigma2d.est.SURE.M.Double)^2)
  # 
  # #*************************************************************
  # 
  # #**********SURE.G Double**************************************
  # 
  # 
  # mu.est.SURE.G.Double.object = estimates.SURE.G.JLPZ(xbar.vec,S2.vec.xbar,n)
  # 
  # mu.est.SURE.G.Double = mu.est.SURE.G.Double.object$mu.est
  # 
  # sigma2d.est.SURE.G.Double = mu.est.SURE.G.Double.object$sigma2.est
  # 
  # avg.loss.mu.est.SURE.G.Double = mean((mud-mu.est.SURE.G.Double)^2)
  # avg.loss.sigma2d.est.SURE.G.Double = mean((sigma2d-sigma2d.est.SURE.G.Double)^2)
  # 
  # 
  # #*************************************************************
  # 
  # #**********SURE.SM Double**************************************
  # 
  # 
  # mu.est.SURE.SM.Double.object = estimates.SURE.SM.JLPZ(xbar.vec,S2.vec.xbar,n)
  # 
  # mu.est.SURE.SM.Double = mu.est.SURE.SM.Double.object$mu.est
  # 
  # sigma2d.est.SURE.SM.Double = mu.est.SURE.SM.Double.object$sigma2.est
  # 
  # avg.loss.mu.est.SURE.SM.Double = mean((mud-mu.est.SURE.SM.Double)^2)
  # avg.loss.sigma2d.est.SURE.SM.Double = mean((sigma2d-sigma2d.est.SURE.SM.Double)^2)
  # 
  # 
  # #*************************************************************
  # 
  # #**********SURE.SG Double**************************************
  # 
  # mu.est.SURE.SG.Double.object = estimates.SURE.SG.JLPZ(xbar.vec,S2.vec.xbar,n)
  # 
  # mu.est.SURE.SG.Double = mu.est.SURE.SG.Double.object$mu.est
  # 
  # sigma2d.est.SURE.SG.Double = mu.est.SURE.SG.Double.object$sigma2.est
  # 
  # avg.loss.mu.est.SURE.SG.Double = mean((mud-mu.est.SURE.SG.Double)^2)
  # avg.loss.sigma2d.est.SURE.SG.Double = mean((sigma2d-sigma2d.est.SURE.SG.Double)^2)
  # 
  # 
  # #*************************************************************
  

  #*********NG mixture method**********************************
  
  prec.grid = 1/sigma2.grid
  
  MCMCoutput_ng1 = DPMM.ng_sensitivity(x.vec, ids, gamma.pi=0.1, 
                          k=10, Burnin=5000, Simsize=5000, mu.grid=mu.grid, prec.grid=prec.grid, 
                          hyperparameters = c('Default'))

  MCMCoutput_ng2 = DPMM.ng_sensitivity(x.vec, ids, gamma.pi=10, 
                          k=10, Burnin=5000, Simsize=5000, mu.grid=mu.grid, prec.grid=prec.grid, 
                          hyperparameters = c('Default'))

  MCMCoutput_ng3 = DPMM.ng_sensitivity(x.vec, ids, gamma.pi=50, 
                          k=10, Burnin=5000, Simsize=5000, mu.grid=mu.grid, prec.grid=prec.grid, 
                          hyperparameters = c('Default'))

  MCMCoutput_ng4 = DPMM.ng_sensitivity(x.vec, ids, gamma.pi=100, 
                          k=10, Burnin=5000, Simsize=5000, mu.grid=mu.grid, prec.grid=prec.grid, 
                          hyperparameters = c('Default'))
  
  avg.loss.mu.est.DPMM_ng1 = mean((MCMCoutput_ng1$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM_ng1 = mean((1/MCMCoutput_ng1$prec.est-sigma2d)^2)
  prob.3_ng1 = MCMCoutput_ng1$pi8.est
  prob.2_ng1 = MCMCoutput_ng1$pi9.est

  avg.loss.mu.est.DPMM_ng2 = mean((MCMCoutput_ng2$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM_ng2 = mean((1/MCMCoutput_ng2$prec.est-sigma2d)^2)
  prob.3_ng2 = MCMCoutput_ng2$pi8.est
  prob.2_ng2 = MCMCoutput_ng2$pi9.est
  
  avg.loss.mu.est.DPMM_ng3 = mean((MCMCoutput_ng3$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM_ng3 = mean((1/MCMCoutput_ng3$prec.est-sigma2d)^2)
  prob.3_ng3 = MCMCoutput_ng3$pi8.est
  prob.2_ng3 = MCMCoutput_ng3$pi9.est

  avg.loss.mu.est.DPMM_ng4 = mean((MCMCoutput_ng4$mu.est-mud)^2)
  avg.loss.sigma2d.est.DPMM_ng4 = mean((1/MCMCoutput_ng4$prec.est-sigma2d)^2)
  prob.3_ng4 = MCMCoutput_ng4$pi8.est
  prob.2_ng4 = MCMCoutput_ng4$pi9.est
  
  
  if(b%%500==0){
    print(b)
  }
  
  #************************************************************
  
  
  #*************************************************************
  
  obj1 = 
    c(q,
      # avg.loss.mu.est.true.oracle,
      # avg.loss.mu.est.sample.mean,
      # avg.loss.mu.est.EBMLE,
      # avg.loss.mu.est.EBMOM,
      # avg.loss.mu.est.JS,
      # avg.loss.mu.est.oracle,
      # avg.loss.mu.est.SURE.G,
      # avg.loss.mu.est.SURE.M,
      # avg.loss.mu.est.SURE.SG,
      # avg.loss.mu.est.SURE.SM,
      # avg.loss.mu.est.gl,
      # avg.loss.mu.est.gl.SURE,
      # avg.loss.mu.est.gl.dynamic,
      # avg.loss.mu.est.gl.dynamicMin,
      # avg.loss.mu.est.gl.dynamicMin2,
      # avg.loss.mu.est.gl.dynamicMin3,
      # avg.loss.mu.est.SURE.M.Double,
      # avg.loss.mu.est.SURE.G.Double,
      # avg.loss.mu.est.SURE.SG.Double,
      # avg.loss.mu.est.SURE.SM.Double,
      avg.loss.mu.est.DPMM_ng1,
      avg.loss.mu.est.DPMM_ng2,
      avg.loss.mu.est.DPMM_ng3,
      avg.loss.mu.est.DPMM_ng4,
      
      # avg.loss.sigma2d.est.sample.var,
      # avg.loss.sigma2d.est.SURE.M.Double,
      # avg.loss.sigma2d.est.SURE.G.Double,
      # avg.loss.sigma2d.est.SURE.SG.Double,
      # avg.loss.sigma2d.est.SURE.SM.Double,
      avg.loss.sigma2d.est.DPMM_ng1,
      avg.loss.sigma2d.est.DPMM_ng2,
      avg.loss.sigma2d.est.DPMM_ng3,
      avg.loss.sigma2d.est.DPMM_ng4,

      prob.2_ng1,
      prob.2_ng2,
      prob.2_ng3,
      prob.2_ng4,

      prob.3_ng1,
      prob.3_ng2,
      prob.3_ng3,
      prob.3_ng4
)
  
  
  obj1
  
}


alloutput_diffq1 = t(alloutput_diffq)

colnames(alloutput_diffq1) = c("q",
                    # "avg.loss.mu.est.true.oracle",
                    # "avg.loss.mu.est.sample.mean",
                    # "avg.loss.mu.est.EBMLE",
                    # "avg.loss.mu.est.EBMOM",
                    # "avg.loss.mu.est.JS",
                    # "avg.loss.mu.est.oracle",
                    # "avg.loss.mu.est.SURE.G",
                    # "avg.loss.mu.est.SURE.M",
                    # "avg.loss.mu.est.SURE.SG",
                    # "avg.loss.mu.est.SURE.SM",
                    # "avg.loss.mu.est.gl",
                    # "avg.loss.mu.est.gl.SURE",
                    # "avg.loss.mu.est.gl.dynamic",
                    # "avg.loss.mu.est.gl.dynamicMin",
                    # "avg.loss.mu.est.gl.dynamicMin2",
                    # "avg.loss.mu.est.gl.dynamicMin3",
                    # "avg.loss.mu.est.SURE.M.Double",
                    # "avg.loss.mu.est.SURE.G.Double",
                    # "avg.loss.mu.est.SURE.SG.Double",
                    # "avg.loss.mu.est.SURE.SM.Double",
                    "avg.loss.mu.est.DPMM_ng1",
                    "avg.loss.mu.est.DPMM_ng2",
                    "avg.loss.mu.est.DPMM_ng3",
                    "avg.loss.mu.est.DPMM_ng4",
                    
                    # "avg.loss.sigma2d.est.sample.var",
                    # "avg.loss.sigma2d.est.SURE.M.Double",
                    # "avg.loss.sigma2d.est.SURE.G.Double",
                    # "avg.loss.sigma2d.est.SURE.SG.Double",
                    # "avg.loss.sigma2d.est.SURE.SM.Double",
                    "avg.loss.sigma2d.est.DPMM_ng1",
                    "avg.loss.sigma2d.est.DPMM_ng2",
                    "avg.loss.sigma2d.est.DPMM_ng3",
                    "avg.loss.sigma2d.est.DPMM_ng4",

                    "prob.2_ng1",
                    "prob.2_ng2",
                    "prob.2_ng3",
                    "prob.2_ng4",
                    
                    "prob.3_ng1",
                    "prob.3_ng2",
                    "prob.3_ng3",
                    "prob.3_ng4"
)

# print(alloutput_diffq1)

write.csv(alloutput_diffq1, file=outputfilecsv)

save(list = ls(), file = outputfile)
