
# The whole code was downloaded from https://github.com/MaZhuang/grouplinear/shuffle_final.R.
# Only modified to add normal-gamma mixture model and to run parallel on different cluster.

setwd("~/mean_estimation")

rm(list = ls())

library(foreach)
library(doParallel)

bat.perm = function(){
  bat = bat.raw
  bat$N1 = bat$AB.4. + bat$AB.5. + bat$AB.6.  # total number at-bats for 1st period
  bat$N2 = bat$AB.7. + bat$AB.8. + bat$AB.9.10.  # total number at-bats for 2nd period
  bat$H1 = bat$H.4. + bat$H.5. + bat$H.6.  # total number hits for 1st period
  bat$H2 = bat$H.7. + bat$H.8. + bat$H.9.10.  # total number hits for 2nd period
  # bat$R1 = bat$H1/bat$N1  # batting avg for 1st period
  # bat$R2 = bat$H2/bat$N2  # batting avg for 2nd period
  # bat$X1 = asin(sqrt((bat$H1+1/4)/(bat$N1+1/2)))  # transformed batting avg for 1st period
  # bat$X2 = asin(sqrt((bat$H2+1/4)/(bat$N2+1/2)))  # transformed batting avg for 2nd period
  bat = bat[bat$N1 > 10,]  # keep only records with N1>=11  
  
  bat$H1.perm = NA
  for(i in 1:dim(bat)[1]){
    bat$H1.perm[i] = rhyper(nn=1,m=bat$H1[i] + bat$H2[i],
                             n=bat$N1[i] + bat$N2[i] -bat$H1[i] - bat$H2[i],k=bat$N1[i])
  }
  bat$H2.perm = bat$H1 + bat$H2 - bat$H1.perm
  # head(cbind(bat$H1,bat$H1.perm,bat$H2,bat$H2.perm))
  bat$H1 = bat$H1.perm 
  bat$H2 = bat$H2.perm 
  
  bat$R1 = bat$H1/bat$N1  # batting avg for 1st period
  bat$R2 = bat$H2/bat$N2  # batting avg for 2nd period
  bat$X1 = asin(sqrt((bat$H1+1/4)/(bat$N1+1/2)))  # transformed batting avg for 1st period
  bat$X2 = asin(sqrt((bat$H2+1/4)/(bat$N2+1/2)))  # transformed batting avg for 2nd period
  bat =  bat[,c('First.Name','Last.Name','Pitcher.','N1','N2','H1','H2','X1','X2')]
}

N = 1000 # num shuffling rounds
source("~/mean_estimation/all_R_functions/DPMMfunction.R")
source("~/mean_estimation/all_R_functions/grouplinearfunction_all.R")

bat.raw = read.table('~/mean_estimation/real_data_examples/Brown_batting_data.txt', header=TRUE, sep=",", quote="")

cl = min(30,N)
registerDoParallel(cl)

alloutput = foreach(j = 1:N,.combine='cbind',.inorder=FALSE) %dopar% {
  
  set.seed(1234+47*j)
  bat = bat.perm()
  index=order(bat$N1,decreasing=TRUE)
  bat=bat[index,]
  n = dim(bat)[1]
  # estimating TSE for various estimators
  # run: functions.R(current folder), functions_XKB.R
  ind = bat$N2 > 10  # indicator for records with N2>=11 (among those with N1>=11)
  tse.zero = sum(((bat$X2 - bat$X1)^2 - 1/(4*bat$N2))[ind])
  
  # grand mean
  tse.delta.gm = sum(((bat$X2 - mean(bat$X1))^2 - 1/(4*bat$N2))[ind])
  tse.gm = tse.delta.gm/tse.zero
  
  # James-Stein
  delta.JS = JS(bat$X1,1/(4*bat$N1))
  tse.delta.JS = sum(((bat$X2 - delta.JS)^2 - 1/(4*bat$N2))[ind])
  tse.JS = tse.delta.JS/tse.zero
  
  # XKB theta.hat.M
  delta.M = thetahat.M(bat$X1,1/(4*bat$N1))
  tse.delta.M = sum(((bat$X2 - delta.M)^2 - 1/(4*bat$N2))[ind])
  tse.M = tse.delta.M/tse.zero
  
  # XKB theta.hat.SG
  delta.SG = thetahat.SG(bat$X1,1/(4*bat$N1))
  tse.delta.SG = sum(((bat$X2 - delta.SG)^2 - 1/(4*bat$N2))[ind])
  tse.SG = tse.delta.SG/tse.zero
  
  # group_linear: num bins = n^1/3
  delta.gl = grouplinear(x=bat$X1, v=1/(4*bat$N1))
  tse.delta.gl = sum(((bat$X2 - delta.gl)^2 - 1/(4*bat$N2))[ind])
  tse.gl = tse.delta.gl/tse.zero
  
  # group_linear: sure(equal-bins)
  delta.gl.sure = grouplinear.sure(bat$X1, 1/(4*bat$N1), kmax=min(ceiling(n^(1/3)/0.8),n))
  tse.delta.gl.sure = sum(((bat$X2 - delta.gl.sure)^2 - 1/(4*bat$N2))[ind])
  tse.gl.sure = tse.delta.gl.sure/tse.zero
  
  # dynamic_group_linear_all_division
  delta.dynamic = GroupSure(bat$X1,1/(4*bat$N1))
  tse.delta.dynamic = sum(((bat$X2 - delta.dynamic)^2 - 1/(4*bat$N2))[ind])
  tse.gl.dynamic = tse.delta.dynamic/tse.zero
  
  # dynamic_group_linear_with_minimum_bin_size_constraint
  delta.dynamicMin = GroupSureMin(bat$X1,1/(4*bat$N1),n^(2/3)*0.8)
  tse.delta.dynamicMin = sum(((bat$X2 - delta.dynamicMin)^2 - 1/(4*bat$N2))[ind])
  tse.gl.dynamicMin = tse.delta.dynamicMin/tse.zero
  
  delta.dynamicMin2 = GroupSureMin(bat$X1,1/(4*bat$N1),n^(2/3))
  tse.delta.dynamicMin2 = sum(((bat$X2 - delta.dynamicMin2)^2 - 1/(4*bat$N2))[ind])
  tse.gl.dynamicMin2 = tse.delta.dynamicMin2/tse.zero
  
  delta.dynamicMin3 = GroupSureMin(bat$X1,1/(4*bat$N1),n^(2/3)*1.2)
  tse.delta.dynamicMin3 = sum(((bat$X2 - delta.dynamicMin3)^2 - 1/(4*bat$N2))[ind])
  tse.gl.dynamicMin3 = tse.delta.dynamicMin3/tse.zero
  
  # DPMM NG mixture
  delta.ngmixture = DPMM.ng.sknown(x.vec=bat$X1, precd=4*bat$N1)$mu.est
  tse.delta.ngmixture = sum(((bat$X2 - delta.ngmixture)^2 - 1/(4*bat$N2))[ind])
  tse.DPMM.ngmixture = tse.delta.ngmixture/tse.zero
  
  ####################pitchers only
  bat_p = bat[bat$Pitcher.==1,]
  n_p = dim(bat_p)[1]
  ind_p = bat_p$N2 > 10  # indicator for records with N2>=11 (among those with N1>=11)
  
  tse.zero_p = sum(((bat_p$X2 - bat_p$X1)^2 - 1/(4*bat_p$N2))[ind_p])
  
  # grand mean
  tse.delta.gm_p = sum(((bat_p$X2 - mean(bat_p$X1))^2 - 1/(4*bat_p$N2))[ind_p])
  tse.gm_p = tse.delta.gm_p/tse.zero_p
  
  # James-Stein
  delta.JS_p = JS(bat_p$X1,1/(4*bat_p$N1))
  tse.delta.JS_p = sum(((bat_p$X2 - delta.JS_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.JS_p = tse.delta.JS_p/tse.zero_p
  
  # XKB theta.hat.M
  delta.M_p = thetahat.M(bat_p$X1,1/(4*bat_p$N1))
  tse.delta.M_p = sum(((bat_p$X2 - delta.M_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.M_p = tse.delta.M_p/tse.zero_p
  
  # XKB theta.hat.SG
  delta.SG_p = thetahat.SG(bat_p$X1,1/(4*bat_p$N1))
  tse.delta.SG_p = sum(((bat_p$X2 - delta.SG_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.SG_p = tse.delta.SG_p/tse.zero_p
  
  # group_linear: num bins = n^1/3
  delta.gl_p = grouplinear(x=bat_p$X1, v=1/(4*bat_p$N1))
  tse.delta.gl_p = sum(((bat_p$X2 - delta.gl_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.gl_p = tse.delta.gl_p/tse.zero_p
  
  
  # group_linear: sure(equal-bins)
  delta.gl.sure_p = grouplinear.sure(bat_p$X1, 1/(4*bat_p$N1), kmax=min(ceiling(n_p^(1/3)/0.8),n_p))
  tse.delta.gl.sure_p = sum(((bat_p$X2 - delta.gl.sure_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.gl.sure_p = tse.delta.gl.sure_p/tse.zero_p
  
  # dynamic_group_linear_all_division
  delta.dynamic_p = GroupSure(bat_p$X1,1/(4*bat_p$N1))
  tse.delta.dynamic_p = sum(((bat_p$X2 - delta.dynamic_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.gl.dynamic_p = tse.delta.dynamic_p/tse.zero_p
  
  # dynamic_group_linear_with_minimum_bin_size_constraint
  delta.dynamicMin_p = GroupSureMin(bat_p$X1,1/(4*bat_p$N1),n_p^(2/3)*0.8)
  tse.delta.dynamicMin_p = sum(((bat_p$X2 - delta.dynamicMin_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.gl.dynamicMin_p = tse.delta.dynamicMin_p/tse.zero_p
  
  delta.dynamicMin2_p = GroupSureMin(bat_p$X1,1/(4*bat_p$N1),n_p^(2/3))
  tse.delta.dynamicMin2_p = sum(((bat_p$X2 - delta.dynamicMin2_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.gl.dynamicMin2_p = tse.delta.dynamicMin2_p/tse.zero_p
  
  delta.dynamicMin3_p = GroupSureMin(bat_p$X1,1/(4*bat_p$N1),n_p^(2/3)*1.2)
  tse.delta.dynamicMin3_p = sum(((bat_p$X2 - delta.dynamicMin3_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.gl.dynamicMin3_p = tse.delta.dynamicMin3_p/tse.zero_p
  
  # DPMM NG mixture
  delta.ngmixture_p = DPMM.ng.sknown(x.vec=bat_p$X1, precd=4*bat_p$N1)$mu.est
  tse.delta.ngmixture_p = sum(((bat_p$X2 - delta.ngmixture_p)^2 - 1/(4*bat_p$N2))[ind_p])
  tse.DPMM.ngmixture_p = tse.delta.ngmixture_p/tse.zero_p
  
  
  
  #########################Nonpitchers only
  bat_n = bat[bat$Pitcher.==0,]
  n_n = dim(bat_n)[1]
  ind_n = bat_n$N2 > 10  # indicator for records with N2>=11 (among those with N1>=11)
  
  tse.zero_n = sum(((bat_n$X2 - bat_n$X1)^2 - 1/(4*bat_n$N2))[ind_n])
  
  # grand mean
  tse.delta.gm_n = sum(((bat_n$X2 - mean(bat_n$X1))^2 - 1/(4*bat_n$N2))[ind_n])
  tse.gm_n = tse.delta.gm_n/tse.zero_n
  
  # James-Stein
  delta.JS_n = JS(bat_n$X1,1/(4*bat_n$N1))
  tse.delta.JS_n = sum(((bat_n$X2 - delta.JS_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.JS_n = tse.delta.JS_n/tse.zero_n
  
  # XKB theta.hat.M
  delta.M_n = thetahat.M(bat_n$X1,1/(4*bat_n$N1))
  tse.delta.M_n = sum(((bat_n$X2 - delta.M_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.M_n = tse.delta.M_n/tse.zero_n
  
  # XKB theta.hat.SG
  delta.SG_n = thetahat.SG(bat_n$X1,1/(4*bat_n$N1))
  tse.delta.SG_n = sum(((bat_n$X2 - delta.SG_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.SG_n = tse.delta.SG_n/tse.zero_n
  
  # group_linear: num bins = n^1/3
  delta.gl_n = grouplinear(x=bat_n$X1, v=1/(4*bat_n$N1))
  tse.delta.gl_n = sum(((bat_n$X2 - delta.gl_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.gl_n = tse.delta.gl_n/tse.zero_n
  
  # group_linear: sure(equal-bins)
  delta.gl.sure_n = grouplinear.sure(bat_n$X1,1/(4*bat_n$N1),kmax=min(ceiling(n_n^(1/3)/0.8),n_n))
  tse.delta.gl.sure_n = sum(((bat_n$X2 - delta.gl.sure_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.gl.sure_n = tse.delta.gl.sure_n/tse.zero_n
  
  # dynamic_group_linear_all_division
  delta.dynamic_n = GroupSure(bat_n$X1,1/(4*bat_n$N1))
  tse.delta.dynamic_n = sum(((bat_n$X2 - delta.dynamic_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.gl.dynamic_n = tse.delta.dynamic_n/tse.zero_n
  
  # dynamic_group_linear_with_minimum_bin_size_constraint
  delta.dynamicMin_n = GroupSureMin(bat_n$X1,1/(4*bat_n$N1),n_n^(2/3)*0.8)
  tse.delta.dynamicMin_n = sum(((bat_n$X2 - delta.dynamicMin_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.gl.dynamicMin_n = tse.delta.dynamicMin_n/tse.zero_n
  
  delta.dynamicMin2_n = GroupSureMin(bat_n$X1,1/(4*bat_n$N1),n_n^(2/3))
  tse.delta.dynamicMin2_n = sum(((bat_n$X2 - delta.dynamicMin2_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.gl.dynamicMin2_n = tse.delta.dynamicMin2_n/tse.zero_n
  
  delta.dynamicMin3_n = GroupSureMin(bat_n$X1,1/(4*bat_n$N1),n_n^(2/3)*1.2)
  tse.delta.dynamicMin3_n = sum(((bat_n$X2 - delta.dynamicMin3_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.gl.dynamicMin3_n = tse.delta.dynamicMin3_n/tse.zero_n
  
  # DPMM NG mixture
  delta.ngmixture_n = DPMM.ng.sknown(x.vec=bat_n$X1, precd=4*bat_n$N1)$mu.est
  tse.delta.ngmixture_n = sum(((bat_n$X2 - delta.ngmixture_n)^2 - 1/(4*bat_n$N2))[ind_n])
  tse.DPMM.ngmixture_n = tse.delta.ngmixture_n/tse.zero_n

    if(j%%50==0){
    print(j)
  }
  
  obj1 = c(tse.gm,tse.JS,tse.M,tse.SG,tse.gl,tse.gl.sure,tse.gl.dynamic,tse.gl.dynamicMin,
          tse.gl.dynamicMin2,tse.gl.dynamicMin3,tse.DPMM.ngmixture,
          tse.gm_p, tse.JS_p, tse.M_p, tse.SG_p, tse.gl_p, tse.gl.sure_p, tse.gl.dynamic_p,
          tse.gl.dynamicMin_p, tse.gl.dynamicMin2_p, tse.gl.dynamicMin3_p,tse.DPMM.ngmixture_p,
          tse.gm_n,tse.JS_n,tse.M_n, tse.SG_n, tse.gl_n, tse.gl.sure_n, tse.gl.dynamic_n,
          tse.gl.dynamicMin_n, tse.gl.dynamicMin2_n, tse.gl.dynamicMin3_n,tse.DPMM.ngmixture_n)
  
  obj1
}

save(list = ls(), file = "shuffle_par.Rdata")

listnames = c("tse.gm","tse.JS","tse.M","tse.SG","tse.gl","tse.gl.sure","tse.gl.dynamic",
              "tse.gl.dynamicMin",
         "tse.gl.dynamicMin2","tse.gl.dynamicMin3","tse.DPMM.ngmixture",
         "tse.gm_p", "tse.JS_p", "tse.M_p", "tse.SG_p", "tse.gl_p", "tse.gl.sure_p", "tse.gl.dynamic_p",
         "tse.gl.dynamicMin_p", "tse.gl.dynamicMin2_p", "tse.gl.dynamicMin3_p","tse.DPMM.ngmixture_p",
         "tse.gm_n","tse.JS_n","tse.M_n", "tse.SG_n", "tse.gl_n", "tse.gl.sure_n", "tse.gl.dynamic_n",
         "tse.gl.dynamicMin_n", "tse.gl.dynamicMin2_n", "tse.gl.dynamicMin3_n","tse.DPMM.ngmixture_n")

tse.average = apply(alloutput,1,mean)

average = cbind(listnames,tse.average)

print(average)

allresult = cbind(listnames,alloutput)

error = t(as.matrix(alloutput))

sd=cbind(listnames,sqrt(diag(cov(error))))

write.table(average, 'baseball_average.txt', sep="\t")
write.table(allresult, 'baseball_allresult.txt', sep="\t")
write.table(sd, 'baseball_sd.txt', sep="\t")
