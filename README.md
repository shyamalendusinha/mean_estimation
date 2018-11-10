# Normal Mean Estimation
---
Store the "mean_estimation" folder in your home folder and run it. In each iteration, we ran 40 simulations
in parallel.

**DPMMfunction.R**: R function for estimating mixture of normal-inverse gamma and all other variantions.
(Our Algorithm)

**SUREmethods.R**: R function for estimating means for several other SURE algoritms.
(SURE Estimates for a Heteroscedastic Hierarchical
Model-Xianchao Xie , S. C. Kou & Lawrence D. Brown, On SURE-Type Double Shrinkage Estimation-
Bing-Yi Jing, Zhouping Li, Guangming Pan & Wang Zhou)

**grouplinearfunction_all.R**: R function for estimating means for group linear algorithms.
(Group-Linear Empirical Bayes Estimates for a Heteroscedastic
Normal Mean-Asaf Weinstein, Zhuang Ma, Lawrence D. Brown, Cun-Hui Zhang)

**bash_script.sh**: Bash script to execute all codes below.

**diffq_sigma2known_exam1.R - diffq_sigma2known_exam8.R**: 6 same example of Xie(2012) plus two
more examples.
q=20,60,...,500, B<sub>q</sub>=1000, X<sub>i</sub> = &mu;<sub>i</sub> + &sigma;<sub>i</sub> &epsilon;<sub>i</sub>,
&sigma;<sub>i</sub> known.

**diffq_exam9.R - diffq_exam14.R**: 8 examples of R codes where, n=4, q=20,60,...,500, B<sub>q</sub>=1000, 
X<sub>ij</sub> = &mu;<sub>i</sub> + &sigma;<sub>i</sub> &epsilon;<sub>ij</sub>,
&sigma;<sub>i</sub> unknown.

**baseball_DPMM.R**: Brown baseball data used to estimate loss for mixture of normal-inverse gamma method.

**baseball_shuffle.R**: 1000 hypergeometric shuffles on 
Brown baseball data used to estimate loss for several algorithms.

**prostatedata_500_cntl.R**: Randomly select 3 columns and 500 rows from control group and tried to estimate overall mean from 3 columns. B=100, q=500, n=3.
