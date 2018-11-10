#!/bin/bash
# cd to mean_estimation folder in home directory and type: nohup ./bash_script.sh > nohup_bash_script.out&
# remember: chmod u+x bash_script.sh
cd ~/mean_estimation/simulation_examples
for i in {1..8}
{
echo "Code diffq_sigma2known_exam $i running"
time R CMD BATCH --slave --no-save --no-restore diffq_sigma2known_exam$i.R diffq_sigma2known_exam$i.Rout >/dev/null 2>&1
}
for i in {9..14}
{
echo "Code diffq_exam $i running"
time R CMD BATCH --slave --no-save --no-restore diffq_exam$i.R diffq_exam$i.Rout >/dev/null 2>&1
}
cd ~/mean_estimation/real_data_examples
echo "Code 15 running"
time R CMD BATCH --slave --no-save --no-restore baseball_DPMM.R baseball_DPMM.Rout >/dev/null 2>&1
echo "Code 16 running"
time R CMD BATCH --slave --no-save --no-restore baseball_shuffle.R baseball_shuffle.Rout >/dev/null 2>&1
echo "Code 17 running"
time R CMD BATCH --slave --no-save --no-restore prostatedata_500_cntl.R prostatedata_500_cntl.Rout >/dev/null 2>&1

