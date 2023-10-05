#!/bin/sh
# Script Syntax is provided below
if [ $# -eq 0 ]
 then
  echo "
%*****************************************************************************************
Script usage: submit-g16.sh partition(compute or amd) jobname(without .com) #procs #hours
%*****************************************************************************************
"
  exit
 else
 echo "#!/bin/sh



#SBATCH -J CERUHELE_PN_Pd                     #Job name
#SBATCH -p amd                         #Queue name (compute, debug)
#SBATCH -N 1                                    #Request 1 node only
#SBATCH -n 32                     #Number of cores
#SBATCH -t 72:00:00                #Requested time (e.g. 1-03:00:00 equals 1 day and 3 hours)
#SBATCH -w dalton-n06     #This ensure we select only the node that we actually want

date

#Load up the g16 module
. /etc/profile.d/modules.sh

module load apps gaussian/g16

export GAUSS_SCRDIR=/home/shared/scratch/

#Run the g16 code
g16 < CERUHELE_PN_Pd_gaussian.com > CERUHELE_PN_Pd_gaussian.log


rm fort.7

date"> run.sh
  chmod u+x CERUHELE_PN_Pd.com
  chmod u+x run.sh
  sbatch run.sh
fi
exit

