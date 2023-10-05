#!/bin/sh

#SBATCH -J SESOPUNE_PN_Ni                     #Job name
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
g16 < SESOPUNE_PN_Ni_gaussian.com > SESOPUNE_PN_Ni_gaussian.log


rm fort.7

  chmod u+x SESOPUNE_PN_Ni.com
  chmod u+x run.sh
fi
exit

