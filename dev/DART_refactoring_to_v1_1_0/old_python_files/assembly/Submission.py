
import random
class SumbitCalc:
    def __init__(self, calc_name: str, queue: str, num_cores: str, requested_time:str, node_preference: list):
        self.calc_name = calc_name
        self.queue = queue
        self.num_cores = num_cores
        self.requested_time = requested_time
        self.node_preference = node_preference


    def gen_submission_string(self):
        sumbit_string = f"""#!/bin/sh

#SBATCH -J {self.calc_name}                     #Job name
#SBATCH -p {self.queue}                         #Queue name (compute, debug)
#SBATCH -N 1                                    #Request 1 node only
#SBATCH -n {self.num_cores}                     #Number of cores
#SBATCH -t {self.requested_time}                #Requested time (e.g. 1-03:00:00 equals 1 day and 3 hours)
#SBATCH -w dalton-n0{str(random.choice(self.node_preference))}     #This ensure we select only the node that we actually want

date

#Load up the g16 module
. /etc/profile.d/modules.sh

module load apps gaussian/g16

export GAUSS_SCRDIR=/home/shared/scratch/

#Run the g16 code
g16 < {self.calc_name}_gaussian.com > {self.calc_name}_gaussian.log


rm fort.7

  chmod u+x {self.calc_name}.com
  chmod u+x run.sh
fi
exit

"""
        return sumbit_string
