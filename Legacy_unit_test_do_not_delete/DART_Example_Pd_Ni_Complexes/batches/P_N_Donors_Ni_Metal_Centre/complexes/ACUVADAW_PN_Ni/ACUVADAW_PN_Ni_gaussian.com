%chk=ACUVADAW_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.4437283851195835  1.4127095773725038  2.561253765118585e-15 
N  1.4437283851195832  -1.4127095773725038  2.220446049250313e-16 
N  1.4672812057301354  -2.7877561977917757  0.005450349972815788 
H  0.7610276042705514  -3.2382794760325497  0.01124139379600854 
C  2.6965559666106733  -3.1662356096709336  0.2008333764765062 
C  3.532717111857112  -2.009794571705762  0.3687854377624034 
H  4.46562089653401  -1.9967312639349428  0.5452202940850974 
C  2.7120234936351704  -0.9435014112943485  0.22411145575593114 
C  3.025468910821683  -4.627943446210579  0.21581255299041066 
H  2.210273286758393  -5.147651942393871  0.058322669178864577 
H  3.4039145467239744  -4.868839370929198  1.0860687187947595 
H  3.6781712102349458  -4.823254396471975  -0.4903278374764185 
C  2.994170387937011  0.5579599321935906  0.24414702530660568 
H  3.630343918874776  0.7946976278649307  -0.47746034160591116 
H  3.39379940324491  0.8173872312437209  1.1111372770577892 
C  0.7082035633874187  3.8722017445186374  1.227533387066405 
H  0.4273874185532782  4.131616804676728  0.3579129069926636 
C  0.5549002984967204  4.792605195324397  2.3231747274322583 
H  0.18190652611991998  5.646822544260653  2.142370588805851 
C  0.9049464096445907  4.513989642415029  3.5753901536857757 
H  0.8023550998056534  5.142244157154524  4.281095657316728 
C  1.4194038475619732  3.2792215570705885  3.7788547943148476 
H  1.6792212079427877  3.019965191950452  4.6565323127630815 
C  1.5801822885735186  2.337972215168959  2.6761820553411804 
H  1.9351810778566922  1.47675575604164  2.862993109981123 
C  1.2507550093001591  2.632414970185631  1.4118975178700002 
C  0.9506840801734746  2.254988881667364  -2.509195772935395 
H  0.1933112142237967  1.6812765822806603  -2.4995947596768024 
C  1.259256487565514  2.9784914659005644  -3.5939604709248716 
H  0.7234348593445604  2.919499050661897  -4.374721437978785 
C  2.3992267389508584  3.8410784544211927  -3.5679110334770465 
H  2.6037086163359087  4.345600375569084  -4.345545936823296 
C  3.2086292886504717  3.9703529372484905  -2.4735505441365637 
H  3.9582240666393433  4.553379543135811  -2.4814559624887353 
C  2.906193246021074  3.2462466691345693  -1.3965105204620343 
H  3.4480387378257493  3.3106938644097186  -0.6185882929661659 
C  1.775677993138351  2.375557534084086  -1.400539553169626 
Br  -1.5839191898578664  1.5839191898578664  0.0 
C  -1.2937473829822488  -1.2398794889499034  -0.4598005675376139 
C  -1.601562203225019  -1.4938654694550033  -1.7919791134336869 
H  -1.1309417375171225  -1.0428217585005106  -2.480879663633267 
C  -2.60218398542082  -2.411889969420291  -2.121542410022959 
H  -2.8231586547736818  -2.5626138747145775  -3.03270050271447 
C  -3.2683208079706407  -3.09759503004543  -1.143064365570095 
H  -3.9596020908062854  -3.707191723557817  -1.3683280760515548 
C  -2.9347449854640226  -2.896087320001541  0.14198826503765638 
H  -3.3668350933519986  -3.405606940434718  0.8171022443470677 
C  -1.9675509241539992  -1.9516175659140573  0.5112421903797603 
H  -1.7755298364628027  -1.8022912010877177  1.4301324877763764 

-Ni 0
lanl2dz
F 1 1.0
3.130 1.0
****
-C 0
6-31g(d)
****
-N 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****
-P 0
6-31+g(d)
****

-Ni 0
lanl2dz


