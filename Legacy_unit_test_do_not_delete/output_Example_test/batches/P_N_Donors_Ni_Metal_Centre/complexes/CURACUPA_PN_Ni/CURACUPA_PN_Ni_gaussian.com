%chk=CURACUPA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4437283851195835  -1.4127095773725038  -2.3882467389131803e-15 
N  1.4437283851195832  1.4127095773725038  -3.95051631130436e-16 
N  1.4672812057301354  2.7877561977917757  -0.00545034997281613 
H  0.7610276042705514  3.2382794760325497  -0.011241393796008936 
C  2.6965559666106733  3.1662356096709336  -0.2008333764765066 
C  3.532717111857112  2.009794571705762  -0.36878543776240363 
H  4.46562089653401  1.9967312639349428  -0.5452202940850976 
C  2.7120234936351704  0.9435014112943485  -0.22411145575593125 
C  3.025468910821683  4.627943446210579  -0.21581255299041122 
H  2.210273286758393  5.147651942393871  -0.05832266917886521 
H  3.4039145467239744  4.868839370929198  -1.0860687187947602 
H  3.6781712102349458  4.823254396471975  0.4903278374764179 
C  2.994170387937011  -0.5579599321935906  -0.24414702530660562 
H  3.630343918874776  -0.7946976278649306  0.4774603416059113 
H  3.39379940324491  -0.817387231243721  -1.1111372770577892 
C  0.7082035633874187  -3.8722017445186374  -1.2275333870664045 
H  0.4273874185532782  -4.131616804676728  -0.3579129069926631 
C  0.5549002984967204  -4.792605195324397  -2.323174727432258 
H  0.18190652611991998  -5.646822544260653  -2.14237058880585 
C  0.9049464096445907  -4.513989642415029  -3.5753901536857753 
H  0.8023550998056534  -5.142244157154525  -4.281095657316727 
C  1.4194038475619732  -3.279221557070589  -3.778854794314847 
H  1.6792212079427877  -3.0199651919504524  -4.6565323127630815 
C  1.5801822885735186  -2.3379722151689593  -2.67618205534118 
H  1.9351810778566922  -1.4767557560416404  -2.862993109981123 
C  1.2507550093001591  -2.632414970185631  -1.41189751787 
C  0.9506840801734746  -2.2549888816673636  2.5091957729353953 
H  0.1933112142237967  -1.6812765822806601  2.4995947596768024 
C  1.259256487565514  -2.978491465900564  3.593960470924872 
H  0.7234348593445604  -2.9194990506618965  4.374721437978785 
C  2.3992267389508584  -3.8410784544211922  3.567911033477047 
H  2.6037086163359087  -4.3456003755690835  4.345545936823297 
C  3.2086292886504717  -3.97035293724849  2.473550544136564 
H  3.9582240666393433  -4.553379543135811  2.481455962488736 
C  2.906193246021074  -3.2462466691345693  1.3965105204620347 
H  3.4480387378257493  -3.3106938644097186  0.6185882929661664 
C  1.775677993138351  -2.375557534084086  1.4005395531696263 
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
-P 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-N 0
6-31+g(d)
****
-Br 0
6-31+g(d)
****
-C 0
6-31g(d)
****

-Ni 0
lanl2dz
