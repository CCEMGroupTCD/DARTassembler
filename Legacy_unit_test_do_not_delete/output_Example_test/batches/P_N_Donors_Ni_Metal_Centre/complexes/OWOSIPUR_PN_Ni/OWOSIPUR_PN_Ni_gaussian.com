%chk=OWOSIPUR_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3436185312431497  -1.508240445850727  -2.5452081209887793e-15 
N  1.3436185312431501  1.5082404458507275  -1.8470618343556704e-16 
H  0.9967571717849236  2.310127869256756  -0.2889518082622693 
H  1.4025579333351355  1.5719537949154367  0.9157851746081156 
C  1.132095141290205  -3.031314725100864  -0.940190339978833 
H  1.0059837367278026  -2.8156836525522446  -1.886932197187463 
H  0.34587912020022626  -3.5139221349488547  -0.6068992417759589 
H  1.9284151856121228  -3.5928733925733205  -0.8356701983902195 
C  1.5586938951768188  -2.0601112614689283  1.7072889886301514 
H  1.6827411098863365  -1.2818319968932474  2.290088555251772 
H  2.3466022744186272  -2.640100723630969  1.7660428795782483 
H  0.7644560537679055  -2.5592995156171763  1.9919284064811398 
C  3.0174690551136063  -1.0154243852253624  -0.5120113906236343 
H  3.6499993078598925  -1.726323103629124  -0.23847804433359548 
H  3.040377379574555  -0.963016166617078  -1.5002304742749868 
C  3.5082554207751637  0.3035659624163891  0.04802900623058778 
H  3.431071809019513  0.28464099585782454  1.0343812687195704 
H  4.466389284666198  0.41453281485721316  -0.17843062929992515 
C  2.7376283042822895  1.4868817559358514  -0.4825364846688089 
H  3.1908365564619383  2.3211825953422878  -0.2027736395267261 
H  2.738101930008514  1.4571642462885617  -1.4722751245627193 
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