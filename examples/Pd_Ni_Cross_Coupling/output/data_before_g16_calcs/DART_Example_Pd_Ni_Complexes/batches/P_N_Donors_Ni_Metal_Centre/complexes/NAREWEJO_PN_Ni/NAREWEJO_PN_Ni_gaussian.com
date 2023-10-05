%chk=NAREWEJO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.4437283851195843  1.4127095773725042  2.2335634654821598e-15 
N  1.4437283851195832  -1.4127095773725038  2.220446049250313e-16 
N  1.4671885090545376  -2.7877561977917757  -0.005836459414371709 
H  0.8405410213018214  -3.2382794760325497  0.31995582388344496 
C  2.651182378862228  -3.1662356096709336  -0.3898279710932818 
C  3.472456054785094  -2.009794571705762  -0.6197908049298853 
H  4.383779161965646  -1.996731263934943  -0.8861156522665066 
C  2.6755320731991583  -0.9435014112943485  -0.3761091609083266 
C  2.951046361875477  -4.62794344621058  -0.5258047790055856 
H  2.1532028312897507  -5.147651942393872  -0.29603822422135667 
H  3.6833319241299725  -4.868839370929199  0.0777884190737117 
H  3.2120273402511295  -4.823254396471976  -1.4513011169962307 
C  2.936022755034461  0.5579599321935909  -0.4863493472668994 
H  3.175254632495608  0.7946976278649311  -1.4181229583198491 
H  3.6857001525906075  0.8173872312437214  0.10471684708106499 
C  1.3456594660949937  3.8722017445186383  1.4276615379242386 
H  0.7006510126179265  4.131616804676729  0.7803118785388986 
C  1.7064760165460042  4.7926051953243975  2.473483346200134 
H  1.292052770626663  5.646822544260654  2.48172130820483 
C  2.5868632924919273  4.51398964241503  3.4302978517494305 
H  2.815837360362768  5.142244157154525  4.105661539662961 
C  3.137619239972067  3.2792215570705894  3.378027384656736 
H  3.767575458379252  3.019965191950453  4.042089166447273 
C  2.780270932046359  2.337972215168959  2.322546895404364 
H  3.1813876134294974  1.4767557560416404  2.327830686125461 
C  1.912775907988649  2.632414970185632  1.345617979232221 
C  -0.1347283502006169  2.254988881667364  -2.0118723737550903 
H  -0.8051937461977483  1.6812765822806608  -1.6594777224011994 
C  -0.3522611893698595  2.9784914659005652  -3.1184937383106854 
H  -1.1841398175259015  2.919499050661898  -3.570898925020054 
C  0.6752859391200083  3.8410784544211936  -3.6128191837078862 
H  0.5044417675632967  4.345600375569085  -4.398529785727705 
C  1.8932981570479646  3.9703529372484914  -3.0051979159201316 
H  2.557603009905228  4.5533795431358115  -3.3525506031096777 
C  2.1127916085029526  3.24624666913457  -1.9082451378815606 
H  2.948748777622324  3.3106938644097195  -1.461104063477619 
C  1.1036659999175968  2.375557534084086  -1.3985918476907855 
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
-H 0
6-31g(d,p)
****
-P 0
6-31+g(d)
****
-N 0
6-31+g(d)
****
-C 0
6-31g(d)
****
-Br 0
6-31+g(d)
****

-Ni 0
lanl2dz

