%chk=UDEWENEW_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
C  1.2225846851497222  2.791253288052885  0.2538606567339347 
H  0.33945680669808476  3.14138329677687  0.2622996905932671 
C  2.2689812442947908  3.6530661300120357  0.5017143160647451 
H  2.1066785332349793  4.572094763072015  0.681827090732574 
C  3.558749915617491  3.154542881242927  0.4842502125594426 
H  4.297061611311233  3.713529152048335  0.6947593312524263 
C  3.7518652430256623  1.8237964263756903  0.15336722239777634 
H  4.633175805066008  1.4735455571017246  0.08448099874612491 
C  2.655178664704339  1.0012543516572896  -0.07736741793146606 
C  2.8416308716970473  -0.43469041835237476  -0.4817055462219427 
H  3.6564607914711047  -0.7951190990113656  -0.04990784127265546 
H  2.9684999060995203  -0.48487757522336683  -1.46272295554505 
C  1.8370772101864548  -2.0251084100267693  1.7583282934564064 
C  2.152584760850531  -0.7353986476824351  2.5438890131896015 
H  1.4121147613871154  -0.10062717762666162  2.4423205107003367 
H  2.976331073047593  -0.3335685501384118  2.1966676920563324 
H  2.269434788272735  -0.9505761491942997  3.493297136145448 
C  3.0742047721388115  -2.9298182990151735  1.8601717626848782 
H  3.3383415351359265  -3.017558115414031  2.799528321746294 
H  3.8114623032864667  -2.5347498687503576  1.3479610461099774 
H  2.862204501357004  -3.8150015155165837  1.4961997604475203 
C  0.6447616818336307  -2.7162549070276762  2.4418705237042517 
H  0.8589586108302898  -2.88302467024641  3.384417770097973 
H  0.45952563871563123  -3.5682871671463783  1.9935998844169807 
H  -0.14555067692627355  -2.138844022956377  2.385830390341014 
C  1.561790365264507  -2.936095413563182  -1.186966184588329 
C  3.007136076012559  -3.4579311495799314  -1.3583730136183232 
H  3.603632053301151  -2.706952991421332  -1.5589116259783173 
H  3.0360475218706187  -4.1047415179326805  -2.094376251172208 
H  3.297381641934523  -3.892820816320601  -0.5285429893212712 
C  0.6682758013814754  -4.090150121399107  -0.7230477052725303 
H  1.0422522014350706  -4.485321487012686  0.09137472442969524 
H  0.6231110675297171  -4.770636938635477  -1.4263652213231595 
H  -0.23272978914463116  -3.7509602390598458  -0.5386583629574602 
C  1.084454735932298  -2.4847975449601547  -2.580504652218291 
H  0.1504378559674291  -2.192609805560413  -2.5274679396016246 
H  1.1577760304087001  -3.2321728909682563  -3.209215762312266 
H  1.6419914988957833  -1.740309849397274  -2.8918086179811615 
N  1.3877186728584434  1.467765950347671  -8.646802187277018e-17 
P  1.3877186728584434  -1.4677659503476705  1.1102230246251565e-16 
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
