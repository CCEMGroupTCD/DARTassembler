%chk=OFOWOREK_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.343618531243151  1.508240445850727  -2.6417132793760433e-15 
N  1.3436185312431501  -1.5082404458507275  0.0 
H  0.9001601341826995  -2.310127869256756  0.08456122196267599 
H  1.8391495565687992  -1.5719537949154365  -0.7723893720725905 
C  0.7028026841849596  3.031314725100864  0.7197604254512712 
H  0.13351360419502667  2.815683652552244  1.4866594891850609 
H  0.17674538882216395  3.5139221349488547  0.04709237213268713 
H  1.4499522582907018  3.5928733925733205  1.0144086681071278 
C  2.3594378092882478  2.0601112614689283  -1.3889579897249993 
H  2.750478784462503  1.2818319968932477  -1.838546692948627 
H  3.0770424567519377  2.640100723630969  -1.0583597384901926 
H  1.8027776691116373  2.5592995156171767  -2.0229633802848808 
C  2.5593731383758134  1.0154243852253624  1.259314089346042 
H  3.2452081605291614  1.726323103629124  1.326733185749236 
H  2.1003110917417094  0.9630161666170778  2.134736150961152 
C  3.260137537900227  -0.3035659624163891  1.0074296730216288 
H  3.67082429579694  -0.28464099585782443  0.10732718856273737 
H  3.9883504874342877  -0.41453281485721316  1.670008247826395 
C  2.3289086208008154  -1.4868817559358514  1.0978652675918936 
H  2.8609252083077044  -2.3211825953422878  1.072898890536961 
H  1.8494880490446826  -1.457164246288562  1.963739805268916 
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
-N 0
6-31+g(d)
****
-C 0
6-31g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****

-Ni 0
lanl2dz


