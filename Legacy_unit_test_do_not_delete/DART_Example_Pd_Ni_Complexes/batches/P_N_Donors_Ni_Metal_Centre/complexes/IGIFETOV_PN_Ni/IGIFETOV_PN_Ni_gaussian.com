%chk=IGIFETOV_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.4038951073709167  -1.4523011146108789  9.267213744023915e-17 
N  1.4038951073709167  1.4523011146108789  -6.683328867811896e-17 
C  -0.1163078442456964  3.24793519404813  -0.6810775056186599 
H  -0.701213766568807  3.150375880366846  0.09780515713960108 
H  -0.48963911396869086  2.7457972195858624  -1.435627343369569 
H  -0.04407732412193943  4.195466219918551  -0.9212461900377611 
C  1.2297996136004177  2.7173850056277025  -0.3563251853834416 
C  2.3401164804502916  3.602640826246192  -0.4946744436026072 
H  2.18179536176219  4.515201180753044  -0.7086518320334323 
C  3.614415430624657  3.1746127102128603  -0.328756700520521 
H  4.343648587740789  3.7800213208339284  -0.38409779509159303 
C  3.8282834220837785  1.810039279081821  -0.07048017568617437 
C  5.1162822767549105  1.2402654273129952  -0.014745911481755677 
H  5.874903201692674  1.810527827904695  -0.0646168119947278 
C  5.308269394767082  -0.10891076988698688  0.11008095031641386 
H  6.188337234431512  -0.46713372707268874  0.1553794109253021 
C  4.186609341263018  -0.9697235155516224  0.17013360450112336 
H  4.311285466926324  -1.9092116260221395  0.23789180676568866 
C  2.9144512991602367  -0.44595772597305494  0.13175061894331108 
C  2.7154647013205597  0.9516638275598317  0.038427165403974665 
C  1.5734045572913415  -2.7598608320106686  1.2414654976361987 
C  1.825805950797224  -2.360268730096259  2.5588638043923257 
H  2.007127491151281  -1.447235756658598  2.7469875360670133 
C  1.8137868304200673  -3.2824380010045413  3.587081483941692 
H  1.992019072252086  -3.0033201163442307  4.477353018044821 
C  1.543748602661183  -4.612358064115904  3.319202665662377 
H  1.5073373557966052  -5.241894174638662  4.029736307828359 
C  1.3235657437194313  -5.026519231485906  2.012628792484534 
H  1.1630261950049312  -5.944997167349365  1.8286537000656333 
C  1.3376722478179515  -4.102691027746298  0.9691176250706982 
H  1.1856420629177529  -4.388381547312669  0.07593702188390238 
C  1.4564394984518074  -2.2162594348503237  -1.6406587859463324 
C  2.4709170887567184  -3.080755187397064  -2.0413833579378626 
H  3.1459605265830772  -3.339908248847989  -1.4254285643140334 
C  2.493660141786843  -3.567125863487372  -3.3382100190527346 
H  3.1910963941937496  -4.148214982341427  -3.612630790456897 
C  1.496598688502418  -3.202982811692155  -4.241887525801492 
H  1.5175530179769454  -3.5308579886850344  -5.133261899669265 
C  0.48299673237978547  -2.367230740176869  -3.8373558465694315 
H  -0.2061821511239208  -2.1348186931128983  -4.44834301248386 
C  0.457377184316791  -1.862370905929305  -2.537346921015221 
H  -0.24024802692802227  -1.2778006529294863  -2.267594768209416 
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

