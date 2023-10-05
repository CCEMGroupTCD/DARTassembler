%chk=AXALAWIF_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.433277970248619  -1.4233110201217436  -1.36708486844814e-14 
N  1.4332779702486187  1.423311020121744  1.3564704601237987e-14 
C  1.4404731912337627  2.724285202027264  0.2324187101122586 
H  0.9322392572137733  3.053855003264288  0.9644083715304869 
C  2.1545641433021547  3.6405169964585795  -0.5429246403120911 
H  2.1098475754632764  4.567002256080332  -0.34365865660106154 
C  2.9127449375372425  3.2080773881800297  -1.5765563244457494 
H  3.3765650511756635  3.8280457629941096  -2.1268459176734598 
C  3.822750937713743  1.2618345285111459  -2.8387764442468475 
H  4.312332355408744  1.832053903283209  -3.4187889614241262 
C  3.9165091308969733  -0.0792422586528633  -2.987539316432324 
H  4.472780224792606  -0.43810836225142946  -3.668061240975658 
C  3.2051807704103332  -0.9419845072830872  -2.1480173357815056 
H  3.2828531495665003  -1.8820317148748278  -2.2593872554869248 
C  2.391967003802179  -0.4336753676980493  -1.165967200394658 
C  2.2638960021941608  0.955409373478819  -1.008756989604421 
C  3.000037566250739  1.833097926675827  -1.8271990530733688 
C  2.5454224071942835  -2.556106064357513  0.8211191382045623 
C  3.567543590108269  -2.0157502397155422  1.6268996419485826 
H  3.702052898505599  -1.0760505036130597  1.664145741390093 
C  4.37160559191909  -2.8689998688064295  2.362698770583612 
H  5.063399385591369  -2.5045343003774563  2.9010318333488683 
C  4.19208179493363  -4.21497866256026  2.3308594382899517 
H  4.744434964519273  -4.782463350433273  2.8532698086678554 
C  3.1928262376436036  -4.756195507726048  1.526887748655408 
H  3.0733876595245793  -5.6974965312264505  1.4900007697332827 
C  2.3691116118287114  -3.9162935375313985  0.7757343251247587 
H  1.6836845231466349  -4.2863178909962825  0.23177054513227485 
C  0.32481656725370156  -2.4471514101690435  -0.9972136550301162 
H  0.8544654848805907  -3.0348225995120197  -1.576145902596981 
H  -0.2425991352352408  -1.8696527760828818  -1.5519012772082896 
H  -0.23897458842767927  -2.9895343321630667  -0.407791744422249 
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
