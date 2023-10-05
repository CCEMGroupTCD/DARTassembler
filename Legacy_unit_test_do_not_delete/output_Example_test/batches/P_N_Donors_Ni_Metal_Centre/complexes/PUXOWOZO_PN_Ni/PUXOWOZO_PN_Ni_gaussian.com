%chk=PUXOWOZO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4805949201250146  1.3740228100362821  -9.51821913277413e-17 
N  1.4805949201250148  -1.3740228100362821  2.220446049250313e-16 
N  1.5116229275351682  -2.681785422763614  -0.3352685811535814 
C  2.773058018717449  -3.1122010284436863  -0.3891061806640286 
C  3.5907163566922446  -2.069042807538932  -0.10040510948635362 
C  2.7549548229246272  -0.9864159296338109  0.13639092492265104 
C  2.9953507004226863  0.45177024570909996  0.4888379523387877 
C  1.4159352549298587  2.893107045577381  0.9663324236655733 
C  0.9287089270223423  4.055754749335529  0.3895495942230951 
C  0.8428034661439094  5.213019038097963  1.1270593960199575 
C  1.2110034513448849  5.217650646069474  2.4387914605459535 
C  1.6746495737870388  4.076431261975992  3.0233312885420185 
C  1.7796510927445917  2.902030052466675  2.2909902886521607 
C  1.7257685129100555  1.8812728133907353  -1.7104749262990784 
C  0.8075256409568976  1.5063211122713076  -2.6860856119390473 
C  1.0188354537391437  1.8604686755036468  -4.004630190813271 
C  2.120289930072075  2.60332094661918  -4.356564389016095 
C  3.0265120634952707  2.9762900605645783  -3.3883029341720263 
C  2.836512367212616  2.619606062001951  -2.083631185960268 
H  0.8189485179347098  -3.1666082729624474  -0.4909390864195117 
H  3.0447471607768835  -3.9783058604068295  -0.589532027642673 
H  4.5197721872317835  -2.0778575811449684  -0.06841552513832028 
H  3.767398788221503  0.7934436091866708  0.010702822050164373 
H  3.154736548347666  0.5439770664944517  1.4410894914196792 
H  0.6596124501023952  4.054254439453521  -0.5004030141665623 
H  0.5319304830691367  5.993904735309693  0.7293917756579948 
H  1.1467806172513062  6.0008607442858075  2.936065760739484 
H  1.922577758369094  4.08458608292813  3.9196526395010878 
H  2.0937260903179435  2.1256754117661822  2.694580461676842 
H  0.051722918741606394  1.017111802869632  -2.4508954193571713 
H  0.41338276978358723  1.5935758973624186  -4.657677227761712 
H  2.2523054866254  2.8530817529852244  -5.242321992670367 
H  3.7741414666578157  3.4743257845726334  -3.62544024279007 
H  3.458876417374194  2.874348026577311  -1.4406092885486945 
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