%chk=XUBOJATI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4805949201250146  -1.3740228100362821  2.634514549543798e-16 
N  1.4805949201250148  1.3740228100362821  -3.9031386855166974e-16 
N  1.5116229275351682  2.681785422763614  0.33526858115358105 
C  2.773058018717449  3.1122010284436863  0.3891061806640282 
C  3.5907163566922446  2.069042807538932  0.10040510948635337 
C  2.7549548229246272  0.9864159296338109  -0.13639092492265115 
C  2.9953507004226863  -0.4517702457091  -0.48883795233878763 
C  1.4159352549298587  -2.893107045577381  -0.9663324236655729 
C  0.9287089270223423  -4.055754749335529  -0.3895495942230946 
C  0.8428034661439094  -5.213019038097963  -1.1270593960199569 
C  1.2110034513448849  -5.217650646069474  -2.438791460545953 
C  1.6746495737870388  -4.076431261975992  -3.023331288542018 
C  1.7796510927445917  -2.9020300524666753  -2.2909902886521603 
C  1.7257685129100555  -1.881272813390735  1.7104749262990786 
C  0.8075256409568976  -1.5063211122713074  2.6860856119390473 
C  1.0188354537391437  -1.8604686755036464  4.004630190813271 
C  2.120289930072075  -2.6033209466191796  4.356564389016095 
C  3.0265120634952707  -2.976290060564578  3.3883029341720268 
C  2.836512367212616  -2.6196060620019503  2.0836311859602685 
H  0.8189485179347098  3.1666082729624474  0.49093908641951134 
H  3.0447471607768835  3.9783058604068295  0.5895320276426725 
H  4.5197721872317835  2.0778575811449684  0.06841552513832003 
H  3.767398788221503  -0.7934436091866708  -0.010702822050164276 
H  3.154736548347666  -0.5439770664944519  -1.4410894914196792 
H  0.6596124501023952  -4.054254439453521  0.5004030141665627 
H  0.5319304830691367  -5.993904735309693  -0.7293917756579941 
H  1.1467806172513062  -6.0008607442858075  -2.936065760739483 
H  1.922577758369094  -4.084586082928131  -3.9196526395010873 
H  2.0937260903179435  -2.1256754117661827  -2.6945804616768414 
H  0.051722918741606394  -1.0171118028696318  2.4508954193571713 
H  0.41338276978358723  -1.593575897362418  4.657677227761712 
H  2.2523054866254  -2.853081752985224  5.242321992670367 
H  3.7741414666578157  -3.474325784572633  3.6254402427900705 
H  3.458876417374194  -2.874348026577311  1.440609288548695 
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
