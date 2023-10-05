%chk=XUJIQILO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.442208277260951  1.4142613920347253  0.0 
N  2.7994598919942777  0.45371381387765686  -0.26882824330863136 
H  3.5348168308484422  0.7311840058875116  -0.3323904571236541 
N  1.442208277260951  -1.4142613920347253  0.0 
N  3.8183506974871397  -1.6081506310002949  -0.10715419644189761 
C  1.3182800232420229  2.5362338533751942  -1.409291965438492 
C  1.3771406702714204  2.023693438935173  -2.7031165044687784 
H  1.485054068469408  1.1122538335978103  -2.8288409822439498 
C  1.2158149687084814  2.864213516549518  -3.790370433936238 
H  1.265016389074601  2.5147123049744353  -4.637746243917574 
C  0.9835272833082686  4.209291088993984  -3.6035728019429163 
H  0.9122917299032938  4.707993729094537  -4.326444590497315 
C  0.9251178415819445  4.721888702902544  -2.3357285081384656 
H  0.8077294450001127  5.67538257793501  -2.198584187718261 
C  1.0830462907651768  3.894588928907673  -1.2367659387637655 
H  1.0016533479169958  4.272882282606793  -0.42438510602483165 
C  1.806629344275919  2.4463328027518183  1.4328507037008962 
C  3.042967882164147  3.088640794835997  1.504095303757413 
H  3.593038286605782  2.9580817722677937  0.8190727772785673 
C  3.3463071679857483  3.8811074324124837  2.5932227049394196 
H  4.188038082283534  4.265653975974387  2.633788500260385 
C  2.4289783817206683  4.043370063841487  3.6093562130628403 
H  2.6018151268795746  4.62165288666772  4.313966279941672 
C  1.2099809351288024  3.408686659447203  3.5451835445112163 
H  0.5657261809029049  3.5204618842326063  4.228118681095603 
C  0.8947842473030124  2.593768115045827  2.4693137410617414 
H  0.08779061160373525  2.201196382460224  2.407963869704001 
C  2.695136730388394  -0.9021860189256093  -0.11948507688771146 
C  1.3636160354553586  -2.763510512987287  0.11077079550455007 
H  0.4873145957265055  -3.127981959286475  0.15616684777491854 
C  2.4698681023716325  -3.55933709521528  0.10870791819340578 
H  2.410861517134454  -4.445366429719824  0.17286110436587487 
C  3.7245015479257924  -2.9429410139040297  0.02016929249205356 
C  5.001686493547853  -3.690700986675762  0.07438937914630056 
C  6.212408995055484  -3.0100579383527943  -0.010948274399921534 
H  6.189942847605692  -2.109374449307423  -0.11506731578920236 
C  7.417669447370674  -3.689872529499035  0.0322385725527312 
H  8.184729844271768  -3.2356205407095207  -0.025039767041686437 
C  7.428456605504415  -5.070251701973867  0.1646919130809796 
H  8.228774534993253  -5.464435830975609  0.19322091905488287 
C  6.232582558640525  -5.7583747253987365  0.2530334586305114 
H  6.25498002538563  -6.724065200081839  0.3317651586611993 
C  5.023615998540757  -5.086796585495247  0.21277332246209157 
H  4.2038107469059245  -5.582019313729635  0.29481269613141353 
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


