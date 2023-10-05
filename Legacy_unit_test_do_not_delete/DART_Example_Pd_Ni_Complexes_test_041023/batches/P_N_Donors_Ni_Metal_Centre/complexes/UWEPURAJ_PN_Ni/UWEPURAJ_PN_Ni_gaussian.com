%chk=UWEPURAJ_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.3583170046053312  -1.4950166938867273  -2.7454500943102698e-18 
C  1.6020287205971684  1.6301113994013023  1.438268994322789 
H  0.7482393268008826  1.7331520447866786  1.9298230467316362 
H  2.173621980695292  2.416946198510994  1.6237800018679516 
C  2.324879041473678  0.3110780246814002  1.857008615070304 
H  3.1362540598771913  0.5023701160469486  2.390463772784363 
H  1.7228722027265704  -0.27421293800687135  2.383374623753582 
C  2.692681506912785  -0.3466141027849034  0.4992390647099371 
H  3.597704625024946  -0.7720583353482285  0.4998478841915979 
C  2.586368549526072  0.8194572441963801  -0.4523508516095756 
H  3.3730001425883938  1.4163302782225866  -0.3775951290956125 
H  2.5014241574199936  0.511408332847639  -1.3892690809479737 
N  1.3583170046053312  1.4950166938867273  2.5080076233192965e-17 
C  1.097679731990093  2.7901812749363537  -0.7225977927715058 
H  0.9185634393149611  2.5830204410363504  -1.6743352316674418 
H  0.2740133814754422  3.1927363115863985  -0.34932219016991733 
C  3.231881922169904  3.8303116436195923  -1.5889208719831127 
H  3.320866701260559  3.1256573315254568  -2.2211694777841373 
C  4.145202493090052  4.879802988041248  -1.581203183331349 
H  4.861088908858569  4.8846602515284685  -2.2066619403326357 
C  4.021405256681252  5.921654096707206  -0.6692641000428372 
H  4.631304067167484  6.649082482254315  -0.6825224886256949 
C  3.0045021131149956  5.882115969647014  0.2522827642639486 
H  2.9380543355138204  6.5733689631592735  0.9019732227225975 
C  2.091935601542716  4.877174923741989  0.25539984817919853 
H  1.3815455853696663  4.889416355610166  0.8868581225469488 
C  2.184744435816641  3.8140511362290033  -0.665241652463034 
C  3.047347109819181  -2.1019133818702964  -2.1494517007358 
H  3.559239916826376  -1.3516581007175745  -1.872000034248181 
C  3.434605886865175  -2.803301673578194  -3.257398199379105 
H  4.2054794247813305  -2.5347321842595547  -3.744443523705118 
C  2.717491109104675  -3.887975748878429  -3.662735646315731 
H  2.97345210736286  -4.3571617037030554  -4.448123192159206 
C  1.5939062356680518  -4.317791200189159  -2.921250184557369 
H  1.1042076553034157  -5.089400955262169  -3.1818429439829456 
C  1.2244518299348996  -3.5800075423141173  -1.7981580038985905 
H  0.47555037859327154  -3.8517385113803266  -1.2794935009405057 
C  1.9378395262398826  -2.457233277074254  -1.433783520604591 
C  0.10625344033554707  -2.758160468616427  2.1528502166638943 
H  -0.6049990121507125  -2.1445973801566947  2.0134322289002995 
C  0.007044123051072715  -3.6815713413144127  3.197881819797077 
H  -0.7729881340707685  -3.726010570837182  3.7399007241988294 
C  1.0806485356688225  -4.534338447670616  3.4218558529451113 
H  1.0424814373261417  -5.16171840191149  4.133779841599003 
C  2.2007625664496553  -4.478354133018984  2.6225893863891 
H  2.9312148505036735  -5.063178756433017  2.7865676488780937 
C  2.2568176975763903  -3.581713884464297  1.5932249431785463 
H  3.0344594797019  -3.5491692278066456  1.048459265096935 
C  1.1970251032058274  -2.7116547872522663  1.3218361367109148 
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


