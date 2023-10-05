%chk=LELIGOZE_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
C  2.4082671143732837  -1.0539917766477398  -1.4953930590458213 
C  2.720165148295062  0.29863304798093115  -1.7729170404131784 
C  3.2674979998733855  0.624943068478184  -3.011829813818254 
H  3.488297492312125  1.5289892181190732  -3.2064774530494633 
C  3.4954490589934837  -0.36382695295291634  -3.967730195325222 
H  3.8185662710420942  -0.12432779514535652  -4.828766849996468 
C  3.2526388704065563  -1.6870034414288066  -3.666037947060568 
H  3.4554657843009022  -2.366013629318436  -4.299023103080972 
C  2.7077288120412426  -2.0246144428346975  -2.425367930488042 
H  2.540518943427334  -2.9375198035823074  -2.2202347489780676 
C  2.580382751297615  1.3468853449757798  -0.7355005506007716 
C  1.4277322616857808  2.299586673202209  1.0413868847741552 
C  3.5263101068914  3.040920699929501  0.6147706421994825 
H  4.2552718043913895  3.6165519469515264  0.8169590964440835 
C  3.646905279094711  2.181007833902723  -0.4663808352923983 
H  4.432107881812289  2.1671801900859684  -1.0016792464165645 
C  2.446374735550421  -1.0572948719051383  1.45963046518713 
C  1.7950151726911807  -0.6475434809962745  2.6141332748502872 
H  0.8496540513407393  -0.5580915088796742  2.623050388929679 
C  2.5336716008453037  -0.3676420279221557  3.759039815889172 
H  2.095931992597662  -0.0849808379168564  4.553136665798573 
C  3.9021611431909333  -0.503867257634475  3.7337013218431094 
H  4.410264139917982  -0.2992410417958813  4.510547971130186 
C  4.549198377610433  -0.9368692998083468  2.5849330298241906 
H  5.491794055989048  -1.0464223358918403  2.587595620083346 
C  3.8397188925788828  -1.2052753667606348  1.4486248658375969 
H  4.286478378802816  -1.4899521631747414  0.6593414483934259 
C  1.129556721135184  -3.1926432192363547  0.042437213558602825 
C  1.6783321717904034  -3.978475071187406  1.0625065127816387 
H  2.2247496741063784  -3.5770588568523443  1.7282825306973464 
C  1.425502824700325  -5.339271844169938  1.1019688225256987 
H  1.8075105924802566  -5.873343384538782  1.7897916086587624 
C  0.6212275071175699  -5.9212402023950235  0.1475929251183967 
H  0.4376261027114374  -6.851862701347842  0.1858865113034185 
C  0.07494604082287593  -5.145658991396184  -0.875307346036382 
H  -0.47372734183027854  -5.553893331885662  -1.5356350468131108 
C  0.3250734237457116  -3.7929069013896846  -0.9342242815916457 
H  -0.0469348039555153  -3.2690133759206974  -1.6339511515161143 
N  1.4429130439842859  1.4135423401865257  3.710105234956919e-15 
N  2.4517894978209798  3.1000242779581413  1.3767537704557309 
N  0.33181564084046133  2.359198177296866  1.7951686295349243 
H  -0.4053203686382527  1.9348192903362929  1.4954835384277279 
H  0.3716286111040059  2.814844362550156  2.543101682196017 
P  1.4429130439842863  -1.4135423401865266  -1.1102230246251565e-15 
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
