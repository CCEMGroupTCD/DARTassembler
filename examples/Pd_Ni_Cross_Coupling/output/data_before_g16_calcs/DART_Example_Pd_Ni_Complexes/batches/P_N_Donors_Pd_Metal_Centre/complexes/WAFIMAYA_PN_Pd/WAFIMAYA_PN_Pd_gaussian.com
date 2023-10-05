%chk=WAFIMAYA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
N  1.5214811706031712  1.4755321235066357  -2.0463224863963083e-15 
N  1.0826635849605468  2.646783984084804  -0.5274505622548636 
H  0.20058075408482545  2.7378648781290558  -0.6677839517756261 
C  2.1180320215255257  3.364242813435215  -0.9939819424986347 
H  2.063897717715223  4.211348416957767  -1.4206486433690777 
C  3.2521902930086153  2.6638941774840497  -0.7491472567454551 
H  4.140610633746094  2.926332547229416  -0.9554020754753432 
C  2.84926081317744  1.4754446330359487  -0.130297747911447 
C  3.6903572050035347  0.3579920891485928  0.2789490258288141 
C  5.046574936481842  0.6283054314647931  0.5268853519081285 
H  5.348137041685358  1.5290711442717397  0.5168077859011275 
C  5.948960072027661  -0.37624908238570454  0.7846895137220944 
H  6.859190931981394  -0.16508347979650484  0.9476404568644579 
C  5.529447983788225  -1.692760169506923  0.8052444476295257 
H  6.155083628231963  -2.391423094953807  0.9594963747083223 
C  4.181055102973304  -1.986805787415186  0.5987008842910413 
H  3.891691136282172  -2.891480839034959  0.6310713300330542 
C  3.2515725348432163  -0.9785153196588534  0.34590788958192287 
P  1.5214811706031728  -1.4755321235066352  -2.220446049250313e-16 
C  1.391115402524241  -3.0911292304911244  0.8118046405487905 
C  1.4969858178063848  -3.1067527398900348  2.1983868491184992 
H  1.5792474728367047  -2.2912094820853937  2.677939332962631 
C  1.4820874963224566  -4.320378943936517  2.8852638227234486 
H  1.5528842103942384  -4.3323549658193885  3.8329813809529734 
C  1.3626271119349032  -5.506760176925056  2.182520522097568 
H  1.367308224701352  -6.334739094182773  2.6484154241516658 
C  1.2364248151842687  -5.4922620411954455  0.8153454805044065 
H  1.1446330036932122  -6.310561053981777  0.3419697217079649 
C  1.2423161209594307  -4.285335067780764  0.11463760031441723 
H  1.1454413143130204  -4.279936164650774  -0.8302278364977466 
C  1.4613927963407898  -1.7061325046026221  -1.7978183700020907 
C  2.3182223533632698  -2.6076630296302046  -2.4381427273836778 
H  2.936104475204493  -3.1217128750496417  -1.9316108923197335 
C  2.2594117344392872  -2.746855180534995  -3.8186597161340297 
H  2.8370486764558516  -3.3613699346056283  -4.253952146859932 
C  1.3672572753880285  -2.0011766673588944  -4.565107755975044 
H  1.3260911853586175  -2.1157434004763394  -5.508027681893726 
C  0.5350657000608989  -1.0873502629594118  -3.94162594022026 
H  -0.06877699914718516  -0.5644055789994153  -4.45630126713142 
C  0.5852254778187038  -0.9355001836351371  -2.5594416796730046 
H  0.01713601239840945  -0.30268630237516536  -2.1334070907448517 
Br  -1.654629867976521  1.654629867976521  0.0 
C  -1.352040782177733  -1.3199206589279684  -0.4819701001757076 
C  -1.6380700204479186  -1.599217458033982  -1.8139382277462257 
H  -1.1579680913617836  -1.159339363685303  -2.503515938957912 
C  -2.63068216563725  -2.52632649780156  -2.142289956238523 
H  -2.8368222554424674  -2.694399528648736  -3.0538898079700485 
C  -3.3101544432602585  -3.1961573187854975  -1.1619874135498987 
H  -3.9960123560788676  -3.811990145619169  -1.3868663427944008 
C  -2.9974710358883443  -2.9701089012107795  0.12423560353005347 
H  -3.438659168087169  -3.468581175155742  0.8017084178296768 
C  -2.0389854250767776  -2.015940632791841  0.4912562833419403 
H  -1.8619130127668426  -1.8492165108145995  1.4101568590026907 

-Pd 0
lanl2dz
F 1 1.0
1.472 1.0
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

-Pd 0
lanl2dz

