%chk=ZUMUVAMI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  2.554762653370124  -0.3030658105118836  -0.43639473020589187 
C  2.0714779495082514  -2.5354505986560025  -0.1154791763978993 
H  1.48114226566289  -3.213378850285834  0.192366037989078 
C  3.2866023284711803  -2.906674905609529  -0.6735929171146599 
H  3.524615949893137  -3.823406044734971  -0.7424957078809202 
C  4.136742850362993  -1.931306079682029  -1.1231404634835915 
H  4.974544760141914  -2.1634059746867935  -1.5084697365859445 
C  3.761010423682917  -0.5965534326720109  -1.006121982717302 
H  4.3310404587149565  0.09721462128044439  -1.3176225067560337 
C  1.8512392117869236  2.326753105055148  -1.4537143694878674 
C  1.6199477278582146  1.7565530347554874  -2.7040113859172967 
H  1.4248061017401126  0.8293293041564544  -2.773459202081614 
C  1.676124135678565  2.541883565984092  -3.8458265626970713 
H  1.5164301814241554  2.1528356799838018  -4.697741325671226 
C  1.960818171000631  3.8838260545634355  -3.743518774149375 
H  2.016385343544881  4.417129377559548  -4.529952360039982 
C  2.1657770644506797  4.455814222620711  -2.5135364663098434 
H  2.3423446864871504  5.385932460939117  -2.4467188967951605 
C  2.1134289390002605  3.68019992039294  -1.3698480174997953 
H  2.259136550966059  4.081754288803508  -0.5221676526064988 
C  2.726532766107498  1.9900118306288725  1.3112254631057454 
C  3.9744736541833774  2.568703557948829  1.0524849084958519 
H  4.290156334246633  2.6527987766784196  0.16047840084306542 
C  4.745427454935171  3.0199150936354684  2.1170908577417076 
H  5.5946493493901  3.412532266876113  1.9514077286821505 
C  4.289566763634576  2.9020686716463056  3.415608228034217 
H  4.820275898959773  3.2204313031918588  4.135894079983147 
C  3.0722144949287946  2.3283989113749266  3.6641759653355117 
H  2.76709931257872  2.238054070060368  4.559223888715806 
C  2.2796776388152633  1.8803004472732752  2.6225771998929552 
H  1.428881074853675  1.4966103584686499  2.80188925332011 
N  1.703126712109231  -1.2615305792964358  1.5449293859618863e-16 
P  1.703126712109231  1.2615305792964358  -1.5449293859618863e-16 
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

-Pd 0
lanl2dz
