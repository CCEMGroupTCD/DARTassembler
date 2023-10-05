%chk=UYIFAXEQ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.4890637184150322  1.5082404458507273  -2.719585859971981e-15 
N  1.4890637184150313  -1.5082404458507273  1.1102230246251565e-16 
H  1.0456053213545842  -2.3101278692567555  0.08456122196267588 
H  1.9845947437406795  -1.5719537949154363  -0.7723893720725896 
C  0.8482478713568373  3.0313147251008634  0.7197604254512727 
H  0.27895879136690427  2.8156836525522446  1.4866594891850653 
H  0.322190575994044  3.513922134948855  0.04709237213268691 
H  1.5953974454625823  3.592873392573321  1.01440866810713 
C  2.504882996460129  2.060111261468928  -1.3889579897249917 
H  2.8959239716343848  1.281831996893248  -1.838546692948634 
H  3.22248764392382  2.6401007236309693  -1.058359738490192 
H  1.9482228562835142  2.559299515617177  -2.0229633802848803 
C  2.704818325547694  1.0154243852253622  1.2593140893460424 
H  3.390653347701039  1.7263231036291238  1.3267331857492373 
H  2.2457562789135883  0.963016166617078  2.1347361509611504 
C  3.405582725072112  -0.3035659624163891  1.0074296730216301 
H  3.8162694829688233  -0.2846409958578242  0.10732718856273793 
H  4.1337956746061675  -0.41453281485721316  1.6700082478263911 
C  2.474353807972696  -1.4868817559358511  1.0978652675918923 
H  3.0063703954795846  -2.3211825953422878  1.0728988905369623 
H  1.9949332362165655  -1.4571642462885617  1.9637398052689186 
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

-Pd 0
lanl2dz


