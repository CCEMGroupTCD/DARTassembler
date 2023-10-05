%chk=GUNUPUZO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5651067886888725  1.4291748458463718  1.004983847035852e-15 
N  1.5651067886888737  -1.4291748458463718  2.248201624865942e-15 
C  3.111518822097024  0.45129836413964886  0.12778270405413142 
H  3.792371243973326  0.8207385145414818  -0.4891502703412288 
H  3.4657511430510466  0.5046573847112994  1.0505541263008293 
C  2.827038066167768  -0.9780961749100516  -0.22595772513204213 
C  3.806331426277332  -1.8161877656493655  -0.7471742927672437 
H  4.682613234672566  -1.4857068161891256  -0.9034073696731113 
C  3.4983836288669794  -3.1323065564793806  -1.0367739172391248 
H  4.157576732591165  -3.711332609453276  -1.4017328354476988 
C  2.2232552059932047  -3.5963997511837764  -0.7873925390077711 
H  1.9913943036451718  -4.499256482639783  -0.9712964933614673 
C  1.2892185073005797  -2.7233284166116523  -0.26761197930602953 
H  0.41444893425284435  -3.048173596576361  -0.09044134651494026 
C  1.382760062128614  1.8599332878866912  -1.759728514835127 
C  0.36818751935415905  2.7528190000224133  -2.112708854956374 
H  -0.2091412961763499  3.1090108553993026  -1.4477577056819793 
C  0.20822274413257302  3.1175289454254327  -3.4475421540662183 
H  -0.4713230753139841  3.7362569845943026  -3.6894187813350943 
C  1.0336312733719801  2.5836014471786384  -4.425700242758824 
H  0.9154930133584134  2.8312303716790637  -5.336138729023411 
C  2.0230219012616684  1.6962564706799967  -4.076574952745754 
H  2.5866793444645753  1.3308784124825377  -4.7484979924682 
C  2.2059859973729097  1.3259534937292763  -2.744422928840214 
H  2.890833194395172  0.7107646856171939  -2.511484939013735 
C  1.9398438558451554  3.0304742786283048  0.7708122644790874 
C  0.9792785819164643  3.673414967565377  1.5411932643950093 
H  0.1471860436222081  3.2488036670210954  1.7108965088105264 
C  1.2309562261102318  4.939398290220701  2.0648061463714535 
H  0.5711464061183001  5.38025308264246  2.5879704450779064 
C  2.4425407720028818  5.549996806235871  1.82090167655957 
H  2.61595127977313  6.41395639353799  2.17678110460175 
C  3.4024281382747024  4.917983289747562  1.064503785864549 
H  4.236889902418009  5.343859341071116  0.9078602954414681 
C  3.159013243581822  3.6573992014983197  0.526300659252962 
H  3.8198020106640866  3.2284903512050698  -0.004210484303944967 
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

-Pd 0
lanl2dz


