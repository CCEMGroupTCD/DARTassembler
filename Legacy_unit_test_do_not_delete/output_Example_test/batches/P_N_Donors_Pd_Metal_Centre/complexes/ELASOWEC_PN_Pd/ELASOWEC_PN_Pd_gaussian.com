%chk=ELASOWEC_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  3.154204805641519  1.1293550938471013  -0.3619854638108013 
H  3.2256673230633686  1.134768766548839  -1.33734864270963 
C  3.5995087486347073  -0.2778138852682408  0.07235583352369138 
H  4.541523726081832  -0.3756890348507911  -0.1416148358669999 
H  3.5130241344968294  -0.3425898312397443  1.0352965681179496 
C  2.8437777712763133  -1.4366762467527754  -0.5554555529574914 
H  3.325743731236539  -2.254477049829581  -0.3120682280480778 
C  4.2200282202713915  2.122189432733769  0.13625418169544543 
H  4.386943539952367  1.9690180646020314  1.0708256920896133 
H  3.9044238451112596  3.0181917151023323  0.010828552854646642 
H  5.030958606256558  1.9941742780401501  -0.3571812439224955 
C  1.2218367473091623  2.0602587020177277  1.7146594262539137 
C  1.9886267166566842  3.0465618423612058  2.3724527055268276 
H  2.6646607140354677  3.4792466457605853  1.9025506888304615 
C  1.7744635367734063  3.3753892723796897  3.6447470989947015 
H  2.3002591079375563  4.040403211439339  4.02859739299033 
C  0.7846184111881385  2.754893891959326  4.42684080836216 
H  0.6685686691657053  2.970378544287155  5.321265557565603 
C  0.0013668181010935143  1.8175334161949388  3.8032944344641195 
H  -0.6807116565921218  1.405607042722644  4.2802001559627865 
C  0.2088897250991375  1.472435154865831  2.473784598810842 
H  -0.34084113167759744  0.8343707560783143  2.0802842073149064 
C  1.2098462180555292  3.0460499742220106  -1.0072065099843577 
C  1.5424218036813013  3.0282444135745408  -2.393396991792581 
H  1.8479698533189757  2.234895033904835  -2.7707735824024264 
C  1.4194717244086579  4.135338869340263  -3.1780562261227785 
H  1.6208183467744162  4.082025761954684  -4.086199656558103 
C  1.0203410460203926  5.290671600600142  -2.664639366593494 
H  0.9842041239845383  6.0586392792077906  -3.1919248750988585 
C  0.6520100783426354  5.329976943670245  -1.3018585411282597 
H  0.31700582276803657  6.127179230666135  -0.9604839926401819 
C  0.7620599897367845  4.252527760548416  -0.45002115946890714 
H  0.5505318924315463  4.322659993017537  0.45079324651011565 
C  2.8071469273553733  -1.3872935629367185  -2.0335991325236913 
H  2.43221686605826  -2.2064417930679734  -2.368071412703341 
H  3.697211723664677  -1.2806738937695255  -2.3731494889033904 
H  2.2636094060696563  -0.6492068750863215  -2.3144706230182077 
C  1.4644171581569774  -1.9585718916420907  1.4238811881809994 
H  2.2619141050170946  -2.490031974951829  1.5769985927998627 
H  1.5420611310726255  -1.1533800601777056  1.9608083042898028 
C  0.2883616894144154  -2.7186055090825105  1.8761281171934112 
H  -0.4992952975461846  -2.193416683967613  1.744032285061634 
H  0.3886891735444471  -2.9272792488761015  2.807987267628009 
H  0.22586120995764958  -3.5321663263901075  1.3707588001494595 
N  1.422381075872425  -1.571283575615809  2.8274792242380584e-16 
H  1.0847955869125179  -2.3783378652929663  -0.4412304626353533 
P  1.4223810758724251  1.571283575615809  -4.144713450681022e-16 
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