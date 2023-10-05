%chk=AQAHOFIY_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  3.154204805641519  -1.1293550938471013  0.3619854638108012 
H  3.2256673230633686  -1.1347687665488393  1.3373486427096297 
C  3.5995087486347073  0.2778138852682408  -0.07235583352369135 
H  4.541523726081832  0.3756890348507911  0.14161483586699997 
H  3.5130241344968294  0.3425898312397444  -1.0352965681179496 
C  2.8437777712763133  1.4366762467527754  0.5554555529574916 
H  3.325743731236539  2.254477049829581  0.3120682280480781 
C  4.2200282202713915  -2.122189432733769  -0.13625418169544568 
H  4.386943539952367  -1.9690180646020312  -1.0708256920896135 
H  3.9044238451112596  -3.0181917151023323  -0.010828552854647011 
H  5.030958606256558  -1.9941742780401501  0.3571812439224953 
C  1.2218367473091623  -2.0602587020177277  -1.7146594262539139 
C  1.9886267166566842  -3.0465618423612053  -2.372452705526828 
H  2.6646607140354677  -3.479246645760585  -1.902550688830462 
C  1.7744635367734063  -3.375389272379689  -3.644747098994702 
H  2.3002591079375563  -4.0404032114393384  -4.0285973929903305 
C  0.7846184111881385  -2.7548938919593255  -4.42684080836216 
H  0.6685686691657053  -2.9703785442871546  -5.321265557565603 
C  0.0013668181010935143  -1.8175334161949384  -3.80329443446412 
H  -0.6807116565921218  -1.4056070427226435  -4.2802001559627865 
C  0.2088897250991375  -1.4724351548658308  -2.473784598810842 
H  -0.34084113167759744  -0.8343707560783141  -2.0802842073149064 
C  1.2098462180555292  -3.0460499742220106  1.0072065099843572 
C  1.5424218036813013  -3.028244413574541  2.3933969917925806 
H  1.8479698533189757  -2.2348950339048352  2.770773582402426 
C  1.4194717244086579  -4.135338869340263  3.178056226122778 
H  1.6208183467744162  -4.082025761954685  4.086199656558102 
C  1.0203410460203926  -5.290671600600142  2.6646393665934935 
H  0.9842041239845383  -6.0586392792077906  3.1919248750988576 
C  0.6520100783426354  -5.329976943670245  1.301858541128259 
H  0.31700582276803657  -6.127179230666135  0.9604839926401811 
C  0.7620599897367845  -4.252527760548416  0.45002115946890664 
H  0.5505318924315463  -4.322659993017537  -0.4507932465101162 
C  2.8071469273553733  1.3872935629367182  2.0335991325236913 
H  2.43221686605826  2.206441793067973  2.3680714127033413 
H  3.697211723664677  1.2806738937695252  2.3731494889033904 
H  2.2636094060696563  0.6492068750863211  2.3144706230182077 
C  1.4644171581569774  1.958571891642091  -1.4238811881809992 
H  2.2619141050170946  2.490031974951829  -1.5769985927998624 
H  1.5420611310726255  1.1533800601777058  -1.9608083042898026 
C  0.2883616894144154  2.718605509082511  -1.876128117193411 
H  -0.4992952975461846  2.193416683967613  -1.7440322850616339 
H  0.3886891735444471  2.927279248876102  -2.8079872676280084 
H  0.22586120995764958  3.5321663263901075  -1.370758800149459 
N  1.422381075872425  1.571283575615809  -9.032118228073497e-17 
H  1.0847955869125179  2.3783378652929663  0.4412304626353536 
P  1.4223810758724251  -1.571283575615809  2.220446049250313e-16 
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
