%chk=UBEJOQAC_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5197065555560383  -1.4773598021470606  1.3151324929422968e-15 
N  1.5197065555560385  1.4773598021470606  6.887578525182374e-17 
N  2.0103670643939147  4.313326595010254  0.47018259734301215 
C  1.474206798510027  2.1078002091802053  1.3537244643163808 
H  1.7495511227988234  1.4369877005687453  2.028251773898692 
H  0.5406050013384809  2.3709827591825348  1.5534090147762112 
C  2.3673876231396376  3.3287535256157414  1.484071783425981 
H  2.261889951827708  3.72412700142787  2.3856027956335475 
H  3.315505237902972  3.0626155378157587  1.373150099228259 
C  2.2103617829694597  3.74410967928136  -0.8540174537784246 
H  3.1537184465112555  3.460625774824678  -0.9587550843860678 
H  2.011261231971934  4.424221557606302  -1.545635660272307 
C  1.2910449635023522  2.5438067859557942  -1.0286579980949673 
H  0.35027493409703125  2.8492187575988988  -0.9720979472658069 
H  1.4303290725280378  2.1587559512347783  -1.9304119621303855 
C  2.856504825134762  0.8225930089838941  -0.25005909905899415 
H  3.5801219793502934  1.4307310662765662  0.044878194949958546 
H  2.9690405410298624  0.6617823251851811  -1.2200444941161555 
C  2.9716533518138393  -0.4946832477355118  0.49785362331371347 
H  2.9704423456446953  -0.34279859196384793  1.4766085695618543 
H  3.8084162674203528  -0.9623267723501261  0.24955723398895768 
C  2.744013748129152  5.560492391942229  0.6406808038066965 
H  2.5816888370206783  5.915229358684087  1.5404800690100164 
H  2.4411027202949453  6.211745487905528  -0.027098451797015293 
H  3.7026639361166795  5.394765830515434  0.5230278229301354 
C  1.2809057094559295  -2.768669334345975  1.2431171859298438 
C  1.1177658246192883  -4.110211853046991  0.9058358664162571 
H  1.2051957632248886  -4.389000902540163  0.001097273386130389 
C  0.8230149493017984  -5.038262665724727  1.8953879977301038 
H  0.7201650218710354  -5.955297296036904  1.6684689169200195 
C  0.6797221028296634  -4.632971230199127  3.211964752970578 
H  0.4717421763017642  -5.272378061647467  3.8828806896263184 
C  0.8381400475091312  -3.302904839368493  3.556640108611663 
H  0.7332170152688979  -3.0281209482608404  4.46051245606638 
C  1.148439567370925  -2.3662687518094647  2.571914229264674 
H  1.2700004057251473  -1.4538221642950875  2.8063878898578967 
C  1.9178155616316652  -2.2483412268106098  -1.5837026368940925 
C  3.0224439974104507  -3.0886124614792996  -1.722380721165819 
H  3.5935259218849493  -3.2479386693928465  -0.9795000207768492 
C  3.2895835706452186  -3.6979949109622345  -2.944625741067375 
H  4.033815956094713  -4.2823970239378655  -3.0347892968856227 
C  2.468650731524549  -3.44796941652061  -4.0296260110132245 
H  2.646587291431229  -3.8699834574427388  -4.8621241490971245 
C  1.3941600350613923  -2.5882930714933354  -3.9113611175972407 
H  0.8492676810468942  -2.4023709964505366  -4.668137725684693 
C  1.1087830715174238  -2.0003142638003317  -2.6933258141752123 
H  0.3583066506910777  -1.424107967431058  -2.609334827705936 
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


