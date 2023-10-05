%chk=ZOZONACU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.4600159065914313  -1.5363767612470576  -1.132704751393463e-15 
C  1.0377003278346741  -3.1527099209495897  0.6930406328530707 
C  0.7757768759226294  -4.241279891015065  -0.1282218681842913 
C  0.3450132014030425  -5.442189109078473  0.424025908347877 
C  0.17436957013785026  -5.550307247278619  1.800796546217826 
C  0.43457242032476295  -4.4605801749678875  2.6215292887276616 
C  0.8665944672868788  -3.262487918934342  2.067913585399854 
C  2.115500319709614  -1.8607060908546875  -1.6512236463045993 
C  3.2797916318608387  -2.6082769757903175  -1.7859597307281334 
C  3.8396190149696157  -2.8010487668447475  -3.0428090286534153 
C  3.236029787453242  -2.2461796998926773  -4.164722239630967 
C  2.0722566654542307  -1.4984263857462765  -4.027816316838066 
C  1.5119847405194515  -1.3072254723963757  -2.7729731908401924 
N  1.4600159065914318  1.5363767612470576  -1.8815188829455864e-16 
C  2.5682364867413714  1.5027184481923697  0.8367541536762645 
C  3.504625984606669  2.5285531325966857  0.7889441279419943 
C  3.3303863068503956  3.5867872493899653  -0.09409562130497884 
C  2.221683224941205  3.619253195411899  -0.9287724081596372 
C  1.2875998728943518  2.5947687917800653  -0.8825207255517601 
C  2.6779578393750616  0.32797172881670034  1.7389787709930054 
C  2.890267082490193  -1.0011521173045501  1.0270660467355044 
H  4.2493352175624794  2.5174171043570293  1.3387917861396599 
H  3.956243825469701  4.283849006644562  -0.13441550425162715 
H  2.1014981921574583  4.338955639428619  -1.5225337069751956 
H  0.54200419976224  2.628388240669722  -1.4431075652680032 
H  3.6722046128792103  -2.9892975625146607  -1.0419343851783094 
H  4.603235991619024  -3.3119209206014304  -3.1397060125531255 
H  3.594160363420113  -2.3854142713830497  -5.0117850027480415 
H  1.6548970004532284  -1.131863427879826  -4.7823287938392465 
H  0.7242813622124338  -0.8089662404755286  -2.682421241057904 
H  0.8768084992304739  -4.162906243979252  -1.0510930198312078 
H  0.1623629851255357  -6.168346131978533  -0.1273378341912074 
H  -0.11624537662290457  -6.348831252552035  2.1699244586294464 
H  0.3226076607310806  -4.530894171329246  3.538996871002283 
H  1.0414630006151924  -2.5295544559952186  2.6126898192387564 
H  3.4092476946176733  0.4946506232515138  2.3397542563195026 
H  1.842938222277277  0.2840504969274414  2.23740366448373 
H  3.6739718234654966  -0.9079160611410025  0.4748142863906978 
H  3.053183292521986  -1.6586804108072617  1.723588337658287 
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

