%chk=BEKAQOVO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5442887739992146  1.4516446474602516  1.3339133067951025e-15 
N  1.544288773999215  -1.451644647460252  0.0 
H  1.249635861351797  -2.1038759953018205  -0.6983141914361897 
C  2.9852701069054923  0.4723216550962237  -0.47916236793268113 
C  4.205320487246788  1.0034274366308953  -0.9054259345784585 
H  4.329145785539281  1.9453006301787226  -0.9246129817402977 
C  5.227061511676119  0.17310269626911667  -1.2974116161298075 
H  6.057196976534058  0.5401917844507951  -1.578293639981537 
C  5.045953552728789  -1.1979607650829134  -1.282526988108667 
H  5.746271955150661  -1.769501542263864  -1.5745376304638603 
C  3.854353537036463  -1.7417068715291273  -0.8461872645827746 
H  3.7414564168267352  -2.684548900667994  -0.8166805138695372 
C  2.8277033317436313  -0.9060279626980895  -0.452931613247416 
C  2.0989120464290982  2.46038439142028  1.4100677411509999 
C  3.2115804395276015  2.096799266835269  2.1496100551399824 
H  3.7503477458997017  1.3663320985408784  1.869257267337639 
C  3.539427327747335  2.800172183744642  3.301647366880818 
H  4.296563288894119  2.5396097377894553  3.8125119654744104 
C  2.778689077622839  3.8697671412409584  3.708166286115837 
H  3.0129515638646076  4.345904583837029  4.496261282482324 
C  1.6845411945753896  4.250065519336358  2.9783603916200234 
H  1.1653699642663407  4.996175813542953  3.2546380825368875 
C  1.3324492779944945  3.5455468709265676  1.8318517233923033 
H  0.565581829311838  3.804466707580525  1.3347206598506427 
C  1.31296373833994  2.571995125688795  -1.4188618233677128 
C  1.841628686997076  3.8592808352247894  -1.4528411418682747 
H  2.2885527358930036  4.2160539862201905  -0.694084953592155 
C  1.706854581914751  4.615871683351122  -2.6114240776699362 
H  2.0555596705002674  5.4992335255509435  -2.639150544085699 
C  1.0780265726725267  4.103191130777818  -3.7127986124244057 
H  1.0069051577183865  4.623521806278539  -4.504255226818344 
C  0.5469988098640395  2.8343463461932807  -3.6773194951970627 
H  0.10001584147904685  2.4841948846153414  -4.439062279822285 
C  0.6666549609147461  2.0645957485143156  -2.523112454459313 
H  0.3015396231497347  1.1877392690535915  -2.4970723976911864 
C  1.6918131644173586  -2.214865206939727  1.2516718027580653 
C  1.4147451719394735  -3.572326847739774  1.2440428250599038 
H  1.1587233785925795  -4.011051921478895  0.44141795564452097 
C  1.5180328336934776  -4.269950124050396  2.4221262572832365 
H  1.3106472859181557  -5.19639323624933  2.4350709254166834 
C  1.9129497316751178  -3.6576886855812845  3.5735060836882013 
H  1.994127211024888  -4.159686489434279  4.375647799205584 
C  2.1959448225520974  -2.2994433474748823  3.570742019334119 
H  2.465493097181846  -1.8680706343969424  4.372736732568503 
C  2.0841137636284945  -1.5696530252679417  2.3898405162979257 
H  2.277510520543649  -0.6395282923573892  2.375852105281174 
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


