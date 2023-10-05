%chk=ILOLEKAJ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.4646129147320792  -1.5319951076945384  -1.0125248334586006e-15 
N  1.4646129147320806  1.5319951076945384  1.1769604127716595e-16 
C  1.0842893312542439  -3.2033173117537728  0.6148239033001048 
C  2.160008503519347  -1.8273345168920614  -1.651438508837109 
C  2.921635830993499  -1.0653188197552745  0.9969543890090321 
C  3.3239823128484347  0.2802058948125513  1.0304678344620168 
C  1.3380207959350712  -4.379232705315982  -0.09317187423982924 
H  1.7165362607768047  -4.3354758619270495  -0.9420484798382549 
C  1.1635313996462757  2.794884297276575  -0.6321900870074407 
C  4.471997307195428  0.6230782364811086  1.7645384714399752 
H  4.770943762354264  1.5035934504167414  1.7551143016698911 
C  3.6545522548358504  -2.0037233047170164  1.7134292791942065 
H  3.402314561681554  -2.898491364428939  1.6852449666859495 
C  0.2407143307369224  -4.539728870587695  2.4505287689899475 
H  -0.1172506448681272  -4.591132853279593  3.3068005522477666 
C  5.164979900138929  -0.31444900024840916  2.4963880026057583 
H  5.900894950872983  -0.06196377489929369  3.0047253995022687 
C  1.2611777748814705  -2.11491033080115  -2.680963619026702 
H  0.35340175479539404  -2.1905592473139492  -2.4901293223524883 
C  0.9856312257350077  2.82842601013317  -2.018973030232893 
C  1.701355304964434  -2.288598274492061  -3.978442859795697 
H  1.0958012023517958  -2.496962079569986  -4.6523433669496885 
C  0.529678297480847  -3.3051173104722382  1.8830597376297222 
H  0.34668407142037827  -2.5301558278689145  2.3640528183593643 
C  1.0460064284758621  3.946707864556424  0.15473771845512602 
C  4.758385370144749  -1.6281227906488407  2.4720004936356554 
H  5.222259695050067  -2.2671917831558392  2.9631795073462395 
C  3.951802305931341  -1.8720570291611165  -3.264090575010843 
H  4.856911483816542  -1.7915324639190997  -3.4607503332242544 
C  0.49218004588685593  -5.689154838827212  1.7312730302983526 
H  0.29700141375821687  -6.52031163796099  2.102086914631309 
C  0.7140081522065139  4.048394801556134  -2.6102626326113976 
H  0.5986801003989011  4.096973507601742  -3.5317785574793215 
C  2.6910624844306352  1.4116619688521932  0.35620875716724015 
H  3.245462671884659  2.132836047314258  0.1620015979136758 
C  3.508360857177751  -1.7108562206470201  -1.9540404038869574 
H  4.1195399741160115  -1.5238474968194737  -1.2782294730557215 
C  1.0326138563393532  -5.608767460707366  0.4592311457839894 
H  1.1899458725881025  -6.385582539308313  -0.026223455925904467 
C  3.050423990569975  -2.151504240088737  -4.273925189638952 
H  3.346813956482582  -2.2477021974188887  -5.149742956706866 
C  1.1319101587013447  1.5770300165225604  -2.8540622475896593 
H  0.9299113760325819  1.7777334120201571  -3.770417429184725 
H  0.5250844695243955  0.9056652681404294  -2.533261970805615 
H  2.033190199428069  1.2515851064860648  -2.7879175673879177 
C  0.6105927380633417  5.204132460969752  -1.8446969497789252 
H  0.4376143401898087  6.0196183876056875  -2.258576772104015 
C  0.7632956250255303  5.151424590301995  -0.48797003964441904 
H  0.6763237386309466  5.928754311538574  0.01425230699024851 
C  1.2448911684980226  3.9033516360238436  1.6504502869768254 
H  0.7824738269787858  3.1435289615560755  2.0126731957318302 
H  0.896535515832515  4.7070170223009695  2.043528098810616 
H  2.182032439863539  3.8312471368350476  1.8473644527957875 
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