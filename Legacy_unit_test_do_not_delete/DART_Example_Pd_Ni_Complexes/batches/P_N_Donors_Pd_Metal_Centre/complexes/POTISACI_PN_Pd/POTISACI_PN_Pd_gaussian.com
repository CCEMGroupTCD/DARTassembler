%chk=POTISACI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
C  2.8311722108828024  -0.8907339210728398  -0.44249366971216564 
C  2.9830689648671656  0.47777892376240355  -0.3434236457394115 
C  4.211917904723654  1.1246152818221753  -0.629853268693486 
H  4.2764137115156275  2.070799019202501  -0.5615205091761759 
C  5.294415686967804  0.40231931337051974  -1.0023428194127264 
H  6.127383372616128  0.8358224804079979  -1.1496926373092886 
C  5.187878485910574  -0.9944649587881745  -1.1706236615813908 
C  6.325223077541768  -1.713006385355819  -1.6202650166560508 
H  7.152225645437384  -1.2574266043625186  -1.7262103603204564 
C  6.255984092036633  -3.0455432371138786  -1.9034505240411301 
H  7.026657154376548  -3.517884940756375  -2.1953875480440694 
C  5.0447049724182635  -3.69815733640681  -1.759005580271257 
H  4.952281517826829  -4.595657981155922  -2.0560950074078077 
C  3.96282366209725  -3.047205932015153  -1.181940705669947 
H  3.20345274055162  -3.5592794570194926  -0.9319437486872428 
C  3.9466400697861808  -1.6748530081382247  -0.9538828102434626 
C  1.7057697980725341  -2.115203987914527  1.3607405565184802 
C  1.8519810949639772  -0.9324723679200623  2.4395159556292736 
H  1.9044648487911249  -1.3160868167765338  3.3400215304109 
H  1.0733456990748165  -0.3398656955172636  2.3819550132096095 
H  2.6671633500677037  -0.41912071053039995  2.25456883664844 
C  2.802028152227603  -3.0664977261594064  1.6008190739275376 
H  3.6609946182592177  -2.6238236799148282  1.4386163642872247 
H  2.710106687076447  -3.830059927199949  0.9937418223469827 
H  2.7660013886108095  -3.379917480032153  2.5289422716689196 
C  0.3486804768078131  -2.802778763964394  1.799040900403029 
H  0.4312930700821551  -3.132549970064962  2.7183605189576894 
H  0.15421740260560646  -3.552754839715522  1.1990774110365618 
H  -0.37951738952865965  -2.148956591905661  1.7504926276486963 
C  0.8698705030633344  -2.441812230114032  -0.9899267433626872 
H  1.1062792250719782  -2.1830642662353315  -1.9050404780136754 
H  -0.10352589235272891  -2.4109401045874876  -0.8827324372801885 
H  1.1883624738064174  -3.352440196366424  -0.8153827040739557 
C  1.795258173412915  2.433034765488565  1.5125168482152445 
C  0.720885676284164  3.1084966912003575  2.073139336461212 
H  -0.12869892005744155  3.0825730400990743  1.6485778356481051 
C  0.8827721618257836  3.817300302447909  3.246020266684518 
H  0.14671734582782214  4.293646561128542  3.613351193353069 
C  2.0990317632391204  3.8415960152248214  3.8883234274353873 
H  2.1963955136389286  4.310588919288311  4.70945282408744 
C  3.1683607903016267  3.1848163843143706  3.339518099197089 
H  4.0150986138997045  3.2214609866745856  3.7688058129815234 
C  3.024451064381303  2.4625802861029533  2.156168911836945 
H  3.7648421077493976  1.9923535235593313  1.7917975289210049 
C  1.4710935258095013  2.6575212164544277  -1.3803048390957997 
C  0.7845764785656998  2.302458073614616  -2.5356682141603075 
H  0.2980722112133376  1.4874939399260596  -2.5655772203157468 
C  0.8078049802661061  3.1309238204823937  -3.6400671606355366 
H  0.32975008421786334  2.891547546750267  -4.424503305862428 
C  1.5265188975666066  4.304022621245631  -3.597557412030676 
H  1.5527347130146283  4.8666115078627525  -4.361786105338717 
C  2.210766526192892  4.671391839268251  -2.453309353099048 
H  2.7011349283722965  5.483178481892739  -2.431773638564971 
C  2.176667062646413  3.8490184689037465  -1.342072633807778 
H  2.6377438917156053  4.102012671550833  -0.5502076761777043 
N  1.523833669565022  -1.4731024904941266  1.8040302497996264e-16 
P  1.523833669565022  1.473102490494127  -1.804030249799627e-16 
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


