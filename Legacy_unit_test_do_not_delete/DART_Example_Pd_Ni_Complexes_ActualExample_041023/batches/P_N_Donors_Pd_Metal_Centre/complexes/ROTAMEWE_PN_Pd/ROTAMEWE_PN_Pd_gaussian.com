%chk=ROTAMEWE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5306830844430195  1.465984070513728  2.4214685298137185e-15 
N  1.5306830844430208  -1.4659840705137284  -5.551115123125783e-17 
C  2.7134258826303617  -0.9842149816091661  0.12263428685356109 
C  3.928170024311562  -1.533119305459024  0.8112174143771527 
C  5.0965015000655285  -1.5557216895274868  -0.21183145786554097 
C  5.458234964576821  -0.054596381099797364  -0.3858244147647441 
C  4.395720750683436  0.6588551331676662  0.45202649205819817 
C  3.0401800764978315  0.43717522781499984  -0.26900924921355 
C  4.272311391559613  -0.26953119269674897  1.6696024157806784 
C  1.22144796802109  -2.8417070306501575  0.32038433127258004 
C  1.4827430081063866  -3.8382433535094176  -0.6393723041558954 
C  1.0266518506623048  -5.128478607796437  -0.35852397374207384 
C  0.3717041551396605  -5.420730736327316  0.8100886992590317 
C  0.16709463716695838  -4.424726816933893  1.7609707569864517 
C  0.5748266500621552  -3.133540941812697  1.5365251699325786 
C  2.2468191721658775  -3.578402307033031  -1.9300978311686063 
C  1.4049021892177975  -3.8561333841908625  -3.177866360425691 
C  3.5277130733259945  -4.4174243944756135  -1.9848003974188633 
C  0.3737737567822932  -2.0598840640483758  2.5915052990976815 
C  -0.9837316449529077  -2.1205673564451515  3.2832650350491535 
C  1.5025751884273233  -2.1246469914973285  3.6277741595506052 
C  1.5752360747144438  2.087647377319398  1.6934060777112476 
C  2.534478465182917  3.045934604483947  2.0905033159056727 
C  2.549936267922045  3.520059800644112  3.4058129694475996 
C  1.6182241337882415  3.051073667828175  4.3022444762199745 
C  0.6596715432064257  2.1081801959259816  3.9357938769120695 
C  0.6567499932164471  1.647233369427929  2.6291705885511623 
C  1.6905063196385877  2.9273452736059675  -1.0658569404889469 
C  0.49972345828298215  3.5794281674296413  -1.4257421675882758 
C  0.5172794739624846  4.67719744908087  -2.2308863103207464 
C  1.7240864567081042  5.166583070268819  -2.7102631142765983 
C  2.9020106129702157  4.538617576975707  -2.3524210942071657 
C  2.8910281150781265  3.4148829620267818  -1.5326985314423842 
H  3.8054104978205174  -2.3332113962202623  1.306031263169823 
H  4.792108143297549  -1.9129667445958365  -1.0290847509055598 
H  5.838865812359087  -2.0131949141614927  0.17727313806701378 
H  5.34734739836214  0.18431740182914225  -1.3004209607814066 
H  6.3162885911406565  0.11137075312338518  -0.022778819421174268 
H  4.593259331026903  1.5627904600608318  0.6641097400303992 
H  3.19642531123403  0.46982463442364564  -1.2057388862328693 
H  5.121483897778969  -0.38591926500384277  2.07243648994252 
H  3.5079031455964005  -0.022153951501409308  2.180454325870123 
H  1.1767159130356377  -5.81658669525096  -0.9936617046034186 
H  0.05555079858325618  -6.302881409729116  0.9742432135279697 
H  -0.2648257655801183  -4.645154403085465  2.5748024809652077 
H  2.4993986016993617  -2.6681079888059873  -1.942809085241403 
H  1.952508038867344  -3.809405966493892  -3.9451788533938523 
H  0.7198392781298844  -3.188271267069153  -3.238217502630907 
H  1.0039322976063296  -4.70607336311803  -3.10741305502151 
H  4.041690186541253  -4.152242485054341  -2.735018940249267 
H  3.290137953026553  -5.33948857797408  -2.0747730806288818 
H  4.017784957586748  -4.293651143012919  -1.1901774236259497 
H  0.43542122781959325  -1.2126071734033568  2.1549364615776057 
H  -1.6739135857317553  -2.0622746664229106  2.6278807348765625 
H  -1.0640234830738406  -1.402795825932368  3.8917313042181827 
H  -1.0623574646761056  -2.9495320221895316  3.743606716216682 
H  2.3457545866665015  -2.086058134266311  3.183212113844251 
H  1.445509857139796  -2.9474959905163134  4.10572248602378 
H  1.429763223831635  -1.3999108286905368  4.230473680150871 
H  3.1676864784210332  3.3686877261015593  1.4614144474163713 
H  3.1953852710187594  4.15869609542453  3.6754473837766963 
H  1.6274130775349662  3.375773809237756  5.191887776358538 
H  0.025909909072152626  1.7913269883483411  4.57273287394959 
H  0.0003656813944192816  1.0100464662491935  2.366642388901258 
H  -0.3306657983296426  3.242989914163244  -1.1009760298192437 
H  -0.30317665547725414  5.111880119116085  -2.4762851001218467 
H  1.7401989002167073  5.925842442445999  -3.2791397660433077 
H  3.733460029614941  4.8762573814972825  -2.6718691180917866 
H  3.7061663422200146  2.994703972780234  -1.2913461176409853 
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


