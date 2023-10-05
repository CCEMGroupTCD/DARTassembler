%chk=PILOSEVO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.3896009157308442  -1.4659840705137281  2.2036345395623093e-16 
N  1.3896009157308442  1.4659840705137286  -1.7953126995556455e-16 
C  2.572343713918185  0.9842149816091663  -0.1226342868535606 
C  3.7870878555993848  1.5331193054590242  -0.8112174143771521 
C  4.955419331353349  1.5557216895274875  0.21183145786554028 
C  5.317152795864647  0.054596381099797635  0.3858244147647423 
C  4.254638581971255  -0.6588551331676664  -0.45202649205819684 
C  2.899097907785655  -0.43717522781499957  0.26900924921355 
C  4.131229222847438  0.26953119269674897  -1.669602415780679 
C  1.0803657993089142  2.8417070306501584  -0.32038433127258104 
C  1.3416608393942107  3.8382433535094176  0.6393723041558932 
C  0.8855696819501273  5.1284786077964375  0.35852397374207434 
C  0.23062198642748433  5.420730736327316  -0.810088699259031 
C  0.02601246845478311  4.424726816933892  -1.7609707569864503 
C  0.43374448134997845  3.133540941812698  -1.5365251699325806 
C  2.1057370034537026  3.578402307033032  1.9300978311686072 
C  1.2638200205056214  3.8561333841908634  3.177866360425686 
C  3.386630904613815  4.417424394475613  1.984800397418861 
C  0.23269158807011592  2.0598840640483753  -2.591505299097683 
C  -1.1248138136651016  2.1205673564451515  -3.283265035049183 
C  1.3614930197151485  2.1246469914973285  -3.627774159550602 
C  1.4341539060022674  -2.0876473773193984  -1.693406077711246 
C  2.3933962964707405  -3.0459346044839477  -2.0905033159056736 
C  2.4088540992098686  -3.520059800644114  -3.4058129694475974 
C  1.4771419650760658  -3.051073667828176  -4.302244476220009 
C  0.5185893744942434  -2.1081801959259825  -3.935793876912077 
C  0.5156678245042685  -1.6472333694279295  -2.6291705885511636 
C  1.5494241509264137  -2.927345273605968  1.0658569404889455 
C  0.358641289570808  -3.5794281674296418  1.4257421675882758 
C  0.3761973052503076  -4.677197449080869  2.230886310320749 
C  1.5830042879959294  -5.1665830702688185  2.7102631142766014 
C  2.760928444258041  -4.538617576975707  2.352421094207162 
C  2.7499459463659512  -3.4148829620267813  1.532698531442385 
H  3.664328329108342  2.3332113962202623  -1.306031263169823 
H  4.651025974585371  1.912966744595837  1.029084750905558 
H  5.6977836436469556  2.0131949141614927  -0.17727313806702297 
H  5.206265229649976  -0.18431740182914186  1.300420960781405 
H  6.175206422428543  -0.11137075312338496  0.022778819421175114 
H  4.452177162314724  -1.562790460060832  -0.6641097400303975 
H  3.055343142521856  -0.46982463442364547  1.2057388862328706 
H  4.980401729066812  0.3859192650038427  -2.0724364899425325 
H  3.3668209768842203  0.022153951501409263  -2.18045432587012 
H  1.0356337443234607  5.816586695250962  0.9936617046034194 
H  -0.0855313701289202  6.302881409729118  -0.9742432135279694 
H  -0.40590793429229444  4.645154403085465  -2.574802480965204 
H  2.358316432987185  2.6681079888059887  1.942809085241404 
H  1.8114258701551695  3.8094059664938924  3.9451788533938488 
H  0.5787571094177058  3.188271267069154  3.238217502630909 
H  0.8628501288941539  4.70607336311803  3.107413055021509 
H  3.9006080178290765  4.152242485054341  2.7350189402492684 
H  3.1490557843143736  5.33948857797408  2.0747730806288813 
H  3.8767027888745735  4.2936511430129185  1.190177423625948 
H  0.2943390591074171  1.2126071734033568  -2.1549364615776057 
H  -1.8149957544439501  2.0622746664229106  -2.627880734876584 
H  -1.2051056517860474  1.4027958259323678  -3.8917313042182307 
H  -1.2034396333883124  2.949532022189532  -3.7436067162167257 
H  2.2046724179543222  2.086058134266311  -3.18321211384425 
H  1.3044276884276225  2.947495990516314  -4.1057224860237955 
H  1.2886810551194623  1.3999108286905366  -4.230473680150915 
H  3.0266043097088557  -3.3686877261015598  -1.4614144474163762 
H  3.0543031023065907  -4.1586960954245304  -3.675447383776708 
H  1.4863309088227852  -3.375773809237757  -5.191887776358598 
H  -0.11517225964004929  -1.791326988348342  -4.572732873949642 
H  -0.14071648731775732  -1.010046466249194  -2.366642388901258 
H  -0.4717479670418192  -3.2429899141632443  1.100976029819246 
H  -0.4442588241894303  -5.111880119116084  2.4762851001218493 
H  1.5991167315045327  -5.925842442445999  3.2791397660433033 
H  3.5923778609027623  -4.876257381497282  2.671869118091784 
H  3.565084173507837  -2.9947039727802336  1.291346117640986 
Br  -1.5839191898578664  1.5839191898578664  0.0 
C  -1.2937473829822488  -1.2398794889499034  -0.4598005675376139 
C  -1.601562203225019  -1.4938654694550033  -1.7919791134336869 
H  -1.1309417375171225  -1.0428217585005106  -2.480879663633267 
C  -2.60218398542082  -2.411889969420291  -2.121542410022959 
H  -2.8231586547736818  -2.5626138747145775  -3.03270050271447 
C  -3.2683208079706407  -3.09759503004543  -1.143064365570095 
H  -3.9596020908062854  -3.707191723557817  -1.3683280760515548 
C  -2.9347449854640226  -2.896087320001541  0.14198826503765638 
H  -3.3668350933519986  -3.405606940434718  0.8171022443470677 
C  -1.9675509241539992  -1.9516175659140573  0.5112421903797603 
H  -1.7755298364628027  -1.8022912010877177  1.4301324877763764 

-Ni 0
lanl2dz
F 1 1.0
3.130 1.0
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

-Ni 0
lanl2dz


