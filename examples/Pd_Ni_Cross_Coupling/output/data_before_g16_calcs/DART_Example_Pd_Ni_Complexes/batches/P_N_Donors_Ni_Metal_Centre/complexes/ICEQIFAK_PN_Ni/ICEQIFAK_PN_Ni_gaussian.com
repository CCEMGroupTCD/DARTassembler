%chk=ICEQIFAK_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.3916249333063841  -1.4640628555495834  1.4557463356492593e-15 
N  1.3916249333063846  1.4640628555495834  -1.792959889799331e-16 
C  2.5940636824809813  1.0120762707580475  0.1402538115416189 
C  3.8987364620383946  1.6976292893292562  0.32087052452773546 
C  4.29336854947092  3.017547906672085  0.29430778444961336 
C  5.6417669459183495  3.2999088882604197  0.4899894149988356 
C  6.545309849990688  2.300148284095274  0.7100035392271281 
C  6.163525025178148  0.9664804585653112  0.7455237966733742 
C  4.830406433418579  0.6718012421883277  0.5325003478690828 
C  4.148152236588758  -0.6313964093112622  0.4582915071959217 
C  2.849136668404457  -0.42260947175505037  0.19615517182671002 
C  1.178644965052933  2.896296626833173  0.05400043997550749 
C  0.9061285152161639  3.4574840662150566  1.2852515648271374 
C  0.699891046853835  4.8436006610782005  1.3205663222027535 
C  0.7621023928119659  5.581693230604088  0.1575435033188022 
C  1.0241385792746494  4.978415781379786  -1.038941472032562 
C  1.2326941915895986  3.614168964770104  -1.1317442195257597 
C  0.8260412725810602  2.6282664199929786  2.5277862238049438 
C  1.5090495055364928  2.944299432678293  -2.4373509483796183 
C  1.690255344630712  -2.551527030304804  -1.408617381000342 
C  2.944813755716284  -2.7151792185232257  -1.9424337487391334 
C  3.1571908757321516  -3.57451957760072  -3.000333359361233 
C  2.120510149254251  -4.2528728882085405  -3.5516067918216683 
C  0.8605974249976901  -4.114402091526864  -3.029721846912541 
C  0.6475487124297485  -3.2615672386603225  -1.9721177906564713 
C  1.2900323148602348  -2.5002103127778104  1.4703303835386063 
C  1.6900052663517888  -3.812421781512075  1.4799428440666362 
C  1.6081813192890646  -4.563788801602932  2.644551487692021 
C  1.155890274249191  -3.997491311807815  3.7855377118106515 
C  0.7640945455362758  -2.690577955767686  3.7954652878616186 
C  0.8163105809378745  -1.944647061571174  2.64496302939788 
H  3.6778124556483203  3.699148150963102  0.14837275664501914 
H  5.932306382141706  4.182073759190893  0.4688319932453395 
H  7.440162283219655  2.519033227310163  0.8411347878981081 
H  6.784993321176773  0.29258333983154083  0.9053874899588583 
H  4.5523428568117374  -1.4604432773464242  0.5737631189663098 
H  0.5198758161826654  5.266691164093166  2.1308368529423234 
H  0.6224366379643091  6.500065505334928  0.18742613922223214 
H  1.0641102390890773  5.498471421145461  -1.809763661619505 
H  0.17967495776019615  1.927937625287569  2.4035387837050495 
H  0.5611706413651067  3.1816124337444798  3.266463228851551 
H  1.6850724083199842  2.2398721493203957  2.7093853932409204 
H  1.2690126973493459  3.533823035936438  -3.1565503726820046 
H  0.9908656657568249  2.1383556027891952  -2.4999079724048574 
H  2.442587875066579  2.7297278459399115  -2.495481119373035 
H  3.659671847054897  -2.2408651633802  -1.5852284067556737 
H  4.016032545150656  -3.688252723939915  -3.337631637290752 
H  2.264941240791645  -4.810525209558858  -4.2810435048337885 
H  0.15174343417545488  -4.594976192108006  -3.39067453570433 
H  -0.2121847554282923  -3.160954574086773  -1.6313876540788281 
H  2.0176728028707887  -4.202130766913404  0.7016601369421277 
H  1.865929570594525  -5.4561072222554365  2.64039395077691 
H  1.1123920526877407  -4.50031041360359  4.56598928980269 
H  0.45856161333522216  -2.3042315582370474  4.583661233331871 
H  0.5316977602571683  -1.0596145200457616  2.6548399087251977 
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

-Ni 0
lanl2dz


