%chk=OGICUGIG_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.3831348858661612  1.4720862364345368  -1.6254086234253677e-16 
P  0.6169366767041696  1.0935252383032426  -2.8740043119263046 
N  1.3831348858661618  -1.4720862364345373  -3.3306690738754696e-16 
N  2.435679813315051  0.898480761360489  2.6822425215977113 
N  -2.2575929392986245  1.718490812486106  -2.1556783567991413 
N  2.4171707448121658  -1.024554509899509  -4.241597815229405 
C  2.8083273497679313  0.48878356423108715  -0.4950709804414092 
C  4.0218288560201945  1.054736940045459  -0.9070317561565945 
C  5.076643991804291  0.23563687297298164  -1.2682108799202019 
C  4.9336567301056995  -1.1362240887131487  -1.2265570110297057 
C  3.7328399489246014  -1.7042331491233695  -0.8309584451503466 
C  2.6649544347638776  -0.8874260217689997  -0.4547135443891504 
C  0.7849222668524556  -2.2623087153957617  -1.1160212839355617 
C  1.6365243872677615  -2.3771016370450333  1.1553868682122475 
C  2.008463000126765  2.732253368010387  1.1331683943753483 
C  2.0737144823141276  4.077898892010291  0.788591120296831 
C  2.604394613209495  4.998589862046591  1.6783420440763703 
C  3.041374921498501  4.591951038732914  2.9192212676681466 
C  2.9623453801743485  3.258347587786377  3.2843064320087185 
C  2.472806743979591  2.306416097417934  2.3851966993745424 
C  1.5011041392547317  0.5922872016056457  3.7739558097001322 
C  3.7759627227986985  0.35933350533946284  2.957994817818561 
C  0.7772960696630352  2.34031635323506  -1.4857868430775056 
C  -0.631873838080113  1.9076947722178392  -3.964727630742573 
C  -0.33713527826603107  2.2121480976463257  -5.300045575650922 
C  -1.312239426369601  2.722522871490909  -6.156569708746055 
C  -2.5840797766443075  2.9457663791511437  -5.681437274891407 
C  -2.904200137564084  2.637152142596394  -4.362878536146584 
C  -1.9453340303699118  2.0863130919141484  -3.5065501769771554 
C  -2.535678065149018  0.28738981251851103  -1.9891868534948514 
C  -3.280193888999264  2.522265710800372  -1.5067867629937859 
C  2.1951052756717786  1.3746857537377641  -3.787370091772076 
C  2.7271650058110164  2.664458205926867  -3.8689296107645026 
C  4.02128416648457  2.8898086621633694  -4.327484645467692 
C  4.793937903784502  1.8122336697892454  -4.702612181175715 
C  4.273119289552869  0.5318312053485563  -4.697995892182726 
C  2.959357778102067  0.280989518319156  -4.258557217190667 
C  1.2252905224694244  -1.2228866984451643  -5.058875227479959 
C  3.351274035416512  -2.1264801782821405  -4.443474272107459 
H  4.124264402284329  2.0298847774961053  -0.936987957132256 
H  5.925614990717301  0.6289259790529735  -1.5598997587679726 
H  5.683854834958026  -1.7143924248708444  -1.4748224708517084 
H  3.6297073154468444  -2.676780340697948  -0.8196151620126678 
H  0.6320008782551799  -1.6763430065598182  -1.8840750235774482 
H  1.398486726841691  -2.9826554340556477  -1.368025010800972 
H  -0.06711832707069676  -2.6453035128782476  -0.8202845912367158 
H  0.7895703971629637  -2.763356939843854  1.460254035525114 
H  2.2433542975438967  -3.096709108931825  0.8811779763087078 
H  2.044713386868482  -1.8684605252216242  1.887126630674297 
H  1.7464019084012155  4.374400736941039  -0.08602877709356094 
H  2.6681273249612323  5.940835496622691  1.4217062991477625 
H  3.4085670110245143  5.244619551772658  3.548163213698157 
H  3.256179381246805  2.976274887340701  4.176045761340479 
H  0.5942899068636991  0.8763071792753401  3.517760297995503 
H  1.501877308308768  -0.369227484129236  3.948502426183428 
H  1.7721730158059326  1.0737767043651671  4.583195133323055 
H  4.163536030779534  0.8257038734660109  3.7283174454716317 
H  3.707922374467469  -0.5971456839574166  3.156901105906568 
H  4.349180719407738  0.49093418552054824  2.1743172506141577 
H  -0.09542102528481156  2.745559962770383  -1.301154047733355 
H  1.4131081440726136  3.0421982161499215  -1.7424322727000878 
H  0.5711168561304221  2.0603850219034903  -5.642083993384066 
H  -1.0909448592475814  2.923758164416068  -7.088873822993981 
H  -3.2684004523442853  3.3203694739641443  -6.2781265693874015 
H  -3.8079218803614463  2.810333902394293  -4.027256934236611 
H  -3.4089039899839424  0.07553038317191407  -2.38455671374364 
H  -1.8404656953948293  -0.23707233575206277  -2.4384057266204198 
H  -2.548671925069656  0.06439890554891203  -1.0354816028994387 
H  -3.0752100260555073  3.4719368580464844  -1.6318420228627923 
H  -4.153906803207475  2.3230041779227473  -1.9057169386650363 
H  -3.3018044810758607  2.3133001030890603  -0.5508582473065714 
H  2.1727214452035377  3.4256450761431005  -3.5959056251673993 
H  4.381904521432382  3.8018679843481302  -4.3766177546017655 
H  5.7234992164879746  1.9541373095555885  -4.980446867538059 
H  4.8270130449356285  -0.2154022941400564  -5.006853537368987 
H  1.4770399937478558  -1.2510765992627353  -6.002842929895116 
H  0.7946196310220623  -2.067104954306336  -4.807448638331821 
H  0.597708986868454  -0.48274521553928074  -4.909494835678936 
H  4.147797705519709  -1.9902911629663147  -3.8891816552770884 
H  2.9204722955127798  -2.972152238239179  -4.192446847502215 
H  3.613995261037184  -2.16466746215768  -5.388172000991468 
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

