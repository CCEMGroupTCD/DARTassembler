%chk=QONINUSO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
C  3.001141627303041  1.1293550938471006  -0.36198546381080177 
H  3.0726041447249073  1.1347687665488384  -1.337348642709637 
C  3.446445570296243  -0.27781388526824125  0.07235583352369515 
H  4.388460547743373  -0.3756890348507913  -0.14161483586699836 
H  3.359960956158361  -0.34258983123974407  1.0352965681179536 
C  2.690714592937833  -1.436676246752775  -0.555455552957488 
H  3.1726805528980613  -2.254477049829581  -0.31206822804807705 
C  4.0669650419329315  2.1221894327337685  0.13625418169544343 
H  4.233880361613919  1.9690180646020308  1.0708256920896322 
H  3.751360666772796  3.018191715102331  0.010828552854644865 
H  4.877895427918083  1.9941742780401488  -0.35718124392248995 
C  1.068773568970681  2.060258702017727  1.714659426253912 
C  1.8355635383182114  3.0465618423612058  2.37245270552684 
H  2.511597535696992  3.479246645760585  1.9025506888304689 
C  1.6214003584349235  3.375389272379689  3.6447470989947264 
H  2.1471959295990795  4.040403211439339  4.028597392990329 
C  0.6315552328496592  2.754893891959326  4.426840808362153 
H  0.515505490827227  2.970378544287155  5.3212655575655905 
C  -0.15169636023739508  1.817533416194939  3.8032944344641266 
H  -0.8337748349306038  1.4056070427226444  4.28020015596279 
C  0.05582654676064691  1.4724351548658314  2.473784598810849 
H  -0.4939043100160936  0.8343707560783148  2.0802842073149206 
C  1.0567830397170475  3.0460499742220093  -1.0072065099843541 
C  1.3893586253428212  3.0282444135745394  -2.393396991792599 
H  1.6949066749805015  2.234895033904833  -2.7707735824024478 
C  1.2664085460701937  4.135338869340262  -3.1780562261228016 
H  1.4677551684359333  4.082025761954682  -4.086199656558109 
C  0.8672778676819086  5.290671600600139  -2.664639366593511 
H  0.8311409456460546  6.058639279207789  -3.1919248750988745 
C  0.49894690000415765  5.3299769436702435  -1.301858541128256 
H  0.16394264442955686  6.127179230666133  -0.9604839926401847 
C  0.6089968113983052  4.252527760548415  -0.4500211594689101 
H  0.39746871409306894  4.3226599930175365  0.45079324651011576 
C  2.654083749016907  -1.3872935629367185  -2.0335991325237117 
H  2.279153687719785  -2.2064417930679734  -2.3680714127033595 
H  3.5441485453262103  -1.2806738937695263  -2.373149488903403 
H  2.110546227731188  -0.6492068750863224  -2.314470623018221 
C  1.3113539798184992  -1.9585718916420898  1.4238811881810052 
H  2.108850926678611  -2.490031974951828  1.576998592799862 
H  1.3889979527341425  -1.1533800601777056  1.9608083042898077 
C  0.13529851107592905  -2.718605509082509  1.8761281171934205 
H  -0.6523584758846785  -2.1934166839676124  1.7440322850616445 
H  0.23562599520595517  -2.927279248876102  2.8079872676280364 
H  0.07279803161916965  -3.532166326390107  1.3707588001494604 
N  1.2693178975339434  -1.5712835756158086  2.3667999921317346e-15 
H  0.9317324085740386  -2.3783378652929663  -0.4412304626353508 
P  1.2693178975339459  1.5712835756158086  2.961786478196049e-17 
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


