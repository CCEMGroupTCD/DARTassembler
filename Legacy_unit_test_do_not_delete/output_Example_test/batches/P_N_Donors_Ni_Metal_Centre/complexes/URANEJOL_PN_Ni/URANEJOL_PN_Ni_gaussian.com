%chk=URANEJOL_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.349771708475176  -1.5027362825858703  -1.941899428857607e-15 
P  2.802843666905349  0.991419710340871  -0.49878863917650973 
N  1.3497717084751777  1.5027362825858703  -2.950544203056535e-16 
C  1.264118310849352  -2.6222497258262796  1.4229595271402495 
C  1.3226390320644974  -2.0772169216734335  2.7131857738062335 
H  1.4437630347941428  -1.1415069862079932  2.8227382390181517 
C  1.20666235600023  -2.8851934702336886  3.825770443656466 
H  1.2493387935823583  -2.509363954073905  4.697338660411992 
C  1.0278218923025437  -4.247308312817876  3.6613530049334035 
H  0.935129372147256  -4.805161656557913  4.424527299302222 
C  0.9828342050377358  -4.80152430510554  2.404990397288929 
H  0.8765039670170062  -5.740644366525469  2.304974994869705 
C  1.0932311188483004  -3.990000018288224  1.2741111980254645 
H  1.0517339694521548  -4.373597354480884  0.405913478460448 
C  1.516249751525352  -2.504863324070907  -1.4924528380845732 
C  2.654895115193776  -3.307317995575182  -1.6701430909894834 
H  3.309550522262888  -3.3615109840214745  -0.9838671202616344 
C  2.824128486763457  -4.0193220227614095  -2.8432738898426035 
H  3.5973740039732096  -4.556872166696156  -2.961705105621108 
C  1.8720861780672975  -3.9522992482618893  -3.845481193822638 
H  1.995393538443735  -4.44037879255684  -4.65081014730553 
C  0.7418360980846104  -3.1748771393142796  -3.674745078080783 
H  0.08542622654947163  -3.1369571990956633  -4.360922277268475 
C  0.5618420825125768  -2.447692121781066  -2.498254804158592 
H  -0.2147593244261261  -1.9111224060272807  -2.386225679031526 
C  2.966810209248646  -0.6495812081679877  0.22354882243341798 
H  3.1891811792676386  -0.5813084971215384  1.186448504964892 
H  3.688026384541211  -1.1544500888836873  -0.2312126664607623 
C  3.003216195511132  0.8026958715100253  -2.281748075741785 
C  1.8705807568424884  0.6947449676289711  -3.0805499882477507 
H  1.006245884204032  0.7851437765046309  -2.695425328825672 
C  2.0020790004300952  0.4573453792021075  -4.433125393668153 
H  1.2277187265144067  0.35729409825394237  -4.974646812784606 
C  3.255429597957619  0.366472911036891  -5.000727863894415 
H  3.3446077028294483  0.22638623885130044  -5.935484538054832 
C  4.382561315360823  0.47741033693914664  -4.209076457555799 
H  5.242560503353473  0.4097379374779394  -4.604921292059409 
C  4.2715508650788845  0.6866120569235946  -2.841113702111135 
H  5.048144192742715  0.7501027878692942  -2.2972007030827775 
C  4.144266643295528  2.045311400022292  0.06563749844737629 
C  4.457184924319135  3.2171494599710275  -0.6450806951123721 
H  3.982941398822708  3.431875392750628  -1.4398233278663124 
C  5.453513841262126  4.055792942932241  -0.18776753353715242 
H  5.644998479396158  4.859863653144133  -0.6565823438173186 
C  6.168008983696633  3.7357800567240975  0.9382775561898399 
H  6.867997604593759  4.307548007599339  1.2329231754793208 
C  5.876971389302111  2.5932447297367496  1.6416421524235194 
H  6.371460297445857  2.3818094042694913  2.425227886639837 
C  4.867977767779304  1.7457443567448399  1.2170471930741373 
H  4.668483329487162  0.9596002616763861  1.7113527775858992 
C  1.0312237303892526  2.888669276375145  -0.04856317929945427 
C  1.1723657096702131  3.662708323331831  -1.1863066149962123 
H  1.4974150642790138  3.265959710877981  -1.9852387995928535 
C  0.8458182912386827  5.015584136313228  -1.1786373147914502 
H  0.9673686513791484  5.537992448468385  -1.962856675264616 
C  0.349519731223374  5.5924788350345676  -0.04018117338152235 
H  0.12453021521987795  6.515906262774674  -0.036506478046668976 
C  0.1746810757661328  4.836300010999906  1.100909866211509 
H  -0.18422804134518134  5.237279003110963  1.884912872289909 
C  0.5204724873165673  3.4925668966804184  1.1079460000369532 
H  0.4107350887572876  2.9793169745631447  1.900167099733281 
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

-Ni 0
lanl2dz