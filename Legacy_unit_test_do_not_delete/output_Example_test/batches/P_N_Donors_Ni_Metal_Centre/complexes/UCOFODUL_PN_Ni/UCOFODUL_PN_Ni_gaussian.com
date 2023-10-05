%chk=UCOFODUL_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3397501455495344  1.5116777260712677  1.1511244570325938e-15 
N  1.3397501455495346  -1.5116777260712677  0.0 
N  3.3318743208337342  -3.4579285533863584  0.23306972652761943 
N  3.154974789156391  -2.5284586103108335  1.1649953891271272 
C  3.468013548788845  -0.32298456481787263  0.027663058792358003 
C  3.0135768817977038  0.9185634633307651  -0.44459429635865266 
C  3.852704646817507  1.7062982790674486  -1.2278623884249582 
C  5.134264453968051  1.293600182614457  -1.5013228094741704 
C  5.5858031869884535  0.0712971988944413  -1.0248322234480922 
C  4.7573676708689  -0.7233566577327788  -0.28037782338100636 
C  2.5294282448753522  -1.2373337783848628  0.8099102594438434 
C  1.4718745634277375  -2.4445964796373385  -0.9247733978900418 
C  2.5995617684684067  -3.3552870397065537  -0.8626876395519976 
C  2.812079434289872  -4.340225253931326  -1.921101006754777 
C  3.9066505461199017  -5.19344297868544  -1.9188975293628354 
C  4.06488387914284  -6.122964063613931  -2.9759496496453868 
C  3.1574984494928025  -6.173181910110414  -3.9894422988304323 
C  2.099830941016194  -5.317545948362427  -4.003366108528858 
C  1.8930405956415268  -4.386729190443443  -2.962554609631458 
C  0.7319143147994079  -3.528132478581322  -2.9492950350780704 
C  0.5036575004745562  -2.617305136712366  -1.9704979212601714 
C  3.847555642273519  -2.6960421571415396  2.3921706931188793 
C  4.128844213186776  -1.601835626561213  3.206870507964171 
C  4.809799984556115  -1.8022788955028592  4.383762254631988 
C  5.205215361806788  -3.052966857885965  4.769163643035951 
C  4.933697844517011  -4.138682147391145  3.943337942334217 
C  4.258191465167277  -3.9752599174148915  2.738983369873356 
C  1.003381942935663  2.8257950248442114  -1.1859212947476905 
C  0.6041643915906334  2.4486679856824765  -2.4652560722907784 
C  0.4953788171460989  3.4624829599779616  -3.4680923833042256 
C  0.6858310793018535  4.762801133685926  -3.158799766249534 
C  1.0283349746296935  5.122368361955304  -1.9258699912663921 
C  1.1968400960572407  4.154497184940274  -0.9298859207329983 
C  1.6052219828755232  2.374621836116652  1.5790568281492612 
C  0.5577427650094181  2.5939896446652613  2.422625552505033 
C  0.7594435201999776  3.3726174233909694  3.5994617380263043 
C  2.0155624133083  3.887173278508018  3.866009760605934 
C  3.0437896474484  3.6464090212043816  3.0364125466655207 
C  2.8464244679503636  2.886201253251981  1.8949040554721601 
H  3.53846861336343  2.5353389822450216  -1.5717156969067565 
H  5.709258938099811  1.8433980004072783  -2.020306405330843 
H  6.473116442220618  -0.21401061345316652  -1.214123399425098 
H  5.073260710227208  -1.5647281505214718  0.030223929362753088 
H  2.264009795379402  -0.7959181555363317  1.6083031856226895 
H  4.542881317030327  -5.1547128849688635  -1.213998014547108 
H  4.808408065847898  -6.715045202711036  -2.9782406123413434 
H  3.2646328895407413  -6.806692104435292  -4.6879063604986495 
H  1.488094452942565  -5.346069835601319  -4.728183290317579 
H  0.10336626972678675  -3.6038302268688467  -3.6568076398701033 
H  -0.2923017670656334  -2.097807844759163  -1.9771656097798835 
H  3.854085481992477  -0.7280505087286795  2.9531125376718492 
H  5.00842974257286  -1.0587299659825418  4.941013632951855 
H  5.663206815860619  -3.179985025639895  5.59343932693917 
H  5.215786483599057  -5.006682609639233  4.206161242743362 
H  4.085787325025177  -4.716563748696697  2.1731246399868267 
H  0.4099734678469831  1.54081030124948  -2.6689943110962293 
H  0.2855695674825809  3.218427872284881  -4.362693329495738 
H  0.5745038719569966  5.426456939879039  -3.829987909627345 
H  1.1580134098362986  6.039828093008199  -1.7205194486844084 
H  1.4508763064966945  4.426746274744352  -0.0565519485972342 
H  -0.29978555649615957  2.2310089308949697  2.2257592546105838 
H  0.03739469888928815  3.538815053128457  4.194427970272741 
H  2.1502630282402664  4.411567834522418  4.646266509349181 
H  3.9027060234013358  3.9932925903384007  3.234892919964081 
H  3.58316280558848  2.716407860405562  1.3195209434503106 
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
