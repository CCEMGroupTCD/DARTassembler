%chk=LAFOZOSU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3547405563428008  -1.4982583305291506  1.8348372687783814e-16 
N  1.3547405563428008  1.4982583305291506  -1.8348372687783814e-16 
N  2.643231283287524  0.8711426850795698  -0.14710830639307443 
C  1.2066223041856057  2.5069000708799463  -1.0890907138526351 
H  1.8833420304208075  3.208508330670994  -0.9794938073400005 
H  0.31250128476490335  2.9066884103103385  -1.0472617013661756 
H  1.3282506918191133  2.0696377499231713  -1.9579760330940073 
C  1.2999967913978634  2.1645154302960856  1.3331801593220092 
H  1.4486722203078413  1.4994915591135352  2.0371328335805994 
H  0.42129715332458006  2.57965688976678  1.4554677493646537 
H  1.996560743611942  2.853374380031066  1.383358619301649 
C  1.5103906365796282  -2.950358879992886  -1.0373238967503873 
C  1.628078557050499  -4.017032745530859  -0.2160157950287876 
C  1.5800824335073078  -3.7027196458441844  1.248252693838681 
C  1.4341687670487264  -2.403754133459783  1.5536493041285488 
C  2.81317151273921  -0.4965946324738444  -0.18633939847192793 
C  4.159518091883143  -0.7223139981593074  -0.3319520585739385 
H  4.5815677366497765  -1.5716931766822468  -0.3868778124974998 
C  4.799978703502813  0.5195327095061693  -0.38441966404812505 
H  5.734073905598084  0.660915017005295  -0.48730910251891896 
C  3.8537166467292776  1.4909432835998744  -0.26167196977378454 
H  4.007517593772348  2.427652879270428  -0.2552679033029371 
C  1.487619787062285  -2.8706427572346596  -2.528009376046636 
H  1.5513505103663539  -3.772570070078461  -2.9050552847281996 
H  2.245556770556406  -2.332192929483607  -2.8373591167071104 
H  0.64994545278452  -2.45332800432474  -2.8193395095653977 
C  1.8169737838858477  -5.429629543342436  -0.6543981054512185 
H  1.8381980769294994  -5.46873658103153  -1.6330854656372975 
H  1.0755234422684459  -5.976719953786212  -0.32130008126626597 
H  2.6616070663660203  -5.77268900413547  -0.2952403622697889 
C  1.7443535121865588  -4.795025933520227  2.2641593266612836 
H  1.6449173124262662  -4.422007747262174  3.164963103794805 
H  2.6336843400305874  -5.1962710911478025  2.1732024020942315 
H  1.0614496718056807  -5.482026021573446  2.1185521709033086 
C  1.3753293739259798  -1.7275026190482703  2.880581463311203 
H  1.237680355858188  -2.3975697460205825  3.5825677617359366 
H  0.6334369038457731  -1.0884576723321397  2.889391623064606 
H  2.2165444532533956  -1.2535789100779886  3.042899941581643 
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
