%chk=UTUSAKIG_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3547405563428008  1.4982583305291506  0.0 
N  1.3547405563428008  -1.4982583305291506  0.0 
N  2.643231283287524  -0.8711426850795698  0.14710830639307432 
C  1.2066223041856057  -2.5069000708799463  1.089090713852635 
H  1.8833420304208075  -3.208508330670994  0.97949380734 
H  0.31250128476490335  -2.9066884103103385  1.0472617013661751 
H  1.3282506918191133  -2.0696377499231717  1.957976033094007 
C  1.2999967913978634  -2.1645154302960856  -1.3331801593220094 
H  1.4486722203078413  -1.499491559113535  -2.0371328335805994 
H  0.42129715332458006  -2.57965688976678  -1.455467749364654 
H  1.996560743611942  -2.853374380031066  -1.3833586193016494 
C  1.5103906365796282  2.950358879992886  1.0373238967503877 
C  1.628078557050499  4.017032745530859  0.2160157950287881 
C  1.5800824335073078  3.7027196458441844  -1.2482526938386806 
C  1.4341687670487264  2.403754133459783  -1.5536493041285486 
C  2.81317151273921  0.4965946324738444  0.18633939847192799 
C  4.159518091883143  0.7223139981593074  0.3319520585739386 
H  4.5815677366497765  1.5716931766822468  0.3868778124975 
C  4.799978703502813  -0.5195327095061693  0.384419664048125 
H  5.734073905598084  -0.6609150170052951  0.4873091025189189 
C  3.8537166467292776  -1.4909432835998744  0.2616719697737844 
H  4.007517593772348  -2.427652879270428  0.2552679033029368 
C  1.487619787062285  2.870642757234659  2.5280093760466364 
H  1.5513505103663539  3.7725700700784603  2.9050552847282 
H  2.245556770556406  2.3321929294836066  2.837359116707111 
H  0.64994545278452  2.4533280043247396  2.819339509565398 
C  1.8169737838858477  5.429629543342436  0.6543981054512191 
H  1.8381980769294994  5.46873658103153  1.6330854656372982 
H  1.0755234422684459  5.976719953786212  0.3213000812662667 
H  2.6616070663660203  5.77268900413547  0.2952403622697896 
C  1.7443535121865588  4.795025933520227  -2.264159326661283 
H  1.6449173124262662  4.422007747262174  -3.1649631037948045 
H  2.6336843400305874  5.1962710911478025  -2.173202402094231 
H  1.0614496718056807  5.482026021573446  -2.1185521709033077 
C  1.3753293739259798  1.7275026190482707  -2.880581463311203 
H  1.237680355858188  2.397569746020583  -3.582567761735936 
H  0.6334369038457731  1.0884576723321402  -2.889391623064606 
H  2.2165444532533956  1.253578910077989  -3.042899941581643 
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
