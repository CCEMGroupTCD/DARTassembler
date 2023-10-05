%chk=HATOPITU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.5899842475634782  -1.2458130246951187  1.5256809330289595e-16 
N  3.7371855678066304  0.8990548784590728  -0.288453580953778 
N  1.5899842475634782  1.2458130246951187  -1.5256809330289595e-16 
C  6.335011151718157  1.1697394782468515  2.7532661321382843 
H  5.8551795692045  1.88922105151051  3.2349580753514693 
H  6.904229155423736  0.6934225685362476  3.409098058611204 
C  5.28764943041345  0.16546518892788736  2.1630330370321427 
H  5.2352556748163215  -0.6265696392851802  2.7536965920886 
H  4.397739369695213  0.5971021475570508  2.156530746830639 
C  5.629042890906114  -0.29373747123052885  0.7434647780053691 
C  5.0286619387106875  0.2541474091406967  -0.389589377850191 
C  5.680820228644046  0.33690393275725716  -1.6084505769496598 
H  5.254268207970801  0.7442056826520126  -2.352803207221252 
C  6.980766526894849  -0.18789162407198745  -1.735157819032929 
C  7.408781534599383  -1.0291743400369218  -0.7074528978622858 
H  8.160132225940309  -1.5929529617702165  -0.8497672130172889 
C  6.764678560970065  -1.0604378837854165  0.5112862211431741 
H  7.102816416199524  -1.6148409292737451  1.2044583455207192 
C  7.935275097168285  0.34637877751005497  -2.76379823592697 
H  7.434148079555091  0.5578919418265755  -3.5907432362377074 
H  8.599291589832392  -0.3546689621487406  -2.983033530050247 
C  8.69784893381006  1.642266878692086  -2.280071833004525 
H  9.671450238897512  1.5031587167410487  -2.393745699880778 
H  8.432230330174574  2.4045787715480915  -2.8547840091527696 
C  8.417597876089614  1.9972329660856774  -0.8329835179915437 
C  7.326773430918676  2.8164209740532096  -0.47604027681769434 
H  6.982865457326848  3.445065533449803  -1.0997772274956572 
C  6.749324784035916  2.707143400050227  0.79247596353295 
H  6.02014410921632  3.2675639335976765  1.0315999043001975 
C  7.227462393395105  1.79839603382562  1.6890172625795232 
C  8.472955183420401  1.2067581031013206  1.441788767208688 
H  8.90659761830213  0.7138579866088985  2.1274959742814716 
C  9.073404269596574  1.3371613672988174  0.19907979339507714 
H  9.93812706278072  0.9723452263604722  0.05241717287190142 
C  2.5481650804547025  0.3058281659828055  -0.022704115738893013 
C  2.1557520829757437  2.4603711766861016  -0.2509226610833934 
H  1.7081124841644548  3.2970014127964316  -0.2896239220711723 
C  3.48677192635005  2.2397635698044356  -0.43577847064814346 
H  4.1361006380299905  2.9052893859299393  -0.6342906925855841 
C  1.9103728221069243  -2.135390873081489  -1.597663660775428 
C  0.8373455360821729  -3.2197664360441625  -1.7671689760187852 
H  0.9464725188609604  -3.8998358792154937  -1.0705753325692984 
H  0.9311088294849683  -3.6375335565375218  -2.6489748094339074 
H  -0.051757613767527744  -2.8151550296709162  -1.690089220180995 
C  1.7418139561137864  -1.0667997493646753  -2.6935324161511436 
H  1.8530516943609086  -1.4821239149846528  -3.5742962568268175 
H  2.4182770592179197  -0.3671929161375972  -2.5743906139123554 
H  0.8478313777872307  -0.6712832711832201  -2.6289059112513984 
C  3.3134659696444184  -2.7365062372295474  -1.7117874020323738 
H  3.42064056638189  -3.15131964803529  -2.592547988133878 
H  3.4359522705000867  -3.4141302371926203  -1.0146474729100023 
H  3.982177749065765  -2.028967820527157  -1.6005788818847912 
C  1.8533532760238653  -2.193651774646363  1.5728130050768105 
C  3.231606809277497  -2.872347857637566  1.6669677447678592 
H  3.310909254802801  -3.3326552541992123  2.5276877415172887 
H  3.9350218623443443  -2.193780112524384  1.5941957143511076 
H  3.3271082186377448  -3.5210258165132737  0.9391093818434189 
C  0.7402038805635679  -3.254014019071673  1.6501214391815557 
H  -0.13246361344068114  -2.818153553466984  1.5642883486441013 
H  0.7905431345387834  -3.7178395800072015  2.512492512752577 
H  0.8555949490512  -3.901536555768076  0.9236033968842186 
C  1.697221513135212  -1.1851199242850434  2.688797825504072 
H  1.7405876726236986  -1.6450875266788578  3.5532765302500557 
H  0.8317210008891158  -0.7343941340827634  2.602327163726295 
H  2.4179337666766476  -0.5229971669781278  2.634562764413327 
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


