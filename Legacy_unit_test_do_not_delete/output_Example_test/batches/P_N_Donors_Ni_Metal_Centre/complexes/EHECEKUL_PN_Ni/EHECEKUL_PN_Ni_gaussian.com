%chk=EHECEKUL_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.5899842475634782  1.2458130246951182  0.0 
N  3.756334303283685  -0.8990548784590728  0.024625476646132627 
N  1.5899842475634782  -1.2458130246951187  0.0 
C  5.964103704973419  -1.169739478246851  -3.3110170171473308 
H  5.4291452334705825  -1.8892210515105092  -3.730641741754216 
H  6.449153037732959  -0.6934225685362473  -4.0313306962575295 
C  4.996480178872574  -0.1654651889278873  -2.5975421447452938 
H  4.872493178756527  0.6265696392851803  -3.177417790638211 
H  4.113995802637279  -0.5971021475570507  -2.4826355668395705 
C  5.5083307959264625  0.29373747123052873  -1.2301605291291506 
C  5.050509558322729  -0.2541474091406969  -0.03238395469649402 
C  5.846348573293333  -0.3369039327572575  1.0979139352642089 
H  5.513689773012367  -0.744205682652013  1.8887018885919606 
C  7.152046997359914  0.187891624071987  1.0652531189492662 
C  7.451625922144475  1.0291743400369213  -0.006953353538004414 
H  8.21471991151305  1.592952961770216  0.04273355687184199 
C  6.663797058192602  1.060437883785416  -1.138111769480101 
H  6.914938053460209  1.614840929273745  -1.8673257323533121 
C  8.224800537700347  -0.3463787775100554  1.9699008734263619 
H  7.828188089104918  -0.5578919418265761  2.851753971702916 
H  8.910585613466512  0.35466896214874  2.1065787669412104 
C  8.922738845489608  -1.6422668786920862  1.3968457209889125 
H  9.902936433605658  -1.503158716741049  1.3910201283097232 
H  8.729139918200994  -2.404578771548092  1.9996448445475106 
C  8.468221033855874  -1.9972329660856771  -0.005302204769358387 
C  7.342026990250418  -2.8164209740532096  -0.22664678624757678 
H  7.076696867179737  -3.445065533449802  0.434352762717759 
C  6.614289316706148  -2.707143400050227  -1.4153344114563344 
H  5.861401966408979  -3.2675639335976756  -1.5638111885661872 
C  6.979602061380189  -1.7983960338256195  -2.3634633441224335 
C  8.245940711418696  -1.2067581031013204  -2.269865041002264 
H  8.59278415425479  -0.7138579866088984  -3.00330880934672 
C  8.993362267841754  -1.3371613672988174  -1.1095953673139392 
H  9.869513225718451  -0.9723452263604722  -1.0694091467568463 
C  2.543789881503833  -0.3058281659828056  -0.09423798626479318 
C  2.1821147152760685  -2.4603711766861016  0.18010256699647736 
H  1.7425282514191904  -3.2970014127964307  0.27306889864283895 
C  3.5257415947117545  -2.2397635698044356  0.20136997501726253 
H  4.194422862785309  -2.9052893859299393  0.3192692532084927 
C  2.10269091557739  2.135390873081489  1.546709372905717 
C  1.058319313699781  3.2197664360441616  1.8457203520369039 
H  1.0817394705497074  3.899835879215493  1.141019778234273 
H  1.2588488076610713  3.6375335565375204  2.7095264674854493 
H  0.16644974467050244  2.8151550296709154  1.8775695541022839 
C  2.06894126748673  1.0667997493646746  2.6549518470815485 
H  2.2866879674593354  1.482123914984652  3.5155941375269686 
H  2.725842383989137  0.36719291613759664  2.4542579955567185 
H  1.1737463100156187  0.6712832711832195  2.6997561281855753 
C  3.509233804973962  2.7365062372295474  1.488988412474553 
H  3.722947252844773  3.1513196480352885  2.350122644992849 
H  3.545847126093382  3.4141302371926194  0.7821175336072336 
H  4.159408198864396  2.028967820527157  1.2971133582869763 
C  1.6597124749955618  2.193651774646363  -1.59318610601213 
C  3.0162181391078784  2.872347857637565  -1.854605883694715 
H  2.990034095370231  3.332655254199212  -2.7185747410756833 
H  3.723258722601102  2.193780112524384  -1.8681010156159965 
H  3.1997113162262614  3.521025816513273  -1.1438115606871087 
C  0.5454387982255542  3.254014019071672  -1.5342595088356283 
H  -0.3102635422246054  2.818153553466983  -1.3427147906338837 
H  0.49030623452745736  3.717839580007201  -2.3963374110248052 
H  0.7485100359735032  3.9015365557680752  -0.8272194556068253 
C  1.3687401572758722  1.1851199242850432  -2.681824869389412 
H  1.3064296198525083  1.6450875266788578  -3.5451448863437776 
H  0.5202290773770828  0.7343941340827635  -2.4905207676677508 
H  2.090689922253836  0.5229971669781279  -2.7158267973429 
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
