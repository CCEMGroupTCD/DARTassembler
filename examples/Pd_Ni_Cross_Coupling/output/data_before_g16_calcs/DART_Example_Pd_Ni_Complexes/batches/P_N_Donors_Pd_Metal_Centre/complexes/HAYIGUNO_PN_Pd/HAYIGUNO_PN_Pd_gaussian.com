%chk=HAYIGUNO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5235755806654299  1.4733694207495955  -3.677796418598735e-15 
N  1.5235755806654294  -1.473369420749596  -4.440892098500626e-16 
C  3.0516315265121343  0.44535836753430225  -0.19700200717187 
H  3.5607555055265427  0.772387999895519  -0.993306957269649 
C  3.8338996733393254  0.789266909997608  1.0642496096280052 
C  4.717386397217659  0.015007858646034933  1.6893544809343752 
H  4.918555058244993  -0.8362985227242691  1.3151087931052832 
C  5.417762337494163  0.4248817174999535  2.971320670859074 
H  5.365118435900774  -0.3180410516200325  3.6234205318835184 
H  6.374120476638915  0.5917655326092044  2.7789910788921004 
C  4.8191738522295475  1.6472661410097778  3.5745164404639147 
H  5.47673243591621  2.041355628563139  4.2011760169645544 
H  4.0262101387803355  1.3805914466211147  4.103512250017793 
C  4.40518667575054  2.7017215532916374  2.588636203710752 
H  3.946703241834604  3.4387048072589685  3.0658536847362994 
H  5.20799702919731  3.076137366617898  2.146373346204132 
C  3.487913969231795  2.1177174889462327  1.561977679807177 
C  2.3358334024667684  2.634198738850853  1.1055995184329392 
C  1.710136339651467  3.9341609296132702  1.480384440985192 
C  0.6135620659327694  3.9740846983379354  2.306416864627946 
H  0.2511232878490699  3.1660319566200528  2.6506114298717987 
C  0.02873896928172237  5.194015595970888  2.643563495374773 
H  -0.7355531183387385  5.212845652852268  3.207394271191979 
C  0.5519685934854709  6.364508871935991  2.163824387863108 
H  0.15774700670839548  7.195890447221298  2.4037871282939753 
C  1.6253012943952734  6.337355023460134  1.3559349637024631 
H  1.973133860378709  7.152003857008844  1.0153625719278068 
C  2.238356612981473  5.129705763918183  1.0072938555683733 
H  3.0081102516692377  5.12713657797833  0.45258588283806334 
C  2.7399965188436055  -1.0074124989270152  -0.38636813813691256 
C  3.7004804545536016  -1.8657038019707746  -0.8881358417234251 
H  4.539130339776821  -1.5233489024484448  -1.1747997089759048 
C  3.4455843157923756  -3.2164021685674715  -0.9753919069431708 
H  4.099032388817592  -3.8097765169744084  -1.325824366216389 
C  2.2356157576029454  -3.6882856420593235  -0.547927749804513 
H  2.0428833142569025  -4.619449687327755  -0.5736822849623189 
C  1.299433121034693  -2.7915469685175527  -0.07999400196975089 
H  0.45534188366272565  -3.1235005933940423  0.19787809796694422 
C  1.2127109043493733  2.2978208874985078  -1.5719023481102545 
C  2.182403417115483  3.100689925867834  -2.1569083975774945 
H  3.005545984781921  3.2522442997100685  -1.7059944977420523 
C  1.9584470896317807  3.679441675423061  -3.3875902922061494 
H  2.627517792034297  4.221739865372934  -3.785909169945735 
C  0.7617106562073628  3.4703317905150026  -4.039571261980871 
H  0.6161256273837368  3.8452093006773858  -4.899442363807706 
C  -0.21825275089826301  2.7169987198209578  -3.438034434070364 
H  -1.0588230810883037  2.6118070633244157  -3.871057511293491 
C  0.0017276491567106422  2.1052863432050084  -2.212725533471588 
H  -0.6720266426049906  1.5622692636235969  -1.8176940845772283 
Br  -1.654629867976521  1.654629867976521  0.0 
C  -1.352040782177733  -1.3199206589279684  -0.4819701001757076 
C  -1.6380700204479186  -1.599217458033982  -1.8139382277462257 
H  -1.1579680913617836  -1.159339363685303  -2.503515938957912 
C  -2.63068216563725  -2.52632649780156  -2.142289956238523 
H  -2.8368222554424674  -2.694399528648736  -3.0538898079700485 
C  -3.3101544432602585  -3.1961573187854975  -1.1619874135498987 
H  -3.9960123560788676  -3.811990145619169  -1.3868663427944008 
C  -2.9974710358883443  -2.9701089012107795  0.12423560353005347 
H  -3.438659168087169  -3.468581175155742  0.8017084178296768 
C  -2.0389854250767776  -2.015940632791841  0.4912562833419403 
H  -1.8619130127668426  -1.8492165108145995  1.4101568590026907 

-Pd 0
lanl2dz
F 1 1.0
1.472 1.0
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

-Pd 0
lanl2dz

