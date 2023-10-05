%chk=REPAMEHA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5848280632295728  -1.4072739640880165  1.8369980547002897e-15 
N  1.584828063229573  1.4072739640880165  -3.94385960489411e-16 
C  2.533361556123682  1.0132205145456807  -0.7850601325770787 
C  2.5136068341395497  -0.40856210281173183  -1.2370057044545297 
C  3.6557126854650175  1.8784025765111558  -1.2477982897668771 
C  3.8436359127661395  2.086787807449498  -2.6124295828887703 
C  4.813937368457836  2.9713492160780635  -3.0396786196900947 
C  5.628356563841757  3.6067043372675855  -2.161800966387245 
C  5.499099208254833  3.346483819198586  -0.790893647126517 
C  4.5065397643503164  2.5108778959680937  -0.3557533233479264 
C  1.5515041587695682  2.7588993821228476  0.49096155319302826 
C  1.4559070580764721  3.851035774339845  -0.3625210741253739 
C  1.3577120533986675  5.116162796819633  0.17184540877541235 
C  1.373594290071233  5.330879559661059  1.5257960785980624 
C  1.4707167990272734  4.230969300891297  2.3846499374746073 
C  1.5539600222886916  2.9561869729438186  1.8637177133476521 
C  2.742341971439701  -1.7504864034036285  1.3578675464857524 
C  2.374123035849483  -1.4496113848890968  2.6548264235156687 
C  3.2304091040783343  -1.7103633915090732  3.7033444016166053 
C  4.4738331697014875  -2.2582457581807684  3.4497147201527363 
C  4.863481924727413  -2.526987140208033  2.1706536710652036 
C  4.009436454139523  -2.291082662137809  1.1226925117683808 
C  1.247825647337053  -2.9991794474327538  -0.7967585228785289 
C  0.4717576945120989  -3.0087988963428143  -1.9533349573612555 
C  0.0979553134028106  -4.2096904804454125  -2.545964346700496 
C  0.5064574407439055  -5.401088774441802  -1.9861269801195365 
C  1.2403746362634664  -5.407562183481173  -0.8301698852859296 
C  1.6124104485492579  -4.209693166489562  -0.23291628545716642 
H  3.4200321510147753  -0.742685125051196  -1.3179939192592391 
H  2.084727183672559  -0.4734518913890245  -2.104422784499488 
H  3.3203310518551095  1.6331040072139476  -3.2323719043596864 
H  4.9110185391440675  3.137502542272461  -3.9508373701242765 
H  6.267482237544371  4.21022898966205  -2.4658629784281643 
H  6.084496098593942  3.7370505560430276  -0.18224097535714728 
H  4.401006236246003  2.36410351139838  0.5562623259038034 
H  1.4590589437657795  3.7277368329624663  -1.2839089113968765 
H  1.277903400113404  5.844762050529604  -0.40074689500981675 
H  1.320889069520296  6.193355908242244  1.8685411265050262 
H  1.4771812496044687  4.357698199848493  3.306144862608688 
H  1.6127109566610731  2.2271748714054755  2.435975307024681 
H  1.5435531318959577  -1.0669653161479844  2.82244304179567 
H  2.973595222633231  -1.5187970036702227  4.576030209538081 
H  5.04763559192299  -2.4464995643056358  4.156475736856579 
H  5.712711615574499  -2.870225878596139  2.007572983509696 
H  4.274827472033515  -2.4931073192089315  0.2544963959301954 
H  0.20295534914049296  -2.20376785127962  -2.3330100363527593 
H  -0.4226876184967119  -4.2083197096863225  -3.315745119681098 
H  0.2819701155587284  -6.205422521022226  -2.395919104181366 
H  1.4910020007776648  -6.21618390820515  -0.44393377814278207 
H  2.108716455053075  -4.221233748078112  0.552964146214757 
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

-Pd 0
lanl2dz

