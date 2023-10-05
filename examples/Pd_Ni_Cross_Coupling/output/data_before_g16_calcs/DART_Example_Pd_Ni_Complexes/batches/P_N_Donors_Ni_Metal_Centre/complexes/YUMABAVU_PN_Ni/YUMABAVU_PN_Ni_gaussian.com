%chk=YUMABAVU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.2795876122798324  -1.5629317139593777  2.092946912165683e-15 
N  1.2795876122798315  1.5629317139593777  -4.1344853700365514e-16 
C  1.1334294628717854  -3.0690578319307615  1.0185523904670246 
H  0.3221622212900439  -3.5551084816264806  0.7580884031928923 
H  1.9005022196288306  -3.6517017963897667  0.8421312564703938 
C  1.0713418867528632  -2.76286894618112  2.5202434835800847 
H  1.8020549309557308  -2.1714826292073113  2.7594227701580434 
H  1.1482001364082122  -3.590217923395366  3.0217675102851653 
H  0.2271976443010857  -2.3349485392775446  2.7284020848701105 
C  2.9026807384261053  -0.8850209834157544  0.4987984214796517 
H  3.4192369773462277  -1.5895851688914997  0.9440251733461565 
H  3.396342693265251  -0.616952464325681  -0.3060677807038313 
C  2.7898382015347996  0.33573399132755605  1.4502202961006947 
H  2.100539474580947  0.16029135838913014  2.1238731277339102 
H  3.6423637007246707  0.46764622918067933  1.9147213668590386 
C  2.441363450969936  1.5645618363615554  0.6946083020932765 
C  3.2923657526898737  2.6436509289537575  0.6311182027397947 
H  4.092494779756845  2.6306471170671117  1.125585725355087 
C  2.99776230748265  3.73296019934217  -0.13833574878551413 
H  3.570151696686565  4.4786434589438056  -0.15119823208604244 
C  1.8489578135441735  3.7210895367698433  -0.8940706044570652 
H  1.6383593747032887  4.442051588365456  -1.4597225802485696 
C  1.002867095790924  2.617766282402239  -0.80321955215999 
H  0.2156221855832976  2.604616421268553  -1.3159910339180616 
C  1.4861897887239217  -2.1728795712365097  -1.6941054312207564 
C  1.630638804965188  -1.2309901803873708  -2.6919403259366987 
H  1.592329604256799  -0.3159035472184624  -2.4801357014907923 
C  1.8350171254912913  -1.6303095104807777  -4.012058473949489 
H  1.96665213278175  -0.9833998560348777  -4.682116392222482 
C  1.8455864245429343  -2.9403982857688944  -4.334584350092266 
H  1.967099972716249  -3.2095640536924788  -5.2270536588017364 
C  1.6809623643565055  -3.8648710487789075  -3.3617891973183256 
H  1.6953992222763274  -4.777240755826173  -3.589711190105089 
C  1.4897340937106724  -3.506675873007748  -2.0450012323975897 
H  1.36384812278758  -4.168529980747023  -1.3893815237299807 
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


