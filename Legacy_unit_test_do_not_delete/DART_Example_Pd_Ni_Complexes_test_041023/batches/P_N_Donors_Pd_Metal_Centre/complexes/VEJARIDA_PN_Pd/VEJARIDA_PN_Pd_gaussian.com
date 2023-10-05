%chk=VEJARIDA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5557815511825557  1.4393205914597345  -1.1630042485209544e-15 
N  1.555781551182556  -1.4393205914597345  -2.220446049250313e-16 
N  2.950525699860518  0.49098465567236316  0.042692535915760596 
C  1.5402418686750512  -2.7872355323780753  0.027677072039335344 
H  0.7068968286898917  -3.2414985602813293  0.019039606420349277 
C  2.7206283510956117  -3.529841398883462  0.0678438387828267 
H  2.669992996166746  -4.478154052852869  0.07993065609443209 
C  3.9676031433540055  -2.9195659152893856  0.08972674526648985 
H  4.7668182078449295  -3.433039629474538  0.09579267731903829 
C  4.0049447459148935  -1.5239688975584023  0.10324039830961679 
C  4.999211242336565  -0.47351488198247305  0.0990997155404676 
H  5.941244574916884  -0.5945581269924745  0.11948483621956824 
C  4.348127622199006  0.7170061042157125  0.061234150573476816 
H  4.766486583587572  1.5704737731206397  0.04868036389910167 
C  2.7659823959158243  -0.8740967977982246  0.06535096392403802 
N  1.5270110974574909  2.475746095861845  1.3068170204039027 
C  2.266850833289965  3.649603567244561  1.4987386506656513 
H  2.8093693836311933  4.0816365616237675  0.8492206722359548 
C  2.070741173912029  4.059184072448145  2.778762439997652 
H  2.4432650265591818  4.835626330435021  3.1808878050533522 
C  1.2080219239121515  3.121028755966123  3.41112299999185 
H  0.9076564175945852  3.1606533784014617  4.31245844360401 
C  0.8875823848423837  2.1607807311683307  2.509322105851302 
H  0.3295033234122484  1.4088571802779826  2.666019313260086 
N  1.8572408988913918  2.4713192745978367  -1.2826990292661078 
C  1.1506330557783346  3.65460386741584  -1.5330665360449771 
H  0.7181908768054136  4.191783040414304  -0.8792083878646144 
C  1.187831413291679  3.903740656069942  -2.8661750210293766 
H  0.8045894743899136  4.653613656153533  -3.305992767388555 
C  1.9043161796602692  2.8375849996406157  -3.4919129180961415 
H  2.0770727968543934  2.7491939276620134  -4.422645892482516 
C  2.293740035528138  1.9744763132428953  -2.5212813783521266 
H  2.779985341649314  1.1701497915012058  -2.6563768952948914 
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
-P 0
6-31+g(d)
****
-N 0
6-31+g(d)
****
-C 0
6-31g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****

-Pd 0
lanl2dz


