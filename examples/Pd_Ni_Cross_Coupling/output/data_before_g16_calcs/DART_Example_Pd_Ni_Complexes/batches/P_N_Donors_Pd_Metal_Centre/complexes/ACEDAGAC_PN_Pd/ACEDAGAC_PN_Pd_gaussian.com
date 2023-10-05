%chk=ACEDAGAC_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5786548854958768  -1.4141954435296422  1.0424873061311846e-16 
N  1.5786548854958768  1.4141954435296422  -1.731889923287348e-16 
H  1.3147072006828435  2.0819702651221887  -0.57482973768664 
H  1.7189558913265621  1.8014633013817956  0.8231281334260898 
C  1.7084679742748907  -2.816280190070151  1.1218338092160176 
C  2.8687313708857576  0.878783136507787  -0.48062341113774204 
H  3.60066115128422  1.5073070997738092  -0.2594424803000261 
H  2.8425030415719768  0.7611745020287481  -1.4624020337596642 
C  3.109669074892956  -0.4551561882376469  0.19623801536361368 
H  3.3154398391048363  -0.32476549800902604  1.1553082007501585 
H  3.8694175824678623  -0.9260730993668704  -0.2297961124363056 
C  1.6337431642404643  -2.0663909687096584  -1.6777625632014124 
H  2.4526697405940094  -2.5901757421573994  -1.798427610527385 
H  1.623952311622743  -1.3238075762904706  -2.3172271396309707 
H  0.8538472788215772  -2.639151207547907  -1.831485582663866 
C  1.7616241179914833  -2.5790251476111146  2.495812138453053 
H  1.6830430496487845  -1.6922989028495175  2.8281757383690493 
C  1.7952873118821366  -4.12664012898701  0.6442770135369597 
H  1.7418358801809999  -4.298587935856677  -0.28891563332746656 
C  1.9303958292810637  -3.6439066722194458  3.376311257144401 
H  1.9723472094439345  -3.4831859698300267  4.3118852774794645 
C  1.9598049537183346  -5.182715399087183  1.539967668510836 
H  2.0180505969400038  -6.073450060808355  1.2170350324175387 
C  2.0372604380592425  -4.932968122912623  2.8926789776928556 
H  2.164991415464799  -5.654371278090174  3.498495010091224 
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


