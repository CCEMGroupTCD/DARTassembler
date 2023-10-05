%chk=RIPEWEZU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5785958048214876  1.4142613920347251  0.0 
N  2.9358474195548165  0.45371381387765686  -0.26882824330862964 
H  3.671204358408979  0.7311840058875114  -0.3323904571236541 
N  1.5785958048214876  -1.4142613920347251  0.0 
N  3.9547382250476772  -1.6081506310002944  -0.10715419644189719 
C  1.454667550802559  2.536233853375194  -1.4092919654384906 
C  1.513528197831957  2.0236934389351724  -2.70311650446878 
H  1.621441596029945  1.1122538335978105  -2.8288409822439484 
C  1.352202496269019  2.864213516549518  -3.790370433936235 
H  1.4014039166351413  2.514712304974436  -4.637746243917579 
C  1.1199148108688055  4.209291088993985  -3.6035728019429127 
H  1.0486792574638304  4.707993729094537  -4.326444590497315 
C  1.061505369142483  4.721888702902544  -2.3357285081384647 
H  0.9441169725606496  5.675382577935012  -2.19858418771826 
C  1.219433818325713  3.894588928907674  -1.2367659387637646 
H  1.1380408754775329  4.272882282606793  -0.4243851060248296 
C  1.9430168718364558  2.4463328027518187  1.4328507037008975 
C  3.1793554097246837  3.0886407948359977  1.5040953037574125 
H  3.729425814166317  2.9580817722677937  0.8190727772785668 
C  3.482694695546283  3.8811074324124837  2.5932227049394228 
H  4.324425609844071  4.265653975974388  2.633788500260385 
C  2.565365909281204  4.043370063841487  3.609356213062837 
H  2.7382026544401112  4.62165288666772  4.313966279941672 
C  1.3463684626893373  3.408686659447201  3.5451835445112168 
H  0.7021137084634402  3.5204618842326054  4.228118681095607 
C  1.0311717748635478  2.593768115045827  2.469313741061741 
H  0.2241781391642712  2.201196382460224  2.4079638697040058 
C  2.8315242579489315  -0.9021860189256097  -0.11948507688771327 
C  1.5000035630158952  -2.7635105129872866  0.11077079550455007 
H  0.6237021232870427  -3.1279819592864744  0.1561668477749159 
C  2.606255629932171  -3.5593370952152794  0.10870791819340547 
H  2.5472490446949916  -4.445366429719824  0.17286110436587795 
C  3.860889075486329  -2.9429410139040293  0.02016929249205356 
C  5.1380740211083875  -3.690700986675762  0.07438937914630109 
C  6.348796522616023  -3.0100579383527943  -0.010948274399922284 
H  6.32633037516623  -2.109374449307423  -0.11506731578920219 
C  7.55405697493121  -3.6898725294990333  0.03223857255272862 
H  8.321117371832315  -3.2356205407095198  -0.02503976704168369 
C  7.564844133064951  -5.070251701973866  0.16469191308098302 
H  8.365162062553782  -5.464435830975609  0.1932209190548832 
C  6.368970086201061  -5.7583747253987365  0.2530334586305123 
H  6.391367552946167  -6.724065200081838  0.3317651586611993 
C  5.1600035261012955  -5.086796585495246  0.21277332246209274 
H  4.340198274466459  -5.582019313729635  0.29481269613141464 
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

-Pd 0
lanl2dz