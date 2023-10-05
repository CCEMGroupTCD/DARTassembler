%chk=SITEGUZI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.6930402793495496  1.2750351416725736  2.8971276275900483e-15 
N  3.8495058551702748  -0.8306112575930586  0.20412486756878806 
N  1.693040279349549  -1.2750351416725736  1.3877787807814457e-15 
C  5.118844443613681  -0.14038283859784406  0.39123330537943707 
H  4.975672455722062  0.8270546066022075  0.3385363272121887 
H  5.746793217203372  -0.41675335458044693  -0.307667113595067 
H  5.487703978881291  -0.36787577233730184  1.2696477055009185 
C  3.696692370207865  -2.201520107765654  0.23213510377086666 
H  4.390646003703806  -2.844695455014704  0.3218511405138285 
C  2.375089094787377  -2.459697011476856  0.10516243431796893 
H  1.9815602272483026  -3.3242202383065287  0.09121644258464087 
C  2.6098828299087504  -0.30313204308470265  0.06420026821918595 
C  2.1320113305925705  2.29090846364304  1.5104086797550327 
C  3.421684873279622  3.1024513664077698  1.3957536743462498 
H  3.567491114047421  3.602232065913888  2.2256675569596065 
H  3.346873714439483  3.728215499428733  0.6459018820982485 
H  4.1769646683268515  2.495825388252083  1.242302822865785 
C  0.9559349583348833  3.248533940065832  1.7834667752644862 
H  1.156838919632506  3.797681624012417  2.5702278823046525 
H  0.14095336283230275  2.7281754665498297  1.9482295988388572 
H  0.8207283300278743  3.8299565971916  1.0071720107306315 
C  2.2263529375095983  1.2956790628789134  2.673863926448706 
H  2.453572204961321  1.776048638572763  3.4973889626476335 
H  2.921033934545602  0.6325092196612596  2.480324211134908 
H  1.3644752847923003  0.8423017432219098  2.7865235366343444 
C  2.0696884673971696  2.0844406719752184  -1.637690327856888 
C  3.5571696666700467  2.283093835893303  -1.9291114188491951 
H  3.665172713289093  2.706682941281814  -2.8064772671229856 
H  4.01001329194796  1.4139376426402515  -1.9295389408135535 
H  3.951466003190876  2.856272100643966  -1.2397343987624407 
C  1.4909515497552388  1.1412832595273543  -2.7095127457522317 
H  1.6630823502233234  1.513788757200898  -3.6000742830607098 
H  0.5247251608341703  1.0489634903282181  -2.576096533965123 
H  1.9168437382945915  0.26129461974119583  -2.6343378989717277 
C  1.346022356028136  3.425333388671056  -1.6892071364545576 
H  1.5338951208215479  3.86456845105949  -2.5446973784898153 
H  1.6579311607924416  3.9938362764025572  -0.9548398345194435 
H  0.38147486178933443  3.2788965659533136  -1.6025745311029345 
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

