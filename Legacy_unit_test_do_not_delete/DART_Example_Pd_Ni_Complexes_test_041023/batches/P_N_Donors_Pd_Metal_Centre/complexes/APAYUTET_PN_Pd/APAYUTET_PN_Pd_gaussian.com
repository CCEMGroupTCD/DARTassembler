%chk=APAYUTET_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
C  1.7262952681381998  3.157485371747408  0.12348807142408363 
C  1.546382756629089  3.7511567321516113  1.3832652001410426 
C  1.6720105819513094  5.129159538098731  1.4951854230977113 
H  1.5310906756914127  5.5442192936733115  2.338379629269738 
C  2.0040990083567385  5.92084050643037  0.3830965680914787 
C  2.1592132554728662  5.316609641855731  -0.841049212703038 
H  2.3681042276270396  5.855751453240702  -1.5947646291362059 
C  2.0222975487497608  3.9317550276998006  -1.019238934358204 
C  1.2452963155899703  2.9349679082036064  2.594349317238119 
H  1.0342864672397623  3.5300103694637577  3.34499129128214 
H  2.0264903407017174  2.387626898466016  2.8211183499325374 
H  0.4790778420249773  2.351973666097387  2.4134061647323173 
C  2.1544153184875228  7.427820090517111  0.5508146352590148 
H  2.008812875573138  7.66844420293228  1.4898650892505934 
H  1.4955091221361438  7.886640564199139  -0.009359888742183587 
H  3.0576198488202815  7.697705844143776  0.2806610745277549 
C  2.149981655524741  3.3457004098463856  -2.395694634251879 
H  2.3051470772286535  4.064282072393226  -3.0439404102958747 
H  1.324408442007293  2.870147667097044  -2.6258190744708374 
H  2.9038475754199893  2.718927701963377  -2.4155701838920693 
C  3.045779332447819  0.5625626703139153  -0.26110676423823187 
C  4.379765535034146  1.1975690961835377  -0.3184204660349634 
C  4.893239624814411  1.7910559250778695  0.8353670070523957 
H  4.370452238338857  1.8013960335247161  1.6283638549019923 
C  6.158061691669545  2.3662182509181147  0.8422494579459552 
H  6.504671863623708  2.743325099828982  1.642486712266459 
C  6.910354722199894  2.3876727783494225  -0.31554618441420296 
H  7.767846218257279  2.794445823203058  -0.31811604235719904 
C  6.402420364151814  1.8063711916428973  -1.4888713453166413 
H  6.914570613914705  1.8256162630050687  -2.2883318189462525 
C  5.157951137587158  1.205398629408686  -1.4826445022018573 
H  4.827388534666612  0.7950022235739863  -2.2738152837834082 
C  2.8599708441032465  -0.8936258733968203  -0.2965976197465764 
C  3.912767797700927  -1.768184576502025  -0.563963285676988 
H  4.765539073488188  -1.427319273616916  -0.8071307310335065 
C  3.7116418253211894  -3.1440073244711626  -0.4741915228283293 
H  4.425495772400197  -3.748364725974614  -0.6437889682328352 
C  2.4682389281062016  -3.610963684515544  -0.1378464425736639 
H  2.3139688623251393  -4.543988136111363  -0.05236568726450255 
C  1.4435102314532235  -2.716400459964525  0.07711233379238608 
H  0.5821551319109004  -3.055140268164091  0.2889893614130667 
N  1.6099263787515259  -1.3784908614133058  -2.577521681806882e-15 
P  1.6099263787515254  1.3784908614133058  -6.129056519584308e-16 
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

