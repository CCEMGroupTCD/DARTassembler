%chk=EDAFEPAP_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5023398700027897  -1.4950166938867278  3.912982468754431e-18 
C  1.7960979874267944  1.6301113994013023  1.4288874234835605 
H  0.9599836880917341  1.7331520447866786  1.9499388545784422 
H  2.373817289150898  2.4169461985109946  1.5943371056392979 
C  2.533121769432068  0.3110780246814  1.8221448469782429 
H  3.362619836614553  0.5023701160469484  2.326958458444424 
H  1.9498515657697966  -0.2742129380068718  2.369199943394695 
C  2.853314705601366  -0.34661410278490384  0.4523662922124896 
H  3.757807755576361  -0.7720583353482291  0.42138988949152967 
C  2.7138565020490852  0.8194572441963799  -0.49493367253534803 
H  3.5026178374903516  1.4163302782225866  -0.4476765358758684 
H  2.5962658811162473  0.5114083328476388  -1.4283166400709166 
N  1.5023398700027897  1.4950166938867273  1.8141182329285736e-17 
C  1.2166430712857803  2.790181274936354  -0.7130614960839029 
H  1.0044207338061912  2.583020441036351  -1.6579680937815706 
H  0.40650560671926184  3.192736311586399  -0.3112677418101046 
C  3.31931092370929  3.8303116436195928  -1.6533394169815614 
H  3.3861763375434997  3.125657331525457  -2.2883083980191605 
C  4.2323444678643805  4.879802988041249  -1.677500857990351 
H  4.925966589143492  4.884660251528469  -2.327562678022681 
C  4.140448860417217  5.921654096707207  -0.7617968421077875 
H  4.749513425938388  6.6490824822543155  -0.7963323155639742 
C  3.156326708522311  5.882115969647015  0.1946780486542701 
H  3.112593279124515  6.573368963159274  0.8462917272349592 
C  2.2444248924719807  4.87717492374199  0.22964134568765032 
H  1.556505202666753  4.889416355610167  0.885507206763008 
C  2.3050472650958316  3.8140511362290037  -0.6936783066329976 
C  3.115326281170587  -2.1019133818702973  -2.207088613423277 
H  3.636590180430316  -1.351658100717575  -1.9476707643245634 
C  3.4636823754570383  -2.803301673578195  -3.3278753174078557 
H  4.217088681351908  -2.5347321842595556  -3.8415270453763424 
C  2.7328583717505914  -3.88797574887843  -3.707938898901886 
H  2.961253815418116  -4.357161703703056  -4.501780917887552 
C  1.6358354252882126  -4.31779120018916  -2.9276925834670773 
H  1.137340599919251  -5.089400955262171  -3.1710363628339344 
C  1.3058014330522172  -3.580007542314118  -1.7923907874052294 
H  0.5754573223507954  -3.8517385113803275  -1.2479059571169049 
C  2.0314710389258948  -2.457233277074255  -1.4531351427561823 
C  0.32617241806418407  -2.758160468616428  2.195235146712509 
H  -0.3895123762496471  -2.1445973801566955  2.0807244412403048 
C  0.26349461340280556  -3.6815713413144135  3.243092500070297 
H  -0.4971462821768975  -3.726010570837183  3.8120039543856135 
C  1.3442615962514985  -4.534338447670618  3.429461840630971 
H  1.3309635371683457  -5.161718401911492  4.142284156949266 
C  2.43579928641526  -4.478354133018985  2.591590849687274 
H  3.171529357504999  -5.063178756433019  2.7299768039706223 
C  2.455895969279168  -3.581713884464298  1.5608971716191187 
H  3.214051985056344  -3.5491692278066465  0.9893240432585052 
C  1.3872776391660473  -2.711654787252267  1.3266599160322818 
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


