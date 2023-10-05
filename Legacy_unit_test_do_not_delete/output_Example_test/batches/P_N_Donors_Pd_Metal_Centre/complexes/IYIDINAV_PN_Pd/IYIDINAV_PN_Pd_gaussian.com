%chk=IYIDINAV_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5818420022239887  -1.4106296041129998  1.0534872228452487e-17 
N  1.5818420022239892  1.4106296041130002  4.9292301979082937e-17 
C  1.6528522750360373  2.6954351084889145  0.32803856035611667 
H  1.1412935121350873  3.00392954156416  1.065476474258602 
C  2.4482661666463508  3.6164565419054817  -0.36439690972080596 
H  2.4846032976078476  4.5239992279991785  -0.08391798746224277 
C  3.1701893846394253  3.202524693107259  -1.4421046539043032 
H  3.659746562956143  3.8296603759403633  -1.960853778699749 
C  3.182858668996012  1.8326928574745163  -1.7800570785455938 
C  3.9729792138823923  1.2797353286337168  -2.8280815800944605 
H  4.4646517331243585  1.8592222312313722  -3.397559484441992 
C  4.036148063622192  -0.06669891212098882  -3.0275996287663753 
H  4.564040684089085  -0.4146908148633597  -3.73705603114946 
C  3.322716443346617  -0.9483260638366959  -2.1870635923297295 
H  3.4027828270962353  -1.8885461018487109  -2.3035439096109864 
C  2.510435092600507  -0.44215619619878316  -1.1992107073688323 
C  2.4127887137469326  0.9545191282488777  -1.004394094176305 
C  2.816065376404336  -2.2616622752689866  1.0139528762482177 
H  3.3618071835454453  -1.5995157576596928  1.4868445225650635 
H  3.3909408026350176  -2.8084856495629293  0.4391171330642235 
H  2.362358038646625  -2.837583365845305  1.66482605060264 
C  0.7575736844527873  -2.6971391986292246  -0.9692382903787754 
H  1.4125609120259859  -3.1377475469779177  -1.5496827365417365 
H  0.05473251900342935  -2.29263445951396  -1.5178126173237323 
H  0.36027723810204537  -3.3580508350231226  -0.3632908152654452 
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
