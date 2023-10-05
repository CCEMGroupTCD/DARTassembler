%chk=FEKISAPE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
C  1.3486888418153813  2.7323717085757018  -0.12353949141234731 
H  0.527011329159025  3.080810304905272  0.201171379644906 
C  2.2452554257675086  3.5928376651443745  -0.6955863060410108 
H  2.065878017148874  4.525157359979606  -0.7196797055438487 
C  3.428120536874201  3.0847414658772685  -1.2439975761813602 
H  4.045907221273119  3.6578290193568144  -1.681387972310756 
C  3.6740343595265337  1.7324357826847325  -1.1353381953912687 
H  4.470967995738304  1.3627089552913512  -1.4972063241066016 
C  2.7522078836089974  0.9138896719817049  -0.49482566393660654 
C  1.063103305822053  -2.200077354121198  -1.5496832588958 
C  1.559299996149044  -1.717339521885326  -2.7714123900197363 
H  2.169611217681496  -0.9891580502051855  -2.8012445400147716 
C  1.1203397150760388  -2.3517499108559887  -3.955487826808461 
H  1.4638023364492243  -2.0668313814933033  -4.794383916527132 
C  0.19368210422281384  -3.3870932211618747  -3.913934462893295 
H  -0.0978255540721662  -3.799207400888589  -4.718776384032829 
C  -0.30176555794866844  -3.8132900397135776  -2.6984824910304552 
H  -0.9369806228713997  -4.519880847881061  -2.668936406544653 
C  0.12276864499935702  -3.2131307300137824  -1.5098582259626345 
H  -0.2318633087220081  -3.501002284821272  -0.6766897584048969 
C  2.2136018828553254  -2.6963860262396766  1.0944316527947964 
C  2.4580150569784522  -2.362324163109622  2.418498540046809 
H  2.233069566710646  -1.497366976550104  2.737540248648086 
C  3.0289458442299244  -3.2863594068083435  3.2775852147025604 
H  3.1999576141268586  -3.055373144489816  4.183448808836252 
C  3.346498094276578  -4.534976924335289  2.8112158831474425 
H  3.7457835143460745  -5.165761144344722  3.400341736515803 
C  3.09791808935572  -4.896074709818602  1.5040213526843393 
H  3.305229581994504  -5.77352833000817  1.2074500974418512 
C  2.5413495221960964  -3.9717384729814844  0.6182420566293213 
H  2.3884612164448473  -4.203053805627422  -0.29079829057614376 
C  6.764738307530508  -1.8734408510458027  0.43452546050657737 
H  7.632198449869141  -2.180373740047693  0.6714073574570925 
C  6.058291238570587  -2.4782259485377347  -0.5674980594816508 
H  6.426831955947241  -3.2184602768026886  -1.0349805285578857 
C  4.790092362135278  -1.9968994346421964  -0.8962164447742981 
H  4.2825164909381375  -2.3791545306602426  -1.6035060082111758 
C  4.303306685049922  -0.9385339779566664  -0.14580540762238015 
C  6.194108609166466  -0.8121443981645263  1.0892260000256757 
H  6.706273515746094  -0.38177474919952137  1.7632382356105105 
N  2.9505804296307945  -0.4676786815393178  -0.367638921176054 
N  1.5820488646056412  1.410397599969597  -1.2368878439247864e-15 
N  4.953729805119481  -0.3391346383532634  0.8408704492637613 
P  1.5820488646056419  -1.410397599969597  0.0 
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

