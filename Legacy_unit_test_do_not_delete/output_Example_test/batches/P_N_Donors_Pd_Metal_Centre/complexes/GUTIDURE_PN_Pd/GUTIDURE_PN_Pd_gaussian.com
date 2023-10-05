%chk=GUTIDURE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5244288619020563  1.47248655172127  -2.4144316457596444e-15 
N  1.524428861902057  -1.47248655172127  1.9984014443252818e-15 
H  1.47933810776954  -1.7692414385412607  -0.9127756356521711 
H  1.4603937958912871  -2.1240450320880284  0.5580524858677423 
C  2.8682245567268003  -0.9243858413571799  0.1699610598604893 
C  3.9640096080959877  -1.765678489192847  0.33592750795630677 
H  3.84173537461054  -2.7062175884335695  0.381819060828292 
C  5.231097460808718  -1.2218867621553495  0.43515180494759753 
H  5.984106187003956  -1.7938538942255806  0.5313578240171217 
C  5.411790435786888  0.15568961885051857  0.39453329576378393 
H  6.286787875506307  0.522630246843621  0.45977382007590606 
C  4.322763619296555  0.9863313612625224  0.2588008430205333 
H  4.447506688444356  1.9274907989361703  0.25160851534364587 
C  3.034995815632181  0.4586189321801222  0.13183621787850724 
C  1.8739626511370102  2.699049230944899  1.3278662372840313 
C  2.0810073115212377  4.058855925042989  1.1026206528824594 
H  2.0199942268236732  4.406920659081399  0.22004544678698523 
C  2.3771267442757877  4.907803681841674  2.164915173141158 
H  2.5128698416983886  5.833946741284802  2.0024992719878556 
C  2.476766745911342  4.421353986146533  3.4439068946281486 
H  2.679390493621316  5.008034478399701  4.1633847308910505 
C  2.2815810464836996  3.0764731906750273  3.680434437867552 
H  2.352680369581172  2.7370016013313587  4.56513935267451 
C  1.9799393846410422  2.2184296937595027  2.632968015855951 
H  1.8453592051481955  1.2948135843898076  2.8059797721052453 
C  1.7649148855793517  2.4633953062320484  -1.5243547889680493 
C  2.8072783676589097  2.228677760189954  -2.4198081126019826 
H  3.4519868026137575  1.5573790791632844  -2.2250144660168214 
C  2.911931352535766  2.9634076792750057  -3.588918915202674 
H  3.6332593788850884  2.8026784354504417  -4.184750064591989 
C  1.980865277192464  3.92008947942679  -3.8889837691978952 
H  2.0478283805421076  4.408541329230946  -4.7011198312051645 
C  0.9484264483807944  4.174308704425771  -3.016201699948681 
H  0.3087178951115406  4.845076772796531  -3.22234918435809 
C  0.8398931432109278  3.4488260922068386  -1.8323450358896392 
H  0.1281380669121004  3.6316263084026006  -1.231688867752118 
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
