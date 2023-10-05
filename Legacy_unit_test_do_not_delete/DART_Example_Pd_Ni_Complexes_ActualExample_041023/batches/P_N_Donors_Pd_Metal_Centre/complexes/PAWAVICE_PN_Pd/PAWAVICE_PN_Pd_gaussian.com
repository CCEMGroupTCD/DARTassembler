%chk=PAWAVICE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.512254027933138  1.4849874595430097  -3.3712151236062787e-16 
N  1.5122540279331376  -1.4849874595430097  4.440892098500626e-16 
C  2.6135645077165597  -1.4394621727363424  -0.6436789769601078 
C  1.2465933854487035  -2.706236493227173  0.7234526596538797 
C  0.9738581722011054  -2.606697329411278  2.0832810084985676 
C  0.6775765904960414  -3.7538292927507055  2.7846040959749367 
C  0.6110982630619858  -4.9677555070197155  2.147214648076638 
C  0.8745807410771705  -5.04965102015585  0.8014734735905585 
C  1.2097169662681408  -3.9097433838173417  0.08857214017764561 
C  2.402699483157704  1.3536174882032959  1.5852790170503757 
C  1.661805328600425  1.0825169159978647  2.722707002240104 
C  2.289892595528242  0.8650179041884352  3.92685492711693 
C  3.645067039584332  0.9409620909728176  4.002354640169198 
C  4.38993245709989  1.237793557236949  2.888942694649978 
C  3.766524717766543  1.4541123839958332  1.6778012595859266 
C  1.1925232259419822  3.2461006381046738  -0.2795003516178913 
C  0.4125162714970776  3.6192553616944125  -1.3731965769383896 
C  0.07398893154428987  4.936261537357243  -1.5683230573823697 
C  0.5124622377839438  5.895777350668212  -0.7189095928543613 
C  1.314100708934041  5.569822618263297  0.3133754908289224 
C  1.6482710838927221  4.228942220786557  0.5530300342087773 
C  2.7556880370350147  1.029746581476576  -1.2587395353869661 
C  3.137593110277913  -0.3074613338085077  -1.4216923086995363 
C  4.115708049677714  -0.6296404518377162  -2.3526515834419395 
C  4.674650393667977  0.345163891254499  -3.1533687112613085 
C  4.312604807733001  1.647362699448538  -2.9945840531767334 
C  3.3640702262947193  1.9937059137934434  -2.047366276900856 
H  3.1413824463714635  -2.2056478214354502  -0.61779219729112 
H  0.9902758661681847  -1.7810965863738306  2.511887719156271 
H  0.5210792466932003  -3.705415654953302  3.7007891215335773 
H  0.3886783110214693  -5.732981824385185  2.6258395828069494 
H  0.828064546966653  -5.871064727834824  0.3681583636271964 
H  1.407783870247505  -3.9671565386908765  -0.8180029759575365 
H  0.7335835435609087  1.0467383983689342  2.67220100799841 
H  1.7902907706873012  0.6684470623782084  4.6865348733875045 
H  4.070861420766434  0.790182645960599  4.815643539586633 
H  5.316152602998827  1.291853343724779  2.950781766569868 
H  4.270100198028503  1.6676532782048565  0.9259038183219159 
H  0.12013857329093192  2.9718476049340565  -1.973871837023402 
H  -0.46229441680610117  5.170817602300167  -2.2912950918903503 
H  0.25909975401024443  6.7807618881163805  -0.8488046285319297 
H  1.6480071655244473  6.236613663289843  0.8683654101244501 
H  2.181771566627102  4.008301307023668  1.282658961880966 
H  4.397324575330268  -1.511105579093957  -2.437546824242171 
H  5.299692182751323  0.11243491244793336  -3.8026901989106325 
H  4.70247173998197  2.3053926401867  -3.5224397268517733 
H  3.1331527782023874  2.887563186102319  -1.9375389628707822 
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


