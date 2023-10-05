%chk=OPEFIQUT_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.3586617156231353  1.62669552851786  8.540646973175044e-16 
N  1.3586617156231346  -1.62669552851786  1.1102230246251565e-16 
C  2.815605834978058  1.0905849536061605  -0.9437600831487944 
C  3.882971800779738  1.9677379333651106  -1.14948656454647 
H  3.846447798929853  2.8259911439533574  -0.7935142850023458 
C  4.999746504998214  1.5750422052456443  -1.8796323206352668 
H  5.699911625002026  2.1707505925939046  -2.0217996937311726 
C  5.062988823242186  0.294935400687512  -2.3926540727151413 
H  5.808186436593825  0.028092984642084318  -2.881071263295156 
C  4.017546930428587  -0.5997784959727261  -2.1825871741486003 
H  4.07921071353263  -1.465457008215922  -2.5169005919076155 
C  2.878460347597083  -0.210705578573942  -1.4764650588171886 
C  1.757318955401707  -1.2263357693726507  -1.381776388721911 
H  0.9780566262063479  -0.8625229109582543  -1.8312236118186664 
H  2.025941882450429  -2.023559326126137  -1.8643648610470418 
C  2.556410480695667  -2.022677046698401  0.7841642831873376 
H  3.0327639902467305  -2.7106449056988278  0.3139248027938317 
H  2.2805633272586294  -2.353000709350616  1.6422377333259197 
H  3.1276326718687164  -1.2611331417190126  0.9018954137357731 
C  0.47312452926818094  -2.8099301266689474  -0.12192435081151237 
H  0.9600024276927223  -3.534226793958313  -0.5220656505785168 
H  -0.28205296876130936  -2.5889691916330624  -0.6726057030089463 
H  0.16754711374412512  -3.07213292831347  0.7499342294246479 
C  0.3747517165754297  2.6713563210315834  -1.1758187753045077 
H  -0.49560311351235176  2.8357442905759815  -0.754961465504664 
C  0.10595177186506066  1.9245114697969288  -2.4782720980412676 
H  -0.47206910231487953  2.4497988422768597  -3.035296958640326 
H  -0.31512133924362384  1.0832936905565798  -2.2848524170592306 
H  0.9348513713681887  1.7690358011384943  -2.9350716497034584 
C  0.998706549870785  4.043788860410441  -1.466981463266889 
H  1.857086659127051  3.9232584467203995  -1.8820481897205035 
H  1.106088667572682  4.527446214357218  -0.6455112965711801 
H  0.4252608676034486  4.538382944487784  -2.0573070695933064 
C  1.9904309789501493  2.7546825505709873  1.3060267478099363 
H  2.4434302580106815  3.509552314747952  0.877255967203259 
C  2.9852722930566853  2.0335968836860796  2.2059036009181527 
H  3.391320604918838  2.665445783483874  2.803789264116593 
H  3.6627148185387255  1.6209410496766181  1.6663899798421986 
H  2.52661243072727  1.3602835710233565  2.7142241333375337 
C  0.7951904501273276  3.288178570438157  2.106622572833572 
H  0.28801192670924447  2.550311130614539  2.455086812938064 
H  0.23677133302934905  3.817422675992694  1.533206086852125 
H  1.1129763503334769  3.829345500307378  2.8337960429533906 
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


