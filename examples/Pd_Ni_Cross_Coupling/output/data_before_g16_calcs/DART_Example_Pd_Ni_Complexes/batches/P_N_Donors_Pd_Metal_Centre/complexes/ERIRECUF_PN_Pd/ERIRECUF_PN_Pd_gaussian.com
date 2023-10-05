%chk=ERIRECUF_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.3689013304471584  1.6180881148750832  1.279050760605033e-15 
N  1.3689013304471584  -1.6180881148750828  1.1102230246251565e-16 
N  2.7289452656250823  -1.495665152751483  0.16136100126940195 
C  1.1646296925356  -2.9338796749437153  -0.22112918938277457 
C  2.372082494995017  -3.640332379522322  -0.15984888237883899 
H  2.4928122175532255  -4.576400487356445  -0.26560942127149684 
C  3.3495694194014822  -2.6999662919081953  0.08586659104463515 
C  -0.1619421124400291  -3.5317613175480096  -0.5373046591468216 
H  -0.8051579001005644  -3.285322497971868  0.15989788582139508 
H  -0.07912812696986404  -4.508333350599492  -0.5755628356877956 
H  -0.4740028918945731  -3.1961348334229935  -1.4030627218179719 
C  4.823308842149342  -2.864405096293419  0.25457830803890463 
H  5.290827647307918  -2.301639174197577  -0.39726052461126404 
H  5.065235990629547  -3.8030340226403947  0.10907451615371616 
H  5.080037493884669  -2.597270846294081  1.1624468568142552 
C  3.38515608180795  -0.22703299599253302  0.46859541594611376 
H  4.365797569877025  -0.352258246791463  0.4181909673326374 
H  3.1597384066243546  0.03565855095873727  1.3963232655201092 
C  2.9795612731366248  0.8991740895472324  -0.47705487441439204 
H  2.921732687250402  0.548138215926808  -1.400318921199462 
H  3.6722592320519842  1.6068063745099057  -0.46211141454089266 
C  1.1985831315582116  2.957749129051276  -1.2494285766796214 
C  0.9766948688036801  4.291109512930347  -0.8821868775853648 
H  0.9001433314814675  4.518758549230473  0.036373887092319694 
C  0.866518747801315  5.286029953418398  -1.8505289243126297 
H  0.7031259575815938  6.184458871865913  -1.589589765500763 
C  0.9961415370714137  4.967698324711167  -3.1977038040763683 
H  0.9162577272588959  5.645021262148728  -3.8585846143219205 
C  1.2433243581917743  3.6540888367852915  -3.5708129920818394 
H  1.3533401200313147  3.437283171645677  -4.489927813913299 
C  1.331293451925082  2.649744901457976  -2.607657050336696 
H  1.4821093436937496  1.7511153140870481  -2.8746823702813966 
C  1.7310192202171126  2.5827143120834464  1.5101406899914802 
C  3.0245988629922866  2.822263803508962  1.9908016504984007 
H  3.775408209931717  2.4807826969360622  1.520319814077343 
C  3.2123633982993303  3.561483419550912  3.159195971153787 
H  4.090944369245214  3.6957311765198027  3.4973544506285057 
C  2.123365085254341  4.1024263335698885  3.8320520909578915 
H  2.259169650153116  4.607969781717118  4.6249300337637695 
C  0.8319668356943434  3.903583920698672  3.3455263815889165 
H  0.08601146653502356  4.287137488823059  3.791411330428238 
C  0.6479534692696937  3.136229220676136  2.197329976207679 
H  -0.23209530528368272  2.9853615251805508  1.8755505442003497 
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

