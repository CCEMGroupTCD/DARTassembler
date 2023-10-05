%chk=ADATOLUP_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
C  3.2261190268552724  1.109803734447016  -0.35856077909678674 
H  3.3017775763556685  1.1294945848748152  -1.3348678590142404 
C  3.6449663620184083  -0.3180461307046971  0.05917829306786142 
H  4.568520363187362  -0.4482868574398733  -0.20895148898611776 
H  3.621868044150336  -0.3614189446996888  1.0296444267983675 
C  2.862959130179253  -1.467019184718456  -0.4580415446236024 
H  3.298778675289708  -2.2870906700954228  -0.14165628426633498 
C  4.2097662668734985  2.1197134596349434  0.13493485834634175 
H  4.094407692210426  2.2397214466695443  1.082115613165415 
H  4.064980252722668  2.9540791259825365  -0.31261827249586177 
H  5.100346869468836  1.810412995766607  -0.041162453761914075 
C  1.2734009134504176  1.9747824753902272  1.6837656241250507 
C  0.11976463983006358  2.6887849456735364  2.0108410139309676 
H  -0.4884846760890187  2.9090502344618  1.3461542740007841 
C  -0.13346326541647802  3.0593830255664365  3.3126185299390487 
H  -0.880999549141043  3.571120031060047  3.5233074516310134 
C  0.7234307777249391  2.680222874781767  4.253475728926086 
H  0.5708067456989809  2.9696707527208477  5.125185767186948 
C  1.7990627050134158  1.895477786068116  4.024416582115675 
H  2.321681397459378  1.575477390006272  4.722385219538854 
C  2.097576391825451  1.5877155042768987  2.6922734529089953 
H  2.8719525222980846  1.1118933663943509  2.494430733419201 
C  1.283875175446191  2.995023980929644  -1.0056569496276295 
C  0.9727361421822384  2.887660585297002  -2.267732202240491 
H  0.7962835878099247  2.0459888304833025  -2.6191302859067442 
C  0.8968026535852256  4.0198966167131545  -3.1295265132408034 
H  0.691679659731528  3.947687190353919  -4.033801419526196 
C  1.1435761878091593  5.172701027853795  -2.546920211814784 
H  1.0506416694810756  5.94300752558891  -3.0719584721827253 
C  1.499063197865417  5.3409627719564625  -1.2777711972570633 
H  1.7070868556492798  6.1886648114396925  -0.9566326756527572 
C  1.559375137540149  4.207024035788335  -0.4128419212176153 
H  1.7747150531472036  4.286045823700366  0.4867323212891111 
C  2.8380831460038616  -1.5212625774188488  -1.9387412743999282 
H  2.3978534605426427  -2.3257601455321812  -2.222716179873016 
H  3.738482917376117  -1.5132991794563186  -2.272844396213711 
H  2.3630544435779717  -0.7590627357107756  -2.2773743400798265 
C  1.4099558835580897  -1.9170429410922525  1.3901006317274816 
H  0.4746629782607372  -1.90320047420267  1.6497608662857493 
H  1.8623286159650907  -1.2440383879764414  1.9221610535406686 
C  1.9586563913050945  -3.2408063387040267  1.7853948638401267 
C  2.95986054801282  -3.362496778796007  2.671145355446776 
H  3.3025947473340223  -2.610837887851561  3.095925922657452 
C  3.4867711096569955  -4.653556781346933  2.968949351123501 
H  4.155689598300909  -4.765920273370507  3.604924878334002 
C  2.986489999087593  -5.690663870339577  2.29168086937592 
H  3.3462492113640545  -6.5336640211243715  2.445148536847095 
C  2.0282558680734253  -5.5833865924066055  1.4405645754425709 
H  1.7324751288359983  -6.340563371679762  0.9909385761657361 
C  1.4349259475055711  -4.36932767705415  1.173911363550853 
H  0.7042467045516012  -4.3076450073973405  0.6002736536849138 
N  1.4867808513698308  -1.5104908804756156  8.659759136643608e-16 
H  1.1481167231237248  -2.3060937043877994  -0.4616303335161925 
P  1.4867808513698315  1.5104908804756148  1.4808512519597407e-16 
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


