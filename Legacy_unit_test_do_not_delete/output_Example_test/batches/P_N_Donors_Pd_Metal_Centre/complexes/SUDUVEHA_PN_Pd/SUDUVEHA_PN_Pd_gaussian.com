%chk=SUDUVEHA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.7146573732090034  -1.2458130246951182  1.525680933028959e-16 
N  3.88100742892924  0.8990548784590728  -0.024625476646135958 
N  1.7146573732090034  1.2458130246951182  -1.525680933028959e-16 
C  6.088776830618955  1.1697394782468529  3.311017017147342 
H  5.55381835911611  1.8892210515105101  3.730641741754227 
H  6.573826163378491  0.6934225685362486  4.031330696257533 
C  5.1211533045181055  0.16546518892788828  2.5975421447453027 
H  4.997166304402066  -0.626569639285179  3.1774177906382235 
H  4.238668928282841  0.5971021475570517  2.4826355668395994 
C  5.633003921572009  -0.2937374712305279  1.2301605291291562 
C  5.175182683968291  0.25414740914069733  0.03238395469649294 
C  5.971021698938877  0.33690393275725783  -1.0979139352642142 
H  5.6383628986579035  0.7442056826520128  -1.8887018885919677 
C  7.276720123005441  -0.1878916240719867  -1.0652531189492753 
C  7.576299047790012  -1.0291743400369209  0.006953353538009647 
H  8.339393037158612  -1.5929529617702152  -0.04273355687184702 
C  6.788470183838139  -1.060437883785415  1.1381117694801042 
H  7.039611179105744  -1.6148409292737433  1.8673257323533163 
C  8.349473663345908  0.34637877751005564  -1.9699008734263779 
H  7.952861214750488  0.5578919418265755  -2.851753971702939 
H  9.035258739112104  -0.35466896214873983  -2.1065787669412317 
C  9.047411971135194  1.642266878692086  -1.3968457209889291 
H  10.027609559251228  1.5031587167410492  -1.3910201283097365 
H  8.853813043846573  2.4045787715480906  -1.999644844547528 
C  8.592894159501444  1.9972329660856767  0.005302204769357477 
C  7.466700115895959  2.8164209740532105  0.22664678624757967 
H  7.201369992825267  3.445065533449803  -0.4343527627177649 
C  6.738962442351673  2.707143400050226  1.4153344114563366 
H  5.986075092054515  3.2675639335976765  1.563811188566198 
C  7.104275187025729  1.7983960338256202  2.3634633441224375 
C  8.37061383706426  1.206758103101321  2.269865041002267 
H  8.71745727990038  0.7138579866089005  3.0033088093467435 
C  9.118035393487355  1.337161367298818  1.1095953673139527 
H  9.994186351364041  0.9723452263604732  1.0694091467568534 
C  2.668463007149359  0.3058281659828058  0.09423798626479388 
C  2.306787840921594  2.4603711766861007  -0.18010256699647723 
H  1.8672013770647164  3.2970014127964298  -0.2730688986428395 
C  3.650414720357279  2.2397635698044347  -0.2013699750172628 
H  4.319095988430853  2.9052893859299385  -0.31926925320849153 
C  2.2273640412229154  -2.135390873081488  -1.5467093729057173 
C  1.182992439345306  -3.219766436044162  -1.8457203520369032 
H  1.2064125961952326  -3.899835879215492  -1.1410197782342708 
H  1.3835219333065951  -3.63753355653752  -2.709526467485485 
H  0.29112287031601314  -2.815155029670916  -1.8775695541023045 
C  2.193614393132256  -1.0667997493646748  -2.6549518470815854 
H  2.4113610931048726  -1.482123914984652  -3.5155941375270063 
H  2.850515509634677  -0.3671929161375969  -2.454257995556751 
H  1.2984194356611392  -0.6712832711832201  -2.6997561281856095 
C  3.633906930619511  -2.7365062372295466  -1.4889884124745718 
H  3.8476203784903253  -3.151319648035289  -2.3501226449928816 
H  3.670520251738923  -3.4141302371926177  -0.7821175336072393 
H  4.284081324509955  -2.028967820527156  -1.2971133582869894 
C  1.7843856006410874  -2.1936517746463613  1.5931861060121302 
C  3.1408912647534235  -2.8723478576375636  1.8546058836947332 
H  3.114707221015779  -3.3326552541992105  2.718574741075719 
H  3.847931848246656  -2.1937801125243817  1.8681010156160207 
H  3.3243844418717874  -3.521025816513272  1.1438115606871078 
C  0.6701119238710829  -3.254014019071672  1.5342595088356283 
H  -0.18559041657909559  -2.818153553466982  1.3427147906338945 
H  0.6149793601729694  -3.717839580007199  2.396337411024839 
H  0.8731831616190284  -3.9015365557680735  0.8272194556068251 
C  1.4934132829213977  -1.1851199242850425  2.681824869389444 
H  1.4311027454980296  -1.645087526678857  3.5451448863438215 
H  0.6449022030225973  -0.7343941340827628  2.4905207676677765 
H  2.215363047899371  -0.5229971669781269  2.7158267973429373 
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
