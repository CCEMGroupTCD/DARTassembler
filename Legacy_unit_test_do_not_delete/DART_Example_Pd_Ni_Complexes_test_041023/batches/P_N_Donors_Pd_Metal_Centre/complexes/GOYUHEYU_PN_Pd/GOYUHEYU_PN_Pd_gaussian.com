%chk=GOYUHEYU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.7146573732090034  1.2458130246951182  0.0 
N  3.88100742892924  -0.8990548784590728  0.024625476646135847 
N  1.7146573732090034  -1.2458130246951182  0.0 
C  6.088776830618955  -1.1697394782468524  -3.311017017147342 
H  5.55381835911611  -1.8892210515105097  -3.7306417417542272 
H  6.573826163378491  -0.6934225685362482  -4.031330696257533 
C  5.1211533045181055  -0.16546518892788797  -2.5975421447453027 
H  4.997166304402066  0.6265696392851794  -3.1774177906382235 
H  4.238668928282841  -0.5971021475570514  -2.4826355668395994 
C  5.633003921572009  0.29373747123052807  -1.2301605291291562 
C  5.175182683968291  -0.25414740914069733  -0.03238395469649297 
C  5.971021698938877  -0.33690393275725794  1.0979139352642142 
H  5.6383628986579035  -0.744205682652013  1.8887018885919677 
C  7.276720123005441  0.18789162407198656  1.0652531189492753 
C  7.576299047790012  1.0291743400369209  -0.006953353538009521 
H  8.339393037158612  1.5929529617702152  0.04273355687184721 
C  6.788470183838139  1.0604378837854151  -1.138111769480104 
H  7.039611179105744  1.6148409292737436  -1.8673257323533161 
C  8.349473663345908  -0.34637877751005586  1.9699008734263779 
H  7.952861214750488  -0.5578919418265759  2.851753971702939 
H  9.035258739112104  0.35466896214873955  2.1065787669412317 
C  9.047411971135194  -1.6422668786920862  1.396845720988929 
H  10.027609559251228  -1.5031587167410494  1.3910201283097363 
H  8.853813043846573  -2.404578771548091  1.9996448445475279 
C  8.592894159501444  -1.9972329660856767  -0.005302204769357721 
C  7.466700115895959  -2.8164209740532105  -0.22664678624758 
H  7.201369992825267  -3.445065533449803  0.43435276271776446 
C  6.738962442351673  -2.707143400050226  -1.4153344114563369 
H  5.986075092054515  -3.2675639335976765  -1.5638111885661985 
C  7.104275187025729  -1.79839603382562  -2.3634633441224375 
C  8.37061383706426  -1.2067581031013208  -2.269865041002267 
H  8.71745727990038  -0.7138579866089002  -3.0033088093467435 
C  9.118035393487355  -1.3371613672988178  -1.109595367313953 
H  9.994186351364041  -0.972345226360473  -1.0694091467568536 
C  2.668463007149359  -0.3058281659828058  -0.09423798626479392 
C  2.306787840921594  -2.4603711766861007  0.18010256699647692 
H  1.8672013770647164  -3.2970014127964298  0.2730688986428391 
C  3.650414720357279  -2.2397635698044347  0.20136997501726253 
H  4.319095988430853  -2.9052893859299385  0.3192692532084912 
C  2.2273640412229154  2.135390873081488  1.5467093729057175 
C  1.182992439345306  3.2197664360441616  1.8457203520369037 
H  1.2064125961952326  3.899835879215492  1.1410197782342713 
H  1.3835219333065951  3.6375335565375195  2.7095264674854853 
H  0.29112287031601314  2.8151550296709154  1.877569554102305 
C  2.193614393132256  1.0667997493646746  2.6549518470815854 
H  2.4113610931048726  1.4821239149846515  3.5155941375270063 
H  2.850515509634677  0.36719291613759664  2.454257995556751 
H  1.2984194356611392  0.6712832711832197  2.6997561281856095 
C  3.633906930619511  2.7365062372295466  1.4889884124745723 
H  3.8476203784903253  3.1513196480352885  2.350122644992882 
H  3.670520251738923  3.4141302371926177  0.7821175336072398 
H  4.284081324509955  2.028967820527156  1.2971133582869896 
C  1.7843856006410874  2.1936517746463613  -1.59318610601213 
C  3.1408912647534235  2.872347857637564  -1.8546058836947328 
H  3.114707221015779  3.332655254199211  -2.7185747410757184 
H  3.847931848246656  2.193780112524382  -1.8681010156160205 
H  3.3243844418717874  3.521025816513272  -1.1438115606871073 
C  0.6701119238710829  3.254014019071672  -1.534259508835628 
H  -0.18559041657909559  2.818153553466982  -1.342714790633894 
H  0.6149793601729694  3.7178395800071993  -2.3963374110248385 
H  0.8731831616190284  3.9015365557680735  -0.8272194556068246 
C  1.4934132829213977  1.1851199242850428  -2.681824869389444 
H  1.4311027454980296  1.6450875266788574  -3.5451448863438215 
H  0.6449022030225973  0.7343941340827631  -2.4905207676677765 
H  2.215363047899371  0.5229971669781273  -2.7158267973429373 
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

