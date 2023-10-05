%chk=ORUQUPOX_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.4935988007159084  -1.5037495211969318  1.9603690905264902e-16 
N  1.4935988007159084  1.503749521196931  -2.0323816202106352e-16 
C  2.8454637664343694  1.42143537029933  -0.10814173581498907 
C  3.5493540132273855  2.4225297638667893  -0.7728470883114856 
H  4.487684179697663  2.337355462432564  -0.895627449392145 
C  2.8907236116365587  3.540117700760909  -1.256509912403922 
H  3.374709000510184  4.229610788130092  -1.6948553536109614 
C  1.522583686244681  3.6430778931450543  -1.094365059469967 
H  1.0520217910418947  4.414694567760138  -1.385359664641829 
C  0.8637867099350385  2.5915418575815083  -0.4974028078019149 
H  -0.08283761803621936  2.633498908347371  -0.43116107183421726 
C  3.6043270608994664  0.28278436435446147  0.4997584590972592 
C  3.1671824776048387  -1.0609711557747272  0.552333789990288 
C  4.022004163038623  -2.0483841883774803  1.0511918250648316 
H  3.7351666802291543  -2.954229524185645  1.0608740794574045 
C  5.282247251545376  -1.7262219012604116  1.533293431192113 
H  5.858305110331774  -2.4082524525876403  1.8580444902922988 
C  5.695974772233483  -0.40683110709358755  1.5389344322658183 
H  6.541440881697618  -0.17646752252248743  1.9064489716396322 
C  4.876365677375476  0.580312878041953  1.0093568418819412 
H  5.184804827539789  1.479585742929473  0.9921552395463553 
C  1.6077966036914795  -1.6028659916589578  -1.808596573887826 
C  2.8416438153224792  -1.6879918142747552  -2.4471979539947437 
H  3.6395459194016677  -1.7349175748521108  -1.9343391183486507 
C  2.9083291424759192  -1.7044833407227649  -3.834768491552336 
H  3.753084984883883  -1.7487039233814188  -4.2679505742692 
C  1.7530254449315321  -1.658460213849273  -4.585793962647102 
H  1.8031690200814487  -1.675937466296929  -5.534519369452516 
C  0.5195615394460267  -1.5858679114337  -3.9608272423375306 
H  -0.2754702828059885  -1.5629225275691478  -4.481208910890496 
C  0.44282606801568813  -1.5461087267043097  -2.5755741771346323 
H  -0.40335487726838126  -1.4800930448366378  -2.1489581230699475 
C  1.1086331595922625  -3.160222896508344  0.6290749415122877 
C  1.0230185238162233  -4.26966568367683  -0.20988057348285613 
H  1.2690697888680367  -4.19496639804674  -1.1250212685089125 
C  0.5772714957347715  -5.483301750903878  0.295649794825745 
H  0.5017799001361272  -6.2343921745942525  -0.2801041668927911 
C  0.24269459624547118  -5.608210768231772  1.6309001480248644 
H  -0.06924250913815944  -6.4392565889514985  1.9670693110204964 
C  0.36344718425350897  -4.516978141474977  2.47795531552627 
H  0.1626035772159513  -4.610030679831645  3.4015963500291573 
C  0.7762811668349114  -3.293346809219988  1.9839924064431282 
H  0.834084610913868  -2.5430461613421813  2.563169166721696 
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


