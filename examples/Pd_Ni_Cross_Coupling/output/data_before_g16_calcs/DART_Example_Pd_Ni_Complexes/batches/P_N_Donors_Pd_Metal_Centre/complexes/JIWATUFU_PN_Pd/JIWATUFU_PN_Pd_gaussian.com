%chk=JIWATUFU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5306830844430201  1.465984070513728  2.70831301769464e-15 
N  1.5306830844430208  -1.4659840705137284  -5.551115123125783e-17 
C  2.7048649880898403  -0.9842149816091661  -0.1876607367683345 
C  4.056436154412422  -1.533119305459024  0.16306057085979236 
C  4.920173168302533  -1.5557216895274868  -1.1275151933131915 
C  5.224548172947795  -0.054596381099797364  -1.3892029938409716 
C  4.415090024591143  0.6588551331676662  -0.3049022501132351 
C  2.91912049676923  0.43717522781499984  -0.6505295513947538 
C  4.611017575325333  -0.26953119269674897  0.9031264725632191 
C  1.3149064657813498  -2.8417070306501575  0.3895034374590285 
C  1.3188947973685692  -3.8382433535094176  -0.6051784164185392 
C  0.9510334659475994  -5.128478607796437  -0.21585468295228258 
C  0.6208617880946186  -5.420730736327316  1.0824514159754546 
C  0.669330556566359  -4.424726816933893  2.053889793349844 
C  1.005078645507899  -3.133540941812697  1.731563193982633 
C  1.7228713490946366  -3.578402307033031  -2.0496810007781594 
C  0.5866957324224179  -3.8561333841908625  -3.037028698836087 
C  2.945961782999803  -4.4174243944756135  -2.4340393586594757 
C  1.0839254130713862  -2.0598840640483758  2.802632064818776 
C  -0.04828351949222931  -2.1205673564451515  3.8221689111598276 
C  2.4424699855812984  -2.1246469914973285  3.5114356114606724 
C  2.0120037123886245  2.087647377319398  1.6241735024576556 
C  3.0413370390157146  3.045934604483947  1.7594697808548032 
C  3.396695318442969  3.520059800644112  3.025960571033952 
C  2.7287440519988637  3.051073667828175  4.132991759791556 
C  1.708008954709478  2.1081801959259816  4.027119328035877 
C  1.3670079623190774  1.647233369427929  2.765774301355954 
C  1.4091963994062473  2.9273452736059675  -1.0709040430662367 
C  0.16584332929612056  3.5794281674296413  -1.1103291953208798 
C  -0.025585499944217327  4.67719744908087  -1.8925825478855187 
C  1.0160286853916105  5.166583070268819  -2.667969614247663 
C  2.2464323792668717  4.538617576975707  -2.6271899706580992 
C  2.4479839118595095  3.4148829620267818  -1.832556297274464 
H  4.0659268051997275  -2.3332113962202623  0.6727865500377761 
H  4.4146310467380125  -1.9129667445958365  -1.8381384577460946 
H  5.737949710026472  -2.0131949141614927  -0.943807037400173 
H  4.880724004232601  0.18431740182914225  -2.243935604271818 
H  6.147327445502299  0.11137075312338518  -1.2606084974521343 
H  4.660788824755782  1.5627904600608318  -0.15117231030902084 
H  2.8275983340629516  0.46982463442364564  -1.595780142498341 
H  5.5355163604665485  -0.38591926500384277  1.0724522913183803 
H  4.004874112279287  -0.022153951501409308  1.5944149382164068 
H  0.9315984783817401  -5.81658669525096  -0.868190057769562 
H  0.35796741056116566  -6.302881409729116  1.3228390106596446 
H  0.462762634326102  -4.645154403085465  2.951780099993576 
H  1.9635544286292124  -2.6681079888059873  -2.127331496135992 
H  0.9170472785006287  -3.809405966493892  -3.919926475727134 
H  -0.0906442510265697  -3.188271267069153  -2.9180160972552227 
H  0.21762321577959898  -4.70607336311803  -2.8651973871370164 
H  3.2482547039023513  -4.152242485054341  -3.2917218901872376 
H  2.6931951946651105  -5.33948857797408  -2.4594573312565395 
H  3.624998432039018  -4.293651143012919  -1.7933324432262405 
H  1.0304799678232583  -1.2126071734033568  2.36498341011828 
H  -0.8845740196746843  -2.0622746664229106  3.3677485203381217 
H  0.031643179189961845  -1.402795825932368  4.4306832518479355 
H  -0.00508503503534552  -2.9495320221895316  4.287174689498034 
H  3.1418576183768545  -2.086058134266311  2.863790763422622 
H  2.5110512377710696  -2.9474959905163134  3.987867838185049 
H  2.528129142920986  -1.3999108286905368  4.112443767049204 
H  3.4901488321885745  3.3686877261015593  0.9879303024124659 
H  4.089937701730053  4.15869609542453  3.1193529208468886 
H  2.9678765195038608  3.375773809237756  4.989942925934404 
H  1.260694167589828  1.7913269883483411  4.80638473601432 
H  0.6650421056133315  1.0100464662491935  2.682076293990991 
H  -0.5521954378799792  3.242989914163244  -0.5817086408857464 
H  -0.8815991451528631  5.111880119116085  -1.9172699047243802 
H  0.8843559989098151  5.925842442445999  -3.221632471405145 
H  2.9668716115552876  4.8762573814972825  -3.150948011152298 
H  3.297813578724611  2.994703972780234  -1.8104010651219116 
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

