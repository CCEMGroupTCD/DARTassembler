%chk=HOJOMIKI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.4331252195813178  -1.5614903473925157  -3.801734943177117e-15 
N  1.4331252195813178  1.5614903473925157  -3.577608692771468e-16 
C  1.3606648862310702  -2.7401838071845823  -1.3869618192056399 
C  0.6024502905059123  -2.443479587543615  -2.5224244750761695 
H  0.13048996269519497  -1.620191398701006  -2.5740473099950165 
C  0.5353775131572183  -3.3496497904928835  -3.5777722735888844 
H  0.01374937958471567  -3.1467656833134723  -4.344465269069968 
C  1.2319052855262873  -4.552412995616876  -3.510575586777426 
H  1.1879426670442987  -5.168301967717103  -4.23194849465625 
C  1.9923057468614669  -4.851745047704489  -2.3853821318576456 
H  2.4638301614881035  -5.675176487512671  -2.337808921471617 
C  2.0658494673020344  -3.9513127924889098  -1.3294487982843832 
H  2.5938564903446686  -4.156413781768029  -0.5676705970903585 
C  1.2481302055460906  -2.6027525349654805  1.4903192951555468 
C  0.3754223054511421  -3.6901349692119236  1.4938538112679898 
H  -0.10535281282684283  -3.910710716333991  0.7041766269538206 
C  0.2020355730245713  -4.454676682193712  2.644161170250495 
H  -0.4033878450423436  -5.1856983864944315  2.6419514875893877 
C  0.9132227772307012  -4.148303651583009  3.791283786489901 
H  0.7942981292614327  -4.6722670794493615  4.575267354151853 
C  1.7915291381552552  -3.087603377152394  3.8008003080335904 
H  2.2857433230282553  -2.889496978021231  4.587827984861811 
C  1.9587055713838257  -2.300950003309139  2.656695903010152 
H  2.5550065477944086  -1.5619314260076667  2.6722206553457637 
C  3.1567053566349097  -1.0473917707727476  0.12640025582212472 
C  4.367342644048727  -1.9539771828208639  0.19074362366110478 
H  4.610932677803528  -2.286928251566287  -0.7089344861518522 
H  4.20168638189182  -2.728061330710994  0.7863148996587154 
C  5.467044898388241  -1.0484685561674243  0.7589939108205964 
H  6.3456385247261595  -1.2669589653907094  0.3577439751995775 
H  5.5306859897965595  -1.1452584948643734  1.7425328115517371 
C  5.029223950651794  0.3801659555509103  0.375318501587115 
H  5.254536401220833  1.0265745111244364  1.091440368939252 
H  5.457612654500441  0.6701735279679863  -0.4681162714425816 
C  3.5261678858454877  0.24849678427276214  0.20904226941921958 
C  2.7235566759016683  1.463010932994095  0.10751149463534833 
H  3.204076310781409  2.282782276529802  0.12211337166205985 
C  0.9868565084789387  2.9494818541086247  -0.1074535935879726 
C  0.7151813982343873  3.654957130878356  1.0693238358099155 
C  0.35584458151652254  5.003832488012115  0.9418581813443092 
H  0.1588703472761217  5.512990470530305  1.7191317946127573 
C  0.2842712084634522  5.607795014309394  -0.3146974318692323 
H  0.06413676758323228  6.529352424128126  -0.38464082574273345 
C  0.5311695727718883  4.87344095511536  -1.4556376070962609 
H  0.4676308135350076  5.293958937244733  -2.304492160468315 
C  0.8754158953470657  3.5125936731911476  -1.3805933223508482 
C  0.7977394020284875  2.9911985929336695  2.4059134098060118 
H  1.0616013178959727  2.0544857099866407  2.2896863721779086 
H  -0.07716923635262285  3.031283387632862  2.8447743487377575 
H  1.4625393990678366  3.4530115437496436  2.9600486752377826 
C  1.1162810260290326  2.6957889570237974  -2.6241180444486143 
H  1.3419867585873304  1.7756261859895035  -2.3724932104341416 
H  1.85733480592329  3.086470568356647  -3.133525920193398 
H  0.3061512226581131  2.694445954805755  -3.175292933240947 
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
