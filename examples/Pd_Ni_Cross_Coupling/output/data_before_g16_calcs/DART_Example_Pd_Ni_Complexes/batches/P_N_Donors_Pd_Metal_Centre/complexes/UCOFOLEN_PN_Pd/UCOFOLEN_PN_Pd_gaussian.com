%chk=UCOFOLEN_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.6504477687888215  1.3297075477337108  3.806601636316281e-16 
N  1.6504477687888217  -1.3297075477337112  -3.3306690738754696e-16 
N  3.943300902101465  -1.7296760753294547  0.26911225537195216 
C  3.095844683486401  0.4190897039351087  -0.7184972035703144 
H  2.942395490450762  0.35257693941727464  -1.6847239147064315 
C  4.459860504274577  1.072062270331248  -0.5277381194520716 
H  4.373844167047332  2.032154011688846  -0.6409903066580852 
H  4.773043602587403  0.9029682538437758  0.37417883931926094 
C  5.482351100143566  0.5367658878950234  -1.5257842027118889 
H  6.3689550212105415  0.8330737946009146  -1.2666604454601864 
H  5.291704348269581  0.8972589119565795  -2.4062366861559745 
C  5.45762822444212  -0.9843621213783793  -1.5872801279507436 
H  4.719124265214777  -1.2628078033864552  -2.1524137089531528 
H  6.278305710355423  -1.294140080225079  -2.0023215720594494 
C  5.318296284866627  -1.63343934626967  -0.255560947064021 
H  5.685386025998404  -2.530462031470572  -0.31229580575567395 
H  5.854052766453246  -1.139621479236327  0.38321210859547683 
C  2.8937874154071386  -0.9725619326626407  -0.13792295959185402 
C  1.365936987157364  -2.6223333532572353  0.6620640769820746 
H  1.435965156473399  -3.342924863498015  0.01756725707435184 
H  0.46391999039867393  -2.6180768609107448  1.0179991129526242 
C  2.3560751424781197  -2.834145116287405  1.780405375700905 
H  2.1764412614541064  -3.673809647712322  2.229488679170577 
H  2.2863196884259978  -2.1163783361960506  2.429104162525159 
C  3.737850772243993  -2.854357714949276  1.1711248903409168 
H  4.400368038194502  -2.8233175286539165  1.8789937617798387 
H  3.8607132076423465  -3.6847533774255212  0.6832723760411209 
C  2.0034619837045637  1.5692299303379351  1.7682156662294104 
C  1.24719326560439  0.8919408779181506  2.6964255484734303 
H  0.5278781115838487  0.36884839702979444  2.4238301419925126 
C  1.5701328703932866  0.9974140928660775  4.052964741636457 
H  1.069010186628725  0.5319794331510077  4.684413427875109 
C  2.6104119340137  1.7751344997048126  4.456973229888454 
H  2.8142229841417783  1.8416043863732305  5.3616274511310085 
C  3.3675424946377883  2.471179548540882  3.5256090131827316 
H  4.078251935825861  3.000568128157314  3.806622308202444 
C  3.0655492276681526  2.3782335280338605  2.1847347595791518 
H  3.5665947977423014  2.850785158707039  1.5588845354382608 
C  1.6899290827788516  2.9688846124310184  -0.7758829012266927 
C  1.4218010232503573  4.127112936865876  -0.059502895803371494 
H  1.2683700628360817  4.080083753564194  0.8564010818625988 
C  1.3817638468376248  5.348184357997013  -0.6991923665637482 
H  1.2078827066547624  6.119518792962101  -0.20991396061699485 
C  1.596531537914983  5.4371289347211  -2.0580716993746857 
H  1.5670699545163083  6.2645969534409165  -2.480526600044036 
C  1.8533686091097694  4.303515782288377  -2.7864507418324287 
H  1.999115105864293  4.361329923549006  -3.7033565678331577 
C  1.8951234779165829  3.057929854147367  -2.144013821545274 
H  2.0618260029888584  2.2878964872276133  -2.6371992349191484 
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


