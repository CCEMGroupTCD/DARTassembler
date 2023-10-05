%chk=RITUSEKU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
N  1.3647190333544847  -1.4891749259237481  -2.586832139131677e-15 
C  1.3561579862224455  -2.0946311583006882  1.362563715628725 
H  2.0992210057482414  -2.7275537677216994  1.4429374022239876 
H  0.5078702407053672  -2.5657566438206625  1.5041557270528638 
H  1.4543053541249922  -1.388427352627793  2.034221172244034 
C  1.212596400162362  -2.5809961563889554  -1.007384536895503 
H  0.39385449950127194  -3.085674127335723  -0.8200977361247529 
H  1.9850592776735203  -3.181723004811466  -0.9583664131905547 
H  1.157634155174473  -2.191615117327806  -1.9050821296546452 
C  2.6839636773907642  -0.8342063772188504  -0.2504312867615225 
H  2.7875820109570553  -0.6727662295136568  -1.220599081350629 
H  3.40847100514552  -1.4455442020450904  0.03519286481955611 
C  2.8176387253937616  0.4791706115769892  0.4924993740965308 
H  3.658547303762955  0.9366694910839658  0.24145324671154716 
H  2.819218746394758  0.32777375679839554  1.470746058679114 
P  1.3647190333544847  1.4891749259237481  -1.6849354283248703e-16 
C  1.8634363447623974  2.245868394490583  -1.5759621340407926 
C  1.1566204754630194  1.9417055778093726  -2.7406007551826623 
H  0.3779470641057482  1.4002269469495292  -2.6920643691675847 
C  1.5838175768555804  2.4266534421793233  -3.970439954185855 
H  1.1003767639056763  2.210640874145747  -4.759543959261164 
C  2.7003191451698942  3.213355406875167  -4.047930319363872 
H  2.9918553060776087  3.5391957440665855  -4.890919193277533 
C  3.4129997131224927  3.539969058188353  -2.888782174889375 
H  4.179959485619552  4.098355709430276  -2.94171850538272 
C  2.9944327332042873  3.0470482285251324  -1.6702947523008826 
H  3.486916156918781  3.2578890737033683  -0.8851079494501177 
C  1.2626973266233588  2.8145479131003146  1.2389494768626985 
C  1.3988930518110902  4.177726532791864  0.9499990914045241 
H  1.4981557198147522  4.468009354826254  0.05024247099399907 
C  1.3867242153557622  5.109549810127286  1.9912754856607364 
H  1.4916258329324514  6.035055428041893  1.8019658666808926 
C  1.2229561081740368  4.689383455518627  3.2854370835785165 
H  1.238652970186369  5.3281370924763225  3.988375985472735 
C  1.0331786747005398  3.346220415918524  3.588917337955012 
H  0.8889777095032224  3.0717636325782647  4.485950735167553 
C  1.0601629826394119  2.411882739545879  2.5657646740599063 
H  0.9399337188982667  1.491405718251881  2.7660212779709705 
Br  -1.5839191898578664  1.5839191898578664  0.0 
C  -1.2937473829822488  -1.2398794889499034  -0.4598005675376139 
C  -1.601562203225019  -1.4938654694550033  -1.7919791134336869 
H  -1.1309417375171225  -1.0428217585005106  -2.480879663633267 
C  -2.60218398542082  -2.411889969420291  -2.121542410022959 
H  -2.8231586547736818  -2.5626138747145775  -3.03270050271447 
C  -3.2683208079706407  -3.09759503004543  -1.143064365570095 
H  -3.9596020908062854  -3.707191723557817  -1.3683280760515548 
C  -2.9347449854640226  -2.896087320001541  0.14198826503765638 
H  -3.3668350933519986  -3.405606940434718  0.8171022443470677 
C  -1.9675509241539992  -1.9516175659140573  0.5112421903797603 
H  -1.7755298364628027  -1.8022912010877177  1.4301324877763764 

-Ni 0
lanl2dz
F 1 1.0
3.130 1.0
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

-Ni 0
lanl2dz


