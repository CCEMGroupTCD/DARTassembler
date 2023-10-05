%chk=UYEJAZAY_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3325013198867792  1.518071221155318  2.15087507499852e-15 
N  1.3325013198867754  -1.518071221155318  -1.6653345369377348e-16 
C  1.098875417723006  -2.7619245290145855  -0.7787377006192355 
H  0.18959601464396547  -3.1048603496434763  -0.5914807646573612 
H  1.7520056328366027  -3.4534420384507225  -0.5047118692658457 
C  1.240907987886799  -2.4906232888220736  -2.247163553729537 
C  0.14458902493804615  -2.282171909472972  -3.0565878700967137 
H  -0.7237882485707301  -2.297279293553771  -2.6721328560837305 
C  0.2754028776060513  -2.0526879727872656  -4.408791593192259 
H  -0.4961083817268599  -1.8973086422835972  -4.940805954252919 
C  1.5163986132043212  -2.048824919513957  -4.98630653762779 
H  1.606473616742746  -1.9065153776496517  -5.920028133196768 
C  2.6436546636996625  -2.2536270497877844  -4.200068199681221 
H  3.5086770514907806  -2.2484187829490367  -4.593653505017594 
C  2.499781192750221  -2.4665920974710924  -2.8289131527145255 
H  3.26947954666538  -2.5973478302943107  -2.2887347356020133 
C  2.4702368815634945  -1.430453039513766  0.5816621149275825 
H  3.013435045187692  -2.2089824052839373  0.5371438976587459 
C  3.0504204531271633  -0.29329925124405287  1.3071802298713022 
C  4.099234957912563  -0.6122654323112967  2.1731954028833957 
H  4.350826792132617  -1.5220537978063664  2.2808140767444156 
C  4.774129920368541  0.3541912526941886  2.8701957799706674 
H  5.465441927715648  0.10785471736654983  3.4721055405332306 
C  4.450438799169602  1.6778325166862524  2.7000902177027517 
H  4.930774098260656  2.349956710058078  3.1703866871041293 
C  3.410789892319156  2.0339439049178125  1.8320341137230052 
H  3.1915595098893235  2.951029213959159  1.7157079113427844 
C  2.6941930033721886  1.0645787496518593  1.1369599058086504 
C  2.180074256620882  1.7984033452806454  -1.5841017982122838 
C  1.6686138068226852  1.2026930173345267  -2.7293656121305387 
H  0.9005955399666257  0.6461704687040086  -2.6665521212131944 
C  2.2669693447839596  1.4114500773351892  -3.963503128816185 
H  1.913486067530487  0.9898621102614626  -4.7383103104561 
C  3.362887724389763  2.2228413630896116  -4.068158977405854 
H  3.7568081700235263  2.3784338489416528  -4.918407459238123 
C  3.9034531048202648  2.8202871102725124  -2.93546462381654 
H  4.66865525477675  3.379298017121907  -3.011615564227046 
C  3.3267023225050547  2.6020336051781694  -1.6962445817687717 
H  3.7072719085263506  2.9967237367413038  -0.9207678314775241 
C  0.8573111889462015  3.1807300607577833  0.5763719480929608 
C  0.9276020320860093  4.324272564434814  -0.21572206126121762 
H  1.29385262617485  4.262677404275702  -1.090812558623958 
C  0.4791453805330643  5.542301952142189  0.2436176195406105 
H  0.529549864656341  6.309972716701638  -0.3126202928862423 
C  -0.047398134233413325  5.6419029180869416  1.5165574706884521 
H  -0.35267366943879064  6.481981556841087  1.8383235010801813 
C  -0.12703885911730795  4.5273939962246414  2.319769521164107 
H  -0.48576153507433206  4.600535024427192  3.1962718895281266 
C  0.3141139719259516  3.2987682051497362  1.8554911723601855 
H  0.24580513456388786  2.5322037490273326  2.4117143625887363 
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

-Ni 0
lanl2dz
