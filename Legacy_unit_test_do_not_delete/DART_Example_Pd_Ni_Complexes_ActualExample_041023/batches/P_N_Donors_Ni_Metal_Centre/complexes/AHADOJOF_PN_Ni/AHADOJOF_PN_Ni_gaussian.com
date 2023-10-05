%chk=AHADOJOF_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
N  1.379458209769328  1.4755321235066357  -2.5242935189902206e-15 
N  1.2396283185760932  2.646783984084804  -0.6717235820045937 
H  0.526677942041543  2.7378648781290558  -1.2097434155975135 
C  2.372827603534954  3.364242813435215  -0.5975703032106398 
H  2.5253377331409004  4.211348416957767  -0.9997091559074057 
C  3.2592869967172504  2.6638941774840497  0.15106094938321168 
H  4.140546370099555  2.926332547229416  0.386036839446898 
C  2.612989136366833  1.4754446330359487  0.5083087009612802 
C  3.1635034484590876  0.3579920891485928  1.2645229922446966 
C  4.244573575384611  0.6283054314647931  2.1201434324137693 
H  4.515568240933903  1.5290711442717397  2.2528203025729088 
C  4.9203006359069565  -0.37624908238570454  2.771415156335877 
H  5.647485948714534  -0.16508347979650484  3.3426198036559462 
C  4.540243490476239  -1.692760169506923  2.59261509007748 
H  5.020230083476162  -2.391423094953807  3.022529601517113 
C  3.4466495715964767  -1.986805787415186  1.7772158350017084 
H  3.175959349995851  -2.891480839034959  1.6699490884929777 
C  2.744644297749097  -0.9785153196588534  1.1176472349894173 
P  1.3794582097693284  -1.4755321235066352  -2.220446049250313e-16 
C  0.8832328753775509  -3.0911292304911244  0.6555779323704616 
C  0.32574998731756133  -3.1067527398900348  1.9295605050434 
H  0.17324646358301443  -2.2912094820853937  2.391599724047534 
C  -0.00987365601426915  -4.320378943936517  2.52904253720817 
H  -0.3922903108839155  -4.3323549658193885  3.3990645179524543 
C  0.21456708073088548  -5.506760176925056  1.8524737781733824 
H  -2.4149781393978387e-05  -6.334739094182773  2.2660322097933348 
C  0.7449868699802509  -5.4922620411954455  0.5860814766992072 
H  0.8861759682471856  -6.310561053981777  0.12502184476046718 
C  1.0791510078024096  -4.285335067780764  -0.029841058882709968 
H  1.437203083707868  -4.279936164650774  -0.9095876887624793 
C  2.1704279241256144  -1.7061325046026221  -1.6155891851611304 
C  3.2275775957388007  -2.6076630296302046  -1.7787049239836552 
H  3.5353328356807996  -3.1217128750496417  -1.0413857741492714 
C  3.823764369486427  -2.746855180534995  -3.025238989276282 
H  4.538144934708103  -3.3613699346056283  -3.1383952752840507 
C  3.386474865211385  -2.0011766673588944  -4.103154637385362 
H  3.7928014563155026  -2.1157434004763394  -4.955029824853874 
C  2.358986334567833  -1.0873502629594118  -3.943343148379446 
H  2.0674503068889405  -0.5644055789994153  -4.6812614650536055 
C  1.7543785847402456  -0.9355001836351371  -2.699398293388799 
H  1.052774234453937  -0.30268630237516536  -2.589933927753111 
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

-Ni 0
lanl2dz


