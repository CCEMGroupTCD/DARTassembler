%chk=CEPEVIJA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
N  1.379458209769328  1.4755321235066357  2.92086356288683e-15 
N  1.2162706842325688  2.646783984084804  -0.6664343933211899 
H  0.4849779963957711  2.7378648781290558  -1.1792448704533836 
C  2.351367566913802  3.364242813435215  -0.6318743714098198 
H  2.489750347905847  4.211348416957767  -1.0390907787482366 
C  3.2634138069163905  2.6638941774840497  0.08536384858302717 
H  4.152336881199696  2.926332547229416  0.28944108908792276 
C  2.6299774204882977  1.4754446330359487  0.4649494445294783 
C  3.206547872301053  0.3579920891485928  1.2014903980793212 
C  4.316820163245059  0.6283054314647931  2.0188608140645368 
H  4.592280102158746  1.5290711442717397  2.141999283622442 
C  5.014864644481096  -0.37624908238570454  2.6461532665287075 
H  5.761541730349633  -0.16508347979650484  3.191631549984624 
C  4.6287989872718205  -1.692760169506923  2.480725923533387 
H  5.123496985482413  -2.391423094953807  2.893627252148888 
C  3.5074142326457505  -1.986805787415186  1.7039892650182604 
H  3.2331453527482736  -2.891480839034959  1.6062348150094115 
C  2.7828179892624623  -0.9785153196588534  1.0693220871159301 
P  1.3794582097693284  -1.4755321235066352  -2.220446049250313e-16 
C  0.9064145023326945  -3.0911292304911244  0.6724965864284882 
C  0.39373256837947257  -3.1067527398900348  1.9651589555535598 
H  0.25744688187113374  -2.2912094820853937  2.4322390089730392 
C  0.07923499910765197  -4.320378943936517  2.575988895695546 
H  -0.27258536862031835  -4.3323549658193885  3.458827031341963 
C  0.27992710344522287  -5.506760176925056  1.8919994153165254 
H  0.07989957723311902  -6.334739094182773  2.3127950442555667 
C  0.7658273213410692  -5.4922620411954455  0.6078671821297931 
H  0.890839661919641  -6.310561053981777  0.1421609867920181 
C  1.0782925088997417  -4.285335067780764  -0.019342310311138883 
H  1.405423754546283  -4.279936164650774  -0.9110488595563261 
C  2.113562837306789  -1.7061325046026221  -1.6422094568192398 
C  3.1643758647140943  -2.6076630296302046  -1.8421198214194352 
H  3.497675695717523  -3.1217128750496417  -1.115990329470961 
C  3.7166960460863434  -2.746855180534995  -3.1087011501841895 
H  4.426692332539215  -3.3613699346056283  -3.246720026620178 
C  3.2420542131553782  -2.0011766673588944  -4.17069897759226 
H  3.6184032657827823  -2.1157434004763394  -5.0362358192398355 
C  2.2207689412964795  -1.0873502629594118  -3.975126008840432 
H  1.9036575316291287  -0.5644055789994153  -4.702420344977329 
C  1.6599425515389679  -0.9355001836351371  -2.710838425273769 
H  0.9625858509336855  -0.30268630237516536  -2.576955103662712 
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