%chk=EPAWOZIY_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3642370019538417  -1.4896165286744105  -4.76752596856313e-15 
N  1.3642370019538403  1.48961652867441  -6.265146212298734e-16 
C  1.1458887505264252  2.5181091141879177  1.0455336927341743 
C  0.7811814481457121  1.9418308281488197  2.3910694821615883 
C  1.363399465745435  2.115918027779162  -1.354289596339089 
C  2.403981402547574  3.1513714953723193  -1.6351361723382156 
C  2.6686632203859633  0.8060638455512505  0.2582667227151879 
C  2.7859832486353704  -0.48744477086740773  -0.48862979637046416 
C  1.3642803784371402  -2.8749365021553475  -1.1378766230102821 
C  1.1394803110603908  -2.6051345146883818  -2.4955754472221665 
C  1.172894408432648  -3.6268503896823145  -3.4278408970175565 
C  1.385060534294417  -4.9025042834337444  -2.9895400204571194 
C  1.5600234829159323  -5.173829414579855  -1.6512100275860768 
C  1.5439556260190352  -4.176392704259387  -0.7399926205531323 
C  1.805832625301331  -2.18214640474796  1.6184512337991388 
C  3.067080722770723  -2.7126424383165593  1.8661585106358527 
C  3.373574575192781  -3.262430191899552  3.1107238019876986 
C  2.3892790459597424  -3.296901731394145  4.090454481356671 
C  1.1409585147592516  -2.7863086053385167  3.8606534130212893 
C  0.8619730532668857  -2.2305128692796825  2.632654430100395 
H  1.9483088533172923  3.028010129233662  1.1256485695825487 
H  0.4466589079282175  3.0973928146500036  0.7440399422200812 
H  0.6536291932421363  2.653257762262569  3.0217493214513347 
H  -0.0332716591932094  1.451880508082573  2.3259401507280977 
H  1.4671690988628971  1.372274741955955  2.7045410402720482 
H  1.459290253738443  1.4178677477347215  -1.989931959077724 
H  0.4970138972877106  2.5192139824330795  -1.4737230913394097 
H  3.254998573023908  2.7906753634099974  -1.5508281303617462 
H  2.277000310075454  3.481811911093287  -2.5358181049941146 
H  2.289096018734573  3.889516037161525  -1.04227973607987 
H  3.3630791979561128  1.3942033251659345  -0.013846879689896455 
H  2.7385244209304718  0.639814240211292  1.1867678166202043 
H  3.6003100571581754  -0.9362469304372449  -0.25846111803535443 
H  2.779792597252808  -0.3413438141375104  -1.4327367586782762 
H  0.9407309902485315  -1.7099457400400138  -2.7919416821754433 
H  1.0880145208397243  -3.440879178523345  -4.3472371953037 
H  1.4069645587833735  -5.613694505619818  -3.6354410744724728 
H  1.6811273668984876  -6.081881308447941  -1.3579928089120574 
H  1.6660807312686932  -4.345527924734017  0.17521044754857398 
H  3.7182133695697233  -2.7106789161994937  1.1828893170866077 
H  4.250360029525609  -3.5750968923788133  3.278268009532414 
H  2.5780356985028243  -3.7716545327941984  4.906963166970917 
H  0.4719613894549066  -2.7767034219073605  4.557950987078295 
H  -0.0053877211638466704  -1.8640011030025958  2.471620544972477 
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