%chk=UGIKIDUV_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.4544958757246442  1.4016211141032375  -1.703669074007909e-16 
N  1.4544958757246444  -1.401621114103237  2.220446049250313e-16 
N  1.4963354289483113  -2.747635490611154  -0.23421069578490789 
H  0.7462748715894192  -3.1636696771203123  -0.4787983405058721 
C  2.7717205705903503  -3.197341244653879  -0.19922905403371466 
C  3.5884883826236997  -2.1133580057369366  0.08184661010925698 
H  4.522767651962047  -2.112108099480265  0.17871432722376346 
C  2.715189306877239  -1.008266451810821  0.19284242868611978 
C  3.1143514068328084  -4.672459252078327  -0.36473453499410574 
C  2.0273558562874436  -5.394602272624637  -1.1934691493133474 
H  1.174336332049636  -5.32936475652279  -0.7376999996750752 
H  2.267533259854721  -6.329687444225058  -1.295918144763251 
H  1.9584133195566846  -4.981850271260746  -2.0694166888379746 
C  3.2149122837588346  -5.331808105845615  1.005900343107163 
H  2.359507333231333  -5.2658044554517  1.4603645867476835 
H  3.8951487794710045  -4.883094920326582  1.5321838760460313 
H  3.4534090103598047  -6.264086869950871  0.898829724359546 
C  4.469646103105145  -4.768684432794365  -1.1012527295095387 
H  4.7066672716709155  -5.701351554348381  -1.2208665498810098 
H  5.1556201208012435  -4.326073477695479  -0.5777751946789105 
H  4.397867718962052  -4.339412007496351  -1.9684091620224793 
C  2.953944834426836  0.4454840978187562  0.4931969650530958 
H  3.733097466678155  0.7671834504205384  -0.006000883691716252 
H  3.126700135460931  0.5668464176985  1.448702419395482 
C  1.8252663992342937  2.1554255369739526  -1.613741020410958 
C  3.085426944019933  2.7295704338385174  -1.8334794530301968 
H  3.7342330728590283  2.726944267634996  -1.153883521294342 
C  3.3667602238871055  3.3089031413937424  -3.0818323236544964 
H  4.210277348722594  3.6923572179569857  -3.2368171334789064 
C  2.42103970295569  3.322000157994941  -4.075198386314904 
H  2.619734888543506  3.7079946928633363  -4.909464664115928 
C  1.188547217654691  2.7712171844565314  -3.852097064631664 
H  0.5459422819693216  2.7901483454764455  -4.5381460845823 
C  0.8626808277179829  2.1807943186241556  -2.6279191913583757 
H  0.01031397661194089  1.8118625047431678  -2.487137417530689 
C  1.3210287236516243  2.808946999217375  1.146355888197586 
C  0.3889730442690682  3.8013418311758973  0.8260093133247453 
H  -0.1619381634729331  3.703384980627279  0.07116564963347466 
C  0.2749055634635087  4.94406160321982  1.6306698581525454 
H  -0.34963042555634494  5.6095082033138475  1.40986163628953 
C  1.0773741017995668  5.091341497852454  2.7463031929081163 
H  0.9819952980636357  5.8489685407922405  3.294766580568803 
C  2.0137352706476133  4.1382944428680855  3.055197029782776 
H  2.584422418826447  4.25765478805455  3.793260866385343 
C  2.117297487363084  2.9874416526459275  2.2693873056721197 
H  2.7368354510588304  2.3222301588846217  2.505915084395059 
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

-Ni 0
lanl2dz

