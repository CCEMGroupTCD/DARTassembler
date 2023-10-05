%chk=KIHOKERU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.300791176553718  1.5453291930847617  7.176011966234525e-15 
N  1.300791176553716  -1.5453291930847608  2.220446049250313e-15 
C  1.6637603426053362  -2.237668178064591  1.267296419300918 
C  0.9591977391872395  -2.636799720884131  -0.9579471177959379 
C  2.5025483590615423  -0.8518650594907876  -0.48812719624054374 
C  2.739360878118445  0.41714147567044635  0.22298909842300102 
C  1.6502227883195302  2.2278831141004423  -1.684052996079664 
C  1.5699785302364995  2.8770587230834326  1.2691614505439552 
C  0.5212041959610937  3.020036064732533  -2.2982229814370574 
C  2.9630321708764074  3.004342534765891  -1.8085728062728559 
C  0.8732983077900569  4.1875247513330685  0.9255631146114024 
C  3.04164650652688  3.1147522589615533  1.6213420584189024 
H  1.8928538134141142  -1.5926325769379348  1.9240037094325122 
H  0.9205266540137353  -2.749149196825522  1.568649951866603 
H  2.4028439981987457  -2.8120880259340613  1.1121264531079988 
H  0.7186809181518401  -2.25532176289433  -1.7918566370886282 
H  1.7140943238510762  -3.2016622006108837  -1.069135406117295 
H  0.23169834155104763  -3.1386254862098917  -0.6126769805204904 
H  2.3964620781794896  -0.6734988471436401  -1.4110598510434653 
H  3.253747674542902  -1.4191450500085852  -0.35675607773098655 
H  2.8614714305692908  0.24251459926923236  1.1449368230899295 
H  3.5172893708937973  0.8298277743916662  -0.1332086120919106 
H  1.7390236910073695  1.431803948246929  -2.1894468735414017 
H  1.1575077047037414  2.529987673497313  2.046790439602053 
H  -0.2851804614929594  2.5386517691863886  -2.2162824894671176 
H  0.4449365893119387  3.8559185716963085  -1.8498331435415312 
H  0.7134279210254826  3.176110787244277  -3.2162464781490265 
H  3.6647504314293413  2.4908859369415075  -1.4188829745384528 
H  3.1540472987683277  3.1507528245685186  -2.727840502842448 
H  2.8843261386243197  3.8311462576991953  -1.3589063474049925 
H  -0.033866261577548906  4.006774127938433  0.7124934634957331 
H  0.9156039703371202  4.764909300847018  1.6776491754442335 
H  1.3055874515106916  4.5871960496970825  0.1863572708267916 
H  3.449270647589143  2.291738984061076  1.8276377175124887 
H  3.47531318908609  3.520119806409187  0.8793586017876961 
H  3.0870483647742377  3.698344408152598  2.3716478725175913 
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

