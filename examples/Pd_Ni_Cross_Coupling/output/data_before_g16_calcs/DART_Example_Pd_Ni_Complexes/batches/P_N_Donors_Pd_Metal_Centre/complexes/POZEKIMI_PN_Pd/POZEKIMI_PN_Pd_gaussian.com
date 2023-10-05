%chk=POZEKIMI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5720736043519077  -1.4215078552368252  1.0933238997477511e-15 
N  1.5720736043519083  1.4215078552368257  9.361385201372968e-16 
C  2.807400693093374  0.9620799156766914  -0.30025686032643917 
C  1.8746344729709732  -2.9781408044310105  0.9682321169345219 
C  3.106878873412841  -0.44747484521921854  0.11456326400243994 
H  3.4411534364410716  -0.4617809779114027  1.0258619752270863 
H  3.785282410679124  -0.8277757440207484  -0.4648877925979794 
C  1.7871150459238998  -2.640409751288137  2.4530025185900697 
H  2.4644841306327643  -1.9959438226440964  2.6737941622219403 
H  1.9221601754787145  -3.4366405043087975  2.970289287150482 
H  0.9218537622447337  -2.2738400393584146  2.6468529435615964 
C  0.13413450589198583  -2.3481033117076273  -2.1937748173338005 
H  -0.5421226200814357  -2.52472518831215  -1.5808611674002395 
C  3.281033324018284  -3.4847392184651174  0.6003142746504028 
H  3.3547142697851937  -3.559023493160459  -0.3540655526531177 
H  3.4254428750264685  -4.345936548814058  0.9999124225158117 
H  3.939445274646264  -2.866822668258182  0.92520264047871 
C  -0.06569656049077732  -2.604546298397467  -3.538014033204995 
H  -0.8786260842094118  -2.945755506783581  -3.8308122667705837 
C  1.3467441693132123  -1.8255357949236681  -1.7574863753832095 
C  2.333112503190434  -1.5518215318848783  -2.6983218562901587 
H  3.145276213257947  -1.189223911969422  -2.4267049678589823 
C  1.263975801374722  2.6913588225389153  -0.3232367625171719 
H  0.4164499056333253  3.0210562619680617  -0.13130507028886693 
C  3.736889546815898  1.7424002571462065  -0.9462711987853464 
H  4.576481343118865  1.3986883999095112  -1.1520197705317112 
C  0.9402511915272955  -2.353488227430605  -4.438690787214021 
H  0.8132587670041189  -2.5544834375150383  -5.338185395111609 
C  2.1100370436557645  -1.8185482746202057  -4.033139684333846 
H  2.7717910339162546  -1.6265385337000429  -4.658996389791356 
C  0.8492075728887096  -4.047959000930979  0.608903383007965 
H  -0.03525772806205629  -3.704723345776112  0.7513391455274301 
H  0.986492278018988  -4.821062161037618  1.1605424258343646 
H  0.9528956816112173  -4.291838821026833  -0.31465305785506187 
C  2.190355137889549  3.513568774242122  -0.9369295685680366 
H  1.9741846693277678  4.400561238163425  -1.1142471014535342 
C  3.412483280469644  3.032811754516651  -1.2825167516324605 
H  4.018786803086351  3.5696478945271917  -1.7409131927630102 
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


