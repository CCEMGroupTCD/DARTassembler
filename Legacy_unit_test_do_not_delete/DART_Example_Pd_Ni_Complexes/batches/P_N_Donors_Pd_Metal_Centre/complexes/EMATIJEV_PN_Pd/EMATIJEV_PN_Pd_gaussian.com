%chk=EMATIJEV_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5336702554656259  -1.4628586902021679  3.298297106091237e-15 
N  1.5336702554656263  1.4628586902021679  -6.2323773110616005e-16 
N  2.877204818926815  -0.4451113746400297  0.1964684922608832 
H  3.66796302721829  -0.8041556647808725  0.3408715868846586 
C  2.7648149440839145  0.9394775016239396  0.13607412365850283 
C  3.9188247418962323  1.7462404158442448  0.23837704048883124 
H  4.779248087709959  1.353250666492207  0.3309269345079501 
C  3.7720619808903235  3.114331328113713  0.20106474026995572 
H  4.532965247308703  3.67853730407571  0.2776181584995804 
C  2.4955634343731816  3.6710997094038023  0.05033411851293109 
H  2.373455685078905  4.612051822016792  0.011979596857209477 
C  1.4184587409897491  2.810360901593194  -0.04083426537198891 
H  0.5508513706679091  3.1850773172466056  -0.13733954400015272 
C  1.7297143083018591  -2.703062988232671  1.3443147096176027 
H  2.6394995984976743  -3.112273439665455  1.2745387669044583 
C  1.6120413108659288  -1.994604941709548  2.6948446655945517 
H  0.7464199613340923  -1.5395286452369208  2.750094404267777 
H  2.33336894454254  -1.3375236655494012  2.782411162694007 
H  1.681493234355451  -2.653486237938369  3.4169598972601416 
C  0.6772561807523585  -3.813285185959473  1.209477172448954 
H  0.77196327259735  -4.4449635026616106  1.952747648764703 
H  0.8075430443177399  -4.28559592904602  0.36046183386704894 
H  -0.21908869789924368  -3.417482666633439  1.2291116882364144 
C  1.853442404853022  -2.3564882928122035  -1.5855345678398929 
H  1.0286683286528873  -2.87628413850316  -1.807330920448897 
C  2.0573490215339  -1.321037346562115  -2.7007335274805335 
H  2.8855227063427322  -0.8232396543603061  -2.537099223531961 
H  1.3000042270594676  -0.6995979648304692  -2.712301462277781 
H  2.1166272924116476  -1.7796236300446917  -3.5651744608030023 
C  3.020103278074575  -3.3486021516015505  -1.534182468641593 
H  3.1532863193601584  -3.7432271887746356  -2.421897720031412 
H  2.8185219322967976  -4.057293053835399  -0.8885990610749146 
H  3.8353991533822116  -2.879250504994376  -1.2596338981964008 
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

-Pd 0
lanl2dz


