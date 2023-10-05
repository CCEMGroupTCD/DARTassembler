%chk=ZEBOXAJI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5720736043519077  1.4215078552368252  -9.192393952598914e-16 
N  1.5720736043519083  -1.4215078552368257  -1.1102230246251565e-15 
C  2.807400693093374  -0.9620799156766914  0.30025686032643906 
C  1.8746344729709732  2.9781408044310105  -0.9682321169345216 
C  3.106878873412841  0.44747484521921854  -0.11456326400243988 
H  3.4411534364410716  0.4617809779114028  -1.0258619752270863 
H  3.785282410679124  0.8277757440207483  0.4648877925979795 
C  1.7871150459238998  2.6404097512881375  -2.4530025185900692 
H  2.4644841306327643  1.9959438226440966  -2.67379416222194 
H  1.9221601754787145  3.436640504308798  -2.9702892871504813 
H  0.9218537622447337  2.273840039358415  -2.646852943561596 
C  0.13413450589198583  2.348103311707627  2.193774817333801 
H  -0.5421226200814357  2.52472518831215  1.5808611674002397 
C  3.281033324018284  3.4847392184651174  -0.6003142746504023 
H  3.3547142697851937  3.559023493160459  0.3540655526531181 
H  3.4254428750264685  4.345936548814058  -0.9999124225158111 
H  3.939445274646264  2.866822668258182  -0.9252026404787097 
C  -0.06569656049077732  2.6045462983974668  3.5380140332049956 
H  -0.8786260842094118  2.9457555067835806  3.830812266770584 
C  1.3467441693132123  1.825535794923668  1.7574863753832097 
C  2.333112503190434  1.551821531884878  2.6983218562901587 
H  3.145276213257947  1.1892239119694217  2.4267049678589823 
C  1.263975801374722  -2.6913588225389153  0.3232367625171716 
H  0.4164499056333253  -3.0210562619680617  0.13130507028886657 
C  3.736889546815898  -1.7424002571462067  0.9462711987853462 
H  4.576481343118865  -1.3986883999095114  1.152019770531711 
C  0.9402511915272955  2.3534882274306046  4.438690787214021 
H  0.8132587670041189  2.554483437515038  5.338185395111609 
C  2.1100370436557645  1.8185482746202053  4.033139684333846 
H  2.7717910339162546  1.6265385337000422  4.658996389791356 
C  0.8492075728887096  4.047959000930979  -0.6089033830079645 
H  -0.03525772806205629  3.704723345776112  -0.7513391455274296 
H  0.986492278018988  4.821062161037618  -1.160542425834364 
H  0.9528956816112173  4.291838821026833  0.31465305785506237 
C  2.190355137889549  -3.513568774242122  0.9369295685680361 
H  1.9741846693277678  -4.400561238163425  1.1142471014535338 
C  3.412483280469644  -3.032811754516651  1.28251675163246 
H  4.018786803086351  -3.5696478945271917  1.7409131927630097 
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

