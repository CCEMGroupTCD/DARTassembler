%chk=HIJUJEVI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  3.135581466050292  0.9539552521268839  2.2236554898534826 
C  2.236731976758807  1.6478390128220908  3.016716459920677 
H  1.4998320707055202  2.060815677473225  2.6274187551852397 
C  2.4274072873065324  1.729546365361125  4.388601587942929 
H  1.8187933812674228  2.1959762074170923  4.915445820914595 
C  3.4991026896616804  1.1302928483549626  4.962123720843334 
H  3.636986859159702  1.1982417065752675  5.8789062678349575 
C  4.384181746966338  0.4188945503844421  4.181997440062226 
H  5.104549989311097  -0.010203755821134086  4.582406254941059 
C  4.2146411122556975  0.3335601654587017  2.8162815999507806 
H  4.825070002829275  -0.13767089157494028  2.29837627999348 
C  3.039406809154178  -0.6327086287124571  -0.2906948386828667 
H  3.777000985178793  -1.1154367669902385  0.1118111935817585 
H  3.2080113124722205  -0.5569438508043094  -1.2422212229147624 
C  4.255380090891358  1.9992725518668624  -0.1910658610829914 
C  4.239498244405947  3.3703839900680284  0.028885499417792948 
H  3.53509992380764  3.7579366302233628  0.49502781357399994 
C  5.273010035733012  4.158446985388203  -0.4456954763860854 
H  5.265981123533652  5.074231357358906  -0.28803816305266094 
C  6.308387178967631  3.5945037278374685  -1.1489226256988028 
H  6.991664833465407  4.134902048008168  -1.4750282564700308 
C  6.345332878142491  2.2349233278911047  -1.3752190052415596 
H  7.047714904881333  1.8600748657558754  -1.8564453646080568 
C  5.32249737361412  1.4213830262013922  -0.8779282570842504 
H  5.354354473704554  0.5005655865197496  -1.0043924731770764 
C  1.7482084927696775  -2.4283840954407028  1.5677213918629795 
C  3.017281876361454  -2.8850488332555027  1.9199561492272448 
H  3.727526838664696  -2.77647871527657  1.3305964432709283 
C  3.224755505813194  -3.4992964365676116  3.144868368493378 
H  4.08142702942454  -3.7650807746541988  3.3860707031862702 
C  2.1759238072633416  -3.719094105163385  4.005863149544311 
H  2.3217813094784923  -4.118980857099204  4.833802588557697 
C  0.9041716466350769  -3.3397440318635607  3.630775234749909 
H  0.18083160180226332  -3.5294785551999066  4.184016411991446 
C  0.7007026753023666  -2.678625697086875  2.43538092945691 
H  -0.15574495921304243  -2.397341997589807  2.2094566728348957 
C  1.5221101838403226  -2.8550957908331425  -1.2610755332828918 
C  2.257680064438643  -2.751423016805064  -2.4149018055968456 
H  2.771503835972343  -1.9910433538509733  -2.568987203908025 
C  2.2412237693966124  -3.7708605316107358  -3.3529913395755235 
H  2.753500368721353  -3.6948010707231065  -4.126257202131096 
C  1.4833909771473568  -4.878475754172208  -3.153970944284537 
H  1.4753509376043952  -5.556102866495534  -3.790605497547645 
C  0.73030829125723  -4.999456475194032  -2.0142385224619543 
H  0.2053834836178854  -5.755972097558779  -1.8787951784287986 
C  0.752052524003828  -3.984959933557903  -1.0536731236729917 
H  0.24989734325921664  -4.0692642022580685  -0.27590059604663114 
N  1.4622056387526345  1.5342928892489858  0.0 
P  2.8997759081637504  0.9982861425824159  0.44015850330675377 
P  1.4622056387526345  -1.5342928892489858  0.0 
H  1.347674643565033  2.2411767558169133  0.06842704383000918 
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

-Pd 0
lanl2dz