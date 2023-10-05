%chk=LIRUHUDI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  3.135581466050292  -0.9539552521268841  -2.2236554898534826 
C  2.236731976758807  -1.6478390128220912  -3.016716459920677 
H  1.4998320707055202  -2.0608156774732254  -2.6274187551852393 
C  2.4274072873065324  -1.7295463653611254  -4.388601587942929 
H  1.8187933812674228  -2.1959762074170928  -4.915445820914595 
C  3.4991026896616804  -1.1302928483549632  -4.962123720843334 
H  3.636986859159702  -1.1982417065752682  -5.8789062678349575 
C  4.384181746966338  -0.4188945503844426  -4.181997440062226 
H  5.104549989311097  0.010203755821133526  -4.582406254941059 
C  4.2146411122556975  -0.33356016545870204  -2.8162815999507806 
H  4.825070002829275  0.13767089157494  -2.29837627999348 
C  3.039406809154178  0.6327086287124571  0.29069483868286666 
H  3.777000985178793  1.1154367669902385  -0.11181119358175864 
H  3.2080113124722205  0.5569438508043095  1.2422212229147624 
C  4.255380090891358  -1.9992725518668624  0.19106586108299164 
C  4.239498244405947  -3.3703839900680284  -0.028885499417792535 
H  3.53509992380764  -3.7579366302233628  -0.4950278135739995 
C  5.273010035733012  -4.158446985388203  0.4456954763860859 
H  5.265981123533652  -5.074231357358906  0.28803816305266156 
C  6.308387178967631  -3.5945037278374685  1.1489226256988032 
H  6.991664833465407  -4.134902048008168  1.4750282564700312 
C  6.345332878142491  -2.2349233278911047  1.3752190052415598 
H  7.047714904881333  -1.8600748657558752  1.856445364608057 
C  5.32249737361412  -1.4213830262013922  0.8779282570842506 
H  5.354354473704554  -0.5005655865197495  1.0043924731770764 
C  1.7482084927696775  2.4283840954407028  -1.5677213918629798 
C  3.017281876361454  2.885048833255502  -1.9199561492272452 
H  3.727526838664696  2.77647871527657  -1.3305964432709287 
C  3.224755505813194  3.499296436567611  -3.1448683684933783 
H  4.08142702942454  3.7650807746541983  -3.3860707031862707 
C  2.1759238072633416  3.7190941051633843  -4.005863149544312 
H  2.3217813094784923  4.118980857099203  -4.833802588557698 
C  0.9041716466350769  3.3397440318635603  -3.6307752347499096 
H  0.18083160180226332  3.529478555199906  -4.184016411991446 
C  0.7007026753023666  2.6786256970868747  -2.4353809294569104 
H  -0.15574495921304243  2.3973419975898067  -2.209456672834896 
C  1.5221101838403226  2.8550957908331425  1.2610755332828913 
C  2.257680064438643  2.7514230168050644  2.414901805596845 
H  2.771503835972343  1.9910433538509735  2.5689872039080246 
C  2.2412237693966124  3.770860531610736  3.352991339575523 
H  2.753500368721353  3.694801070723107  4.126257202131095 
C  1.4833909771473568  4.878475754172208  3.1539709442845365 
H  1.4753509376043952  5.556102866495535  3.790605497547644 
C  0.73030829125723  4.999456475194032  2.014238522461954 
H  0.2053834836178854  5.755972097558779  1.878795178428798 
C  0.752052524003828  3.984959933557903  1.0536731236729913 
H  0.24989734325921664  4.0692642022580685  0.27590059604663064 
N  1.4622056387526345  -1.5342928892489858  1.878966875773315e-16 
P  2.8997759081637504  -0.9982861425824159  -0.44015850330675366 
P  1.4622056387526345  1.5342928892489858  -1.878966875773315e-16 
H  1.347674643565033  -2.2411767558169133  -0.0684270438300089 
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
