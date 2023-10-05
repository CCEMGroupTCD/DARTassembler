%chk=MOKIPULA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.521253710759652  -1.475766630433145  -3.4500155054505903e-15 
N  1.521253710759649  1.4757666304331454  -1.5129969175750304e-15 
H  1.551108641814845  1.8630570313769907  -0.8126737896772983 
H  1.374547301790202  2.1267444342326756  0.6034598151721752 
C  2.8506925454484877  0.8711260411970928  0.2775967447026531 
H  2.9488590322390102  0.7068313281984435  1.2288491392604382 
H  3.555000734382754  1.4763348232508346  -0.002039867930900011 
C  2.951624564284591  -0.4424655423387301  -0.49283568532448885 
H  3.7827643789375704  -0.8945939386856256  -0.2787393945686403 
H  2.930909344444886  -0.2732835525503314  -1.448617513985326 
C  2.182132540209862  -2.298245415336789  1.5125909794680075 
C  3.297561611342547  -3.1409715275507373  1.4454212164426403 
H  3.681993298888633  -3.3295842521237975  0.61867746669933 
C  3.837172146454531  -3.698528598588788  2.590220711084228 
H  4.584863937222517  -4.247369379574787  2.5303279473741966 
C  3.275313216124835  -3.446152204977907  3.8015105023219737 
H  3.6369419811648287  -3.8284723620845886  4.56928226122277 
C  2.1735734388471784  -2.629145260833817  3.8938128662383273 
H  1.7930859867432836  -2.4483242695625362  4.724149929067151 
C  1.6309168071638016  -2.0752929828756823  2.743652792904779 
H  0.8755367259180725  -1.539697582017479  2.809985284466246 
C  1.4838324019033264  -2.865383848839888  -1.1840600035264264 
C  1.7710198407588449  -2.7049311457383833  -2.522927530592147 
H  2.0887794129972375  -1.8846977531153777  -2.8268222095140247 
C  1.5969732490613802  -3.7341187209813684  -3.425630493301358 
H  1.8295886141861246  -3.612019866539109  -4.317647966029201 
C  1.0873482881502226  -4.936352785915645  -3.018186555651141 
H  0.9510607332724167  -5.621424855679998  -3.631923406112478 
C  0.7738693157203691  -5.123698911175877  -1.6929516771505864 
H  0.413322967825017  -5.933662143403975  -1.4131674862856642 
C  0.9952497588610241  -4.100508320698296  -0.7662020373054939 
H  0.8165239995675555  -4.246359389262452  0.1344775776318532 
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