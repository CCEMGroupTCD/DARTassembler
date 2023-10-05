%chk=ZUSENEJE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.521253710759652  1.475766630433145  3.1287390105014262e-15 
N  1.521253710759649  -1.4757666304331454  1.3322676295501878e-15 
H  1.551108641814846  -1.8630570313769907  0.8126737896772965 
H  1.3745473017902048  -2.1267444342326756  -0.6034598151721796 
C  2.850692545448481  -0.8711260411970928  -0.2775967447026516 
H  2.9488590322390085  -0.7068313281984434  -1.2288491392604388 
H  3.5550007343827494  -1.4763348232508346  0.0020398679308968326 
C  2.951624564284592  0.44246554233873003  0.492835685324492 
H  3.7827643789375784  0.8945939386856256  0.2787393945686405 
H  2.9309093444448893  0.2732835525503312  1.448617513985325 
C  2.1821325402098637  2.298245415336789  -1.512590979468011 
C  3.297561611342533  3.1409715275507373  -1.4454212164426283 
H  3.681993298888637  3.3295842521237975  -0.618677466699333 
C  3.8371721464545225  3.6985285985887884  -2.5902207110842252 
H  4.584863937222503  4.247369379574787  -2.5303279473741873 
C  3.2753132161248413  3.4461522049779076  -3.801510502321985 
H  3.6369419811648287  3.828472362084589  -4.56928226122276 
C  2.173573438847172  2.6291452608338175  -3.893812866238332 
H  1.7930859867432805  2.4483242695625367  -4.724149929067142 
C  1.6309168071637978  2.0752929828756828  -2.743652792904768 
H  0.8755367259180694  1.5396975820174794  -2.809985284466249 
C  1.4838324019033273  2.865383848839888  1.1840600035264297 
C  1.7710198407588416  2.704931145738383  2.5229275305921455 
H  2.0887794129972397  1.8846977531153772  2.8268222095140247 
C  1.596973249061385  3.734118720981368  3.425630493301347 
H  1.8295886141861248  3.6120198665391086  4.317647966029206 
C  1.0873482881502334  4.936352785915645  3.0181865556511327 
H  0.9510607332724177  5.621424855679997  3.6319234061124686 
C  0.7738693157203708  5.123698911175877  1.6929516771505906 
H  0.41332296782501876  5.933662143403975  1.4131674862856662 
C  0.9952497588610218  4.100508320698296  0.7662020373054964 
H  0.8165239995675536  4.246359389262452  -0.13447757763185073 
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

