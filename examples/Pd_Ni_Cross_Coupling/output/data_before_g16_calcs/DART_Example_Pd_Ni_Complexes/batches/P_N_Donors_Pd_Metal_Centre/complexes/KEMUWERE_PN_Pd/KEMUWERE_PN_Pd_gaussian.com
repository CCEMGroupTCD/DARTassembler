%chk=KEMUWERE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5492020099070365  -1.4464000596308062  2.6830139175038353e-15 
N  1.549202009907036  1.4464000596308066  4.49116845936905e-17 
N  1.0820612980241866  -2.4003919761910186  1.258103504660645 
N  2.0913700320375357  -2.7009041907099727  -0.9713220891570603 
C  1.3535195735858483  2.740755585637813  0.03531637595403396 
C  2.310873879757343  3.6589163608378485  0.4844533409646931 
C  3.5062166745871917  3.212320624271798  0.9030925542849412 
C  3.776117079213659  1.819731128310264  0.9180897980338169 
C  4.99075848453904  1.2612378541837457  1.37294578538002 
C  5.171475677941742  -0.09331240800311541  1.3907963657111568 
C  4.157760885713498  -0.9544149824304936  0.9714528116504755 
C  2.960832193247069  -0.4389187267193881  0.49680955034949614 
C  2.760233859472414  0.9493915624218855  0.46279911737212626 
C  -0.18300511318501056  -2.28746344102372  1.9719765544217733 
C  -0.8802193208955662  -3.6168309885408245  1.7304210904494248 
C  0.2840144984709234  -4.593115762312501  1.6877185806296549 
C  1.4149311775252715  -3.8425890579115833  0.9981855930084798 
C  1.4870906903453704  -3.977910875479788  -0.5106005417368159 
C  2.415999387028074  -2.5244372939475737  -2.347362766690059 
C  1.9944179496427692  -3.4269723887907033  -3.2997785919621125 
C  2.331994806142257  -3.2300864663214477  -4.648242547538781 
C  3.075728491880739  -2.1455628177256374  -5.023557615627257 
C  3.489867373766388  -1.2538900495917646  -4.081861315815579 
C  3.1802292470927114  -1.4380971527551916  -2.7480563137053737 
H  0.49171105730960507  3.0458258924743813  -0.2577866761631457 
H  2.1242419070580483  4.6005957295096565  0.49193220346108213 
H  4.1786148961719585  3.8342891308477918  1.190681198770908 
H  5.6933568741711715  1.8444920751600178  1.6694677349887066 
H  6.006917369392157  -0.44850809648444523  1.7025312900791636 
H  4.288185160869501  -1.9048649086776612  1.0086277514188435 
H  -0.7268128039244701  -1.5741162635745136  1.6297504742965592 
H  -0.0511216483768282  -2.165127121399143  2.915002713449676 
H  -1.3538822486106203  -3.6415932766515846  0.8957552890678553 
H  -1.4738246704085618  -3.8655036483662166  2.442694638522753 
H  0.06590830615898535  -5.375065896004206  1.1752635510956397 
H  0.5601133448052806  -4.854866138689445  2.5691202016398393 
H  2.26218662440785  -4.055463986911927  1.3963162152054223 
H  2.0375335041210416  -4.732071574476461  -0.7341413408129078 
H  0.6057575532823647  -4.1078046028403605  -0.8684663002372369 
H  1.4706443075096685  -4.189804281428783  -3.0441251441448878 
H  2.0364343565833884  -3.8660216309222997  -5.303995402721765 
H  3.3073109394431857  -2.0064783585813615  -5.944794856084659 
H  4.006947363742645  -0.49255232862881015  -4.354990531764953 
H  3.4905880846482535  -0.8085053403541045  -2.09303803021347 
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

