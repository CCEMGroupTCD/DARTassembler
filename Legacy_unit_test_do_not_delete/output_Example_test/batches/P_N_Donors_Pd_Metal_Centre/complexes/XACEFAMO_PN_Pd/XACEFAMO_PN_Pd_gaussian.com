%chk=XACEFAMO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  2.3829943284474915  -1.0539917766477398  -1.5881081168897482 
C  2.664174597012761  0.29863304798093115  -1.8967140150314057 
C  3.0790074536350645  0.6249430684781834  -3.1860517565055715 
H  3.2782511647440753  1.5289892181190716  -3.402712927237653 
C  3.2057909750865052  -0.3638269529529172  -4.1605429896152035 
H  3.4371352789092926  -0.12432779514535741  -5.050637741031183 
C  2.9958463532091786  -1.6870034414288075  -3.8351228671744084 
H  3.131397194417026  -2.366013629318437  -4.485841649898513 
C  2.5836066994130396  -2.0246144428346984  -2.5442907597460818 
H  2.4387550801833022  -2.9375198035823074  -2.322803128658464 
C  2.633597494060358  1.3468853449757792  -0.8503693621468562 
C  1.6729966574537698  2.2995866732022083  1.0372688823021186 
C  3.7154847334041787  3.0409206999295018  0.3936285709349787 
H  4.4615875508653975  3.6165519469515264  0.5185121696512439 
C  3.722408170134295  2.181007833902722  -0.694205873668367 
H  4.447355430532653  2.1671801900859675  -1.3086478855030306 
C  2.7297772599963883  -1.0572948719051374  1.3467441482078142 
C  2.202664317507645  -0.6475434809962733  2.5630080848223753 
H  1.2634140755262258  -0.5580915088796724  2.670693495265567 
C  3.0569496297813163  -0.36764202792215417  3.6244320867181856 
H  2.70461372831583  -0.08498083791685462  4.459935039557586 
C  4.415293849382422  -0.5038672576344739  3.4561862907344927 
H  5.00181599117432  -0.2992410417958802  4.175666067364391 
C  4.938707561915801  -0.9368692998083459  2.2460772638260287 
H  5.8764179190109385  -1.0464223358918396  2.150197190393506 
C  4.11433813355729  -1.2052753667606346  1.1901547151459004 
H  4.476147641862395  -1.4899521631747412  0.3584959923013412 
C  1.2720358796990423  -3.192643219236354  0.07495939294507026 
C  1.9244313572905432  -3.978475071187405  1.0320779912793463 
H  2.537448071389415  -3.577058856852344  1.6370906366879763 
C  1.6771119704123358  -5.339271844169936  1.0977519854819082 
H  2.1289241185052012  -5.873343384538781  1.7418761215021044 
C  0.7774831112363416  -5.921240202395022  0.23267392208963023 
H  0.5988902642510072  -6.85186270134784  0.2899493046591619 
C  0.12727203855250546  -5.145658991396184  -0.7275208321685577 
H  -0.4874186937170859  -5.553893331885662  -1.3268792031846015 
C  0.36987070076639905  -3.792906901389684  -0.812260445544355 
H  -0.07324100135990963  -3.2690133759206983  -1.4692686900972096 
N  1.5792397071059223  1.413542340186525  3.2426370845947084e-15 
N  2.726499386091694  3.100024277958143  1.2637554639025765 
N  0.661875229944141  2.359198177296866  1.901475812042043 
H  -0.1025482933844939  1.934819290336293  1.6804841215105588 
H  0.7796503932126335  2.814844362550156  2.6411500205833627 
P  1.579239707105922  -1.413542340186526  -5.551115123125783e-16 
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