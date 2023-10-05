%chk=BERUYOSO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5842275901839362  -1.4079499076671722  -2.1794196824537714e-15 
N  1.5842275901839356  1.4079499076671718  -2.1708255791037233e-15 
N  2.9515549847908984  -0.4567420431636664  -0.20581526781087872 
H  3.73616167456081  -0.8286028118947686  -0.35353869536049964 
C  1.9127120093754768  -2.465167775926635  1.4379150654268051 
C  0.8713490661927107  -2.917679283282482  2.2462790281788956 
H  -0.008911890229316821  -2.5830338939584627  2.1191815380260466 
C  1.1217970868317575  -3.862605042185614  3.243665657886423 
H  0.4130541025667125  -4.180479589825981  3.789440679500229 
C  2.415087828140557  -4.334185335194841  3.431137057545915 
H  2.583658364802176  -4.983717118974878  4.102542066485492 
C  3.4634730527693054  -3.869238616256069  2.650057945465566 
H  4.345704308650077  -4.18853355533871  2.7960665857005966 
C  3.2149747429379727  -2.932044300737906  1.6504979747927662 
H  3.9289695545074173  -2.6101585308450677  1.1132522207914826 
C  1.7363628421399813  -2.6156912773991774  -1.345465933775422 
C  1.1620053500031906  -3.8869008532892115  -1.2342206849248536 
H  0.6925229046789738  -4.1296152873319745  -0.4435909883919509 
C  1.2789474551370197  -4.792234241258949  -2.2790051510584637 
H  0.8938587382555471  -5.656647128658137  -2.1994236819291455 
C  1.9587852154422878  -4.4432971645031385  -3.4426284690304993 
H  2.030270476410237  -5.064083294918949  -4.158301535274952 
C  2.5298855686467174  -3.1844009599238245  -3.55246590648048 
H  3.001639387580581  -2.9458012603389085  -4.342795899381743 
C  2.4160090658085904  -2.270934981485032  -2.511684321162844 
H  2.8030862690230167  -1.4074258513808107  -2.595621148469813 
C  2.8311325024554206  0.9066194724320742  -0.14053392636750484 
C  3.96111096279391  1.7341274637713642  -0.22694465163169752 
H  4.826939837112754  1.3606060926372432  -0.34323323964368774 
C  3.790484629573065  3.103736080170976  -0.13991895257937645 
H  4.534193988383857  3.68853018791322  -0.22490072219597626 
C  2.5054543153373636  3.616062325637455  0.07606241166657104 
H  2.3699559195932354  4.551703798268176  0.17103519431605158 
C  1.4418019955253567  2.7485409434145778  0.14920695312165316 
H  0.57444296269843  3.101035646739861  0.30976062732804416 
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


