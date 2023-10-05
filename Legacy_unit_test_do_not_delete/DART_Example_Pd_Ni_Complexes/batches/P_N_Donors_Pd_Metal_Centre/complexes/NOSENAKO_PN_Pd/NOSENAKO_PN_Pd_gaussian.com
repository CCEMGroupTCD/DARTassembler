%chk=NOSENAKO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5076944642400167  -1.4896165286744096  3.0938944959541946e-15 
N  1.5076944642400203  1.4896165286744096  -5.154923187673577e-16 
C  1.111108392944973  2.5181091141879177  0.9917339106829651 
C  0.5182919761605735  1.9418308281488201  2.253497229626242 
C  1.742039572526123  2.115918027779162  -1.3338603309348225 
C  2.815581227659574  3.1513714953723193  -1.429745059346883 
C  2.747455971632868  0.8060638455512511  0.48085430660667516 
C  2.9926908644814905  -0.4874447708674075  -0.2343227669773627 
C  1.705327383732611  -2.8749365021553466  -1.1205821880645939 
C  1.7197044611472485  -2.60513451468838  -2.4966906384433685 
C  1.91449711975817  -3.6268503896823128  -3.408990584150213 
C  2.0473298169488623  -4.902504283433743  -2.9405062216429005 
C  1.9872361208621532  -5.173829414579854  -1.5921264713873755 
C  1.8131811286261401  -4.176392704259385  -0.6975426583308146 
C  1.6615401504168679  -2.182146404747959  1.67054559817755 
C  2.8606131380590822  -2.7126424383165575  2.1335030785953197 
C  2.946334165344144  -3.26243019189955  3.412382725587727 
C  1.8068638498233869  -3.296901731394144  4.206307969557186 
C  0.6174126491951817  -2.7863086053385158  3.7632295104228275 
C  0.5559053893008725  -2.2305128692796816  2.5054412743675103 
H  1.887426128977542  3.028010129233662  1.2099704511078095 
H  0.4748552631715004  3.0973928146500036  0.5734005397544437 
H  0.2831611218334078  2.653257762262569  2.8524464084059566 
H  -0.2724781686955793  1.4518805080825725  2.0479290611974603 
H  1.1394241682411934  1.3722747419559553  2.6813269558512878 
H  1.9468517019215494  1.417867747734721  -1.9431945973128018 
H  0.909555756370266  2.5192139824330795  -1.6019256379000988 
H  3.639029597239281  2.7906753634099974  -1.1989402651506849 
H  2.846931039570289  3.481811911093287  -2.3387936449303535 
H  2.5994927711983724  3.889516037161524  -0.8658450820361698 
H  3.478574241335526  1.394203325165935  0.3334587903059259 
H  2.6550233006454786  0.6398142402112921  1.4073806527288277 
H  3.7546778474515223  -0.9362469304372445  0.1337554982729401 
H  3.1505367165311577  -0.34134381413750997  -1.1651616184620486 
H  1.575438045727468  -1.7099457400400133  -2.823066861727918 
H  1.990558240132964  -3.440879178523344  -4.3291584245942385 
H  2.181060611070672  -5.613694505619815  -3.572790993379894 
H  2.0555835290447577  -6.081881308447939  -1.2823344123585712 
H  1.7745275341423907  -4.345527924734015  0.2249632206902886 
H  3.62050206718187  -2.710678916199493  1.5736822769297174 
H  3.7807055321806633  -3.5750968923788116  3.729633756499842 
H  1.8509676193743114  -3.7716545327941966  5.043189302288538 
H  -0.16250535957399093  -2.7767034219073596  4.333763435837276 
H  -0.2703149853220872  -1.864001103002595  2.1962382379361736 
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


