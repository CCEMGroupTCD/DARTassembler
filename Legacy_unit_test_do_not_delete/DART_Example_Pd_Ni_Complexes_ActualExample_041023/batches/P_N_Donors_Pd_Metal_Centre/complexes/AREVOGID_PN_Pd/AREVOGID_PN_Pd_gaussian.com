%chk=AREVOGID_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.361660904006574  -1.624185821419457  1.9890539674218528e-16 
N  3.980908048725577  -1.7211624160447856  -1.1073120178957392 
N  1.361660904006574  1.624185821419458  -1.9890539674218538e-16 
C  3.0066497612201966  -1.0109226271074299  -0.46953198415581215 
C  3.38985929919052  0.31415204514842715  -0.3564107529275427 
C  4.680368515601714  0.4157628047201164  -1.0015831375636481 
C  5.006414146097969  -0.8447097058764907  -1.4699139361073823 
C  6.175138513737111  -1.133931698708235  -2.184945824451469 
C  7.030475454741155  -0.06910212244182275  -2.387560979023305 
C  6.734385441797554  1.2046515070486514  -1.924778410644888 
C  5.572827919388746  1.475465350636813  -1.2671186093368318 
C  2.690312506324288  1.4231599531387884  0.24560278895121623 
C  0.7734374796117134  2.7690741774665995  0.41740381548432315 
C  1.3957073969215616  3.770591933962213  1.087068905984946 
C  2.7589787530630647  3.5529098342065266  1.4476782277863989 
C  3.454655392411965  4.542741877312885  2.222822781985244 
C  4.730373982217611  4.25084793651634  2.6772523881406993 
C  5.315651801719648  3.021993479915007  2.407347648482418 
C  4.679467212840503  2.0856370329224556  1.6241783481084777 
C  3.378254567177441  2.3725822142273216  1.1026236877971298 
C  1.0445971187137282  -3.1198321526237227  -0.9843730481941387 
C  1.0082435090089783  -4.388730161595907  -0.42591030289426146 
C  0.6665474293939517  -5.474913784636161  -1.1812360182084374 
C  0.3738386080481795  -5.329382462491369  -2.5241948523434092 
C  0.4255801262601503  -4.066735579385513  -3.1023048798262454 
C  0.7371461908518744  -2.9543811700722635  -2.34504987952832 
C  1.5160890023718108  -2.208197873791038  1.7260360990199781 
C  2.7256606019896497  -2.167460290610946  2.424681398663926 
C  2.761964085065362  -2.5103385885592058  3.759149434966308 
C  1.6029559710920882  -2.9172039430556946  4.411268929426188 
C  0.40822216005748124  -2.950796192341751  3.7214215973601354 
C  0.34992478163763097  -2.61918520738112  2.3925763918064966 
C  4.122855715416309  -3.191576812910292  -1.2701618610387502 
H  6.36284709778871  -1.9874638264474414  -2.502846318931007 
H  7.82659805461553  -0.20945935065648982  -2.8472727848943222 
H  7.344454966889337  1.8910760345239919  -2.0671345929292317 
H  5.373252449247488  2.342794560399809  -0.9984643080316228 
H  -0.13064482160616153  2.8773287704332695  0.230520521083182 
H  0.9552478335251613  4.560661621849606  1.301004015993049 
H  3.057243718278137  5.362062938025635  2.416222887297213 
H  5.199539455359306  4.885593097078698  3.1700270830616217 
H  6.154150654998486  2.8269342411124403  2.76091733950229 
H  5.091047422976981  1.2739604515766976  1.4368043960263288 
H  1.2181682439013675  -4.502473766276893  0.47321294676317727 
H  0.630720464706652  -6.316699319252589  -0.7901261841313182 
H  0.1430927088768703  -6.0710085800911955  -3.036311413010442 
H  0.2490811453837154  -3.968836073120256  -4.010270273672609 
H  0.7421280087500896  -2.108010473523139  -2.730303673167503 
H  3.5056996105937674  -1.9088707225571249  1.9895226341311798 
H  3.5653656821334714  -2.4714460251794854  4.226137585075234 
H  1.6324510743878906  -3.1640438980120944  5.306783912840051 
H  -0.3694409416017195  -3.2013878146995287  4.165124244548463 
H  -0.4580750284324515  -2.666202647745955  1.9335187404943213 
H  5.051651479537032  -3.4162083083885277  -1.3595663899544803 
H  3.647961794351609  -3.4734640815726134  -2.0550119623331957 
H  3.7589187083856084  -3.6353395249687543  -0.49875607268863636 
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


