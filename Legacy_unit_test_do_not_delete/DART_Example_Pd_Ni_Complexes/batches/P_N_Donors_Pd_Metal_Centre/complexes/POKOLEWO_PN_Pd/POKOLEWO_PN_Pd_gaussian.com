%chk=POKOLEWO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5898296300232928  -1.4016211141032369  1.275843202386186e-16 
N  1.5898296300232928  1.4016211141032373  -1.7164908110038768e-16 
N  1.582059746856553  2.747635490611154  0.23779156240199478 
H  0.7975371815038228  3.1636696771203123  0.3210880215587424 
C  2.8368477554493285  3.197341244653879  0.46874183469597885 
C  3.694206147726028  2.113358005736938  0.3636239249952456 
H  4.628209104295599  2.112108099480266  0.46312058249524796 
C  2.863068080375606  1.0082664518108213  0.0734845439347712 
C  3.1375807614712263  4.672459252078326  0.7018675802951525 
C  1.902035056794639  5.394602272624637  1.2864932721113067 
H  1.1624157903061991  5.329364756522791  0.6633310403193041 
H  2.1156636639778754  6.329687444225059  1.4366392012801115 
H  1.6524793458927227  4.981850271260746  2.128965296885037 
C  3.520915156963933  5.3318081058456155  -0.6179078552498656 
H  2.778691286261154  5.2658044554517  -1.240289654384201 
H  4.295707352359217  4.883094920326582  -0.9912617103202643 
H  3.731938924490807  6.264086869950873  -0.4635907287306923 
C  4.310128273777502  4.768684432794365  1.704072697015525 
H  4.517100849494632  5.701351554348382  1.8703521403426009 
H  5.089949212724419  4.326073477695478  1.3346564201591553 
H  4.059626459456115  4.339412007496352  2.5373561157275875 
C  3.159053446311021  -0.44548409781875553  -0.1706664597571786 
H  3.8173906553679524  -0.7671834504205376  0.4796176594605176 
H  3.5266943841294696  -0.5668464176984992  -1.0693739806745137 
C  1.6169823039214553  -2.1554255369739526  1.6555644337689779 
C  2.8039191282794214  -2.7295704338385174  2.1325031639926935 
H  3.579843225784841  -2.7269442676349955  1.6026524131566762 
C  2.819557444919358  -3.3089031413937415  3.4120690071634923 
H  3.612418542893418  -3.692357217956984  3.739044098701718 
C  1.687970768716878  -3.3220001579949403  4.187101285457268 
H  1.7088702753671068  -3.7079946928633354  5.0444478954543674 
C  0.5287965743121064  -2.7712171844565305  3.712625666393151 
H  -0.24240351347005085  -2.7901483454764446  4.250077790537466 
C  0.4645720383296643  -2.1807943186241556  2.4474475846670445 
H  -0.3398983753918954  -1.811862504743167  2.1325251971599735 
C  1.6976198464403491  -2.8089469992173726  -1.1490546428836705 
C  0.7193280218724944  -3.801341831175896  -1.0294936815052345 
H  0.02351472338116878  -3.703384980627278  -0.4056860436285961 
C  0.775051523610276  -4.9440616032198195  -1.8402864256353069 
H  0.11825453362459681  -5.609508203313845  -1.7541517266513533 
C  1.7919374120097387  -5.091341497852453  -2.764697904691602 
H  1.812674814254752  -5.848968540792237  -3.321006419775145 
C  2.7720594826465508  -4.138294442868085  -2.8721612362807245 
H  3.4837278476067755  -4.257654788054549  -3.4754440773360913 
C  2.709979588052575  -2.9874416526459275  -2.081991544426446 
H  3.3651560512047225  -2.322230158884621  -2.1845414381334054 
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


