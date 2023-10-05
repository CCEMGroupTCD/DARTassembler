%chk=LELECEWO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5436714263404625  -1.4523011146108789  9.99830105446962e-17 
N  1.5436714263404625  1.4523011146108789  -6.683328867811896e-17 
C  -0.030106509954346716  3.24793519404813  -0.5459913827115307 
H  -0.5449025916691383  3.150375880366846  0.2809053064614066 
H  -0.46778049303358515  2.7457972195858624  -1.2651319664213543 
H  0.02091707115807595  4.195466219918551  -0.7915414574079461 
C  1.339182632292836  2.7173850056277025  -0.33979583840695465 
C  2.4332164758918413  3.602640826246192  -0.5743891271453673 
H  2.2568484586389053  4.515201180753044  -0.7737536723197707 
C  3.717127017972791  3.1746127102128603  -0.520165222652659 
H  4.438761928562776  3.7800213208339284  -0.6388525850956276 
C  3.952691459512118  1.810039279081821  -0.28151134161669134 
C  5.240646650876964  1.2402654273129952  -0.3382456599511961 
H  5.9920342587857025  1.810527827904695  -0.45404495680814694 
C  5.442782577797482  -0.10891076988698692  -0.23062658191676155 
H  6.323449514611923  -0.4671337270726888  -0.2622034618607401 
C  4.330624713114884  -0.9697235155516224  -0.07304333113809712 
H  4.460731924924896  -1.9092116260221395  -0.016409209624988253 
C  3.059962318822109  -0.44595772597305494  -0.0004043787958392673 
C  2.853599250153912  0.9516638275598317  -0.07602988368779301 
C  1.82073668916971  -2.7598608320106686  1.2219676245984483 
C  2.1869964470910848  -2.360268730096259  2.5123546021424237 
H  2.384024067804203  -1.4472357566585983  2.6839592526962828 
C  2.264638138663027  -3.2824380010045413  3.537707138357337 
H  2.5197844298069385  -3.0033201163442307  4.409056957074424 
C  1.972280310514045  -4.612358064115904  3.294383062167001 
H  1.9979347067668078  -5.241894174638662  4.005386358573414 
C  1.639059897454079  -5.026519231485906  2.0119712976533823 
H  1.463096764363815  -5.944997167349365  1.8426882296133151 
C  1.5621647311842073  -4.102691027746298  0.9712015424339072 
H  1.3328672481581754  -4.388381547312669  0.09467006478587353 
C  1.453023035065082  -2.2162594348503237  -1.6389951293692542 
C  2.4287147841505714  -3.080755187397064  -2.126612371266197 
H  3.154873475432095  -3.339908248847989  -1.5718353838046282 
C  2.3383454021331738  -3.567125863487372  -3.4204864030920663 
H  3.0092103528807916  -4.148214982341428  -3.7546484952087824 
C  1.2663173843945292  -3.202982811692155  -4.233825512573851 
H  1.2095035806984435  -3.5308579886850344  -5.123634227984906 
C  0.2918297486981605  -2.367230740176869  -3.742491967176383 
H  -0.4479776413395862  -2.1348186931128983  -4.2910882449846675 
C  0.3796109342355509  -1.862370905929305  -2.4451970773270584 
H  -0.2918491532260996  -1.2778006529294863  -2.115669369457541 
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


