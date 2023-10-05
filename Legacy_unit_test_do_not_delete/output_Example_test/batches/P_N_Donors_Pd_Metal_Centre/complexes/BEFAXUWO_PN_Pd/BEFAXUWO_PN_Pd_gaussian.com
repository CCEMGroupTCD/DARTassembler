%chk=BEFAXUWO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.47682973798607  1.5202216696916258  -1.7597714260673189e-15 
N  1.4768297379860698  -1.5202216696916262  4.440892098500626e-16 
C  1.5238647249797102  2.924049813670721  -1.1273266395196877 
C  0.6973394215977784  4.02595453151348  -0.9234063470883118 
C  0.6974278553817014  5.064101047554219  -1.8339491815623168 
C  1.4814091982448918  5.015566934752783  -2.9577539394882693 
C  2.3377611039486554  3.933354986455978  -3.157014980741096 
C  2.3445508646319877  2.8927190966117653  -2.232261816962792 
C  1.606792977447453  2.0157356924280645  1.7699313403802759 
C  2.929605649013599  2.7554309667548003  1.7304751991250544 
C  3.800876943514564  2.1680333340259255  0.947775133682911 
C  3.195433506765079  0.9307166633712083  0.3095572700838972 
C  2.880442656672577  -0.08233746926222163  1.4592588394976054 
C  1.854512953099074  0.6106939293848657  2.3737962820415817 
C  3.1304562907401916  3.942332539053801  2.6649920301089454 
C  5.265657434768546  2.464886446961451  0.7216067853617293 
C  2.500292532947574  -1.418377949735046  0.8618653862611177 
C  3.3235066789038505  -2.530193123599033  1.0987444462543554 
C  3.109992056469837  -3.714238267071516  0.4399591593618102 
C  2.1188600150310313  -3.8143784196789188  -0.47695466225347416 
C  1.28855328325815  -2.6970794567292957  -0.6564968817015736 
H  0.1460935615523351  4.063739715177946  -0.17612615085786182 
H  0.15865730269967027  5.806745349045468  -1.6827206028304162 
H  1.4415375107505355  5.701580014813371  -3.584217380179064 
H  2.8988143802341027  3.907815954370033  -3.8996275910632723 
H  2.9107730973281525  2.1688177919926623  -2.362826530309592 
H  0.8459757544610488  2.5139207302414155  2.1337214739822876 
H  3.6730284036277947  0.5680377948928761  -0.46549682018979777 
H  3.701282487078098  -0.21143062647251765  1.9785179458995126 
H  2.20173569534062  0.6873527892889222  3.2767568495463184 
H  1.0278147960677733  0.10314644444701959  2.401816412203919 
H  2.303636530065077  4.145075097027284  3.1074454723416074 
H  3.417026897263322  4.704950516493348  2.1568808533346338 
H  3.7996385836457156  3.723739858377562  3.3187217898847594 
H  5.622651157115892  1.8444470506520132  0.08074634992218321 
H  5.741865103629238  2.3763883037746805  1.5504367044469118 
H  5.362852918218599  3.3593024108358573  0.3890772805602349 
H  4.021192142854372  -2.4640015792958847  1.7104346479536372 
H  3.6489399027153633  -4.449379445677854  0.6225341481871953 
H  1.9922623458727624  -4.594063148314865  -0.9684649714999725 
H  0.5771821634332319  -2.7642947990949445  -1.251332865506968 
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