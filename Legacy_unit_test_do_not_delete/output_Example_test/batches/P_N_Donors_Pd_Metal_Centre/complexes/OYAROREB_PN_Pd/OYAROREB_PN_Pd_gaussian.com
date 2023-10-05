%chk=OYAROREB_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.4278016922878298  1.566359577970526  -1.1508594778551925e-15 
N  1.4278016922878327  -1.566359577970525  -3.3861802251067274e-15 
C  1.769398388126053  -2.5204897094032974  -0.7679781570360666 
H  2.4159058353947263  -3.13036898516815  -0.43435951304461634 
C  1.2543126093547252  -2.7714419846843645  -2.136979974433982 
H  0.5568984015767591  -2.09061712812  -2.3586514599109387 
C  0.6277264500395998  -4.179909932292189  -2.1821437377420096 
H  0.28482573070136286  -4.355146773730368  -3.083767027192814 
H  1.3080011794458033  -4.847695346772338  -1.9562318240701135 
H  -0.10798343544252487  -4.2305989861446065  -1.5374769290544104 
C  2.396065216899419  -2.6735337315879066  -3.147162124910085 
H  2.0494969152202263  -2.840841684809953  -4.04929685490326 
H  2.7899759518432887  -1.777027444175002  -3.110641589247074 
H  3.081437592755891  -3.3405873473060628  -2.932813016747234 
C  2.1193796451390607  -1.4762380490538414  1.295838687937612 
H  2.6173822485667193  -2.316265679685639  1.4597859256527062 
H  1.4500964145543447  -1.3720520452171985  2.0178004901561524 
C  3.0950105077504952  -0.29075458400815335  1.329526579283736 
H  3.7163012682823116  -0.412208605597814  2.0899856850333354 
H  3.6322241025616826  -0.29759350538397134  0.49816867638963985 
C  2.4177975976075827  1.0659951207776985  1.4618799052686724 
H  3.113125524738974  1.7514982262595256  1.6246653349788756 
H  1.825170193473114  1.0470600209340035  2.255374903276069 
C  0.8650430282988124  3.2460865493527713  0.4513547692407142 
C  0.5649786717243404  3.575651904115591  1.7701543896734924 
H  0.6679827793063312  2.9198544043289028  2.4497815544445523 
C  0.11920044730584367  4.842046720413212  2.107967420637082 
H  -0.07780057920657768  5.052325815732279  3.013032405941914 
C  -0.03930951026093665  5.806412311331458  1.118674107000988 
H  -0.3429159569771567  6.678104513558154  1.3439235023154454 
C  0.2460705837314876  5.4863136398312475  -0.19575351841500122 
H  0.13394379507496668  6.143533140690567  -0.8732629950688762 
C  0.6918326403328313  4.228665474170355  -0.5371436768911066 
H  0.8837380882988594  4.024721740884096  -1.4455114287924158 
C  2.716317476865671  1.842488399910886  -1.2556555883631777 
C  3.7375172926660793  2.7629828223142763  -1.0228563518235 
H  3.724355200108591  3.301070906847226  -0.24078434580156843 
C  4.78074387530159  2.891295786225416  -1.944176046921258 
H  5.471090702447886  3.5262205978632153  -1.7939649893530611 
C  4.808291614543633  2.0954981401861534  -3.073934312490072 
H  5.529614949033928  2.1722371161470493  -3.687041673959729 
C  3.787499832048857  1.1846615129740807  -3.3179585435002887 
H  3.805231903231319  0.6465043798002434  -4.100007578833572 
C  2.742106693011449  1.0658997212269878  -2.4109783443368387 
H  2.039387933257904  0.449633170062091  -2.580440235645084 
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
