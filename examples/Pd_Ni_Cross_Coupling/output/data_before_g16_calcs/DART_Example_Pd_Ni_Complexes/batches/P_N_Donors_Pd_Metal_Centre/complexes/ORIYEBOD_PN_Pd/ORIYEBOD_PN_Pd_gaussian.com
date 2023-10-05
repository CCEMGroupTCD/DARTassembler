%chk=ORIYEBOD_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.3138773173321785  1.6630773869546776  -4.960002975792329e-15 
N  1.3138773173321785  -1.663077386954678  -8.881784197001252e-16 
C  0.7702194155227764  -2.7380046417070867  0.8699289124926664 
H  1.4823709448632918  -3.3768519336815714  1.0846695858520496 
H  0.04545213897080913  -3.202846573335772  0.402371308085283 
H  0.42284637649210555  -2.3452232834107383  1.698613554759847 
C  1.8053328748766695  -2.3448087867640948  -1.2206032580132535 
H  2.4665439960586744  -3.0256788796907195  -0.9708484586933193 
H  2.2215880591643087  -1.6892773012471731  -1.8177208871698993 
H  1.053048976498066  -2.7738902658326605  -1.6805460376724553 
C  2.496357715058177  -1.082211963867588  0.7132003289789437 
H  3.042731054030134  -1.8245248800816571  1.0701217893626003 
H  2.1691264739402225  -0.5536888765508146  1.4837652587607002 
C  3.3745726554707542  -0.2019732530998286  -0.13472020415023578 
C  4.644836369872662  -0.6396413710776934  -0.4781455444785825 
H  4.9347761231577785  -1.4968369088093665  -0.1908634365117331 
C  5.494218718063422  0.14688989334845615  -1.2270401409754763 
H  6.365689728354766  -0.1658037606445506  -1.4374264926580262 
C  5.081123475326914  1.3854732185489067  -1.6705337113503922 
H  5.658163417190504  1.9199698342642542  -2.2039262903041554 
C  3.807932040260909  1.8472882735935636  -1.3332209712874985 
H  3.524460841813268  2.705073438126545  -1.628912011323862 
C  2.956761079682577  1.071413870440991  -0.5723826662575568 
C  1.7385481149366129  2.185570790939419  1.7231893820629607 
C  3.05755117667531  2.2453793246734617  2.174665065147394 
H  3.759695871554051  1.9431603576292678  1.6117285470457619 
C  3.3565314129445754  2.739640753184123  3.4293142585903342 
H  4.258397380497708  2.7457579880643537  3.7252146573851705 
C  2.372124469024284  3.2188227451052946  4.250546751481387 
H  2.5812753432589792  3.582445466298894  5.102778109997726 
C  1.046805543983473  3.15973291214211  3.802221201719676 
H  0.35140687214106114  3.488286423413497  4.358989822810268 
C  0.7345361835201332  2.6380401083040863  2.5767843143446187 
H  -0.17520099866596928  2.583893815590136  2.305608267589993 
C  1.1017890155742978  3.3065510559731135  -0.8143898215845702 
C  1.3110229686652066  4.515602845007388  -0.16614214954325943 
H  1.6856245709325932  4.517276774929012  0.7055509429711574 
C  0.9814656588222669  5.7218929074761835  -0.7709548103880226 
H  1.1220868925442447  6.537775367693278  -0.30286061074802073 
C  0.4475041735026357  5.743917943885988  -2.053931353418501 
H  0.20951562143627367  6.567437213490081  -2.462881376864038 
C  0.2695085997708633  4.543939677297734  -2.7289187869665352 
H  -0.07654476502960206  4.547432199681575  -3.614708542494152 
C  0.5911280339842102  3.3429568543291777  -2.1206487659472737 
H  0.4635102265333144  2.52955749263317  -2.5962685969194017 
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


