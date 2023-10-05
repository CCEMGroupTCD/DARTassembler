%chk=AYUKEGUQ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5872363056583605  -1.404557193566712  2.7322208351371497e-15 
N  1.5872363056583607  1.404557193566712  -5.475193438312842e-18 
N  2.756574860355988  1.2693458964619357  0.7142650135859557 
C  3.7335684021995506  2.0510842585800098  0.1772522089071956 
C  5.134859223390258  2.0803802674491894  0.6915142416405751 
H  5.152184398949148  2.5105572675510377  1.5489641275670043 
H  5.690834483929907  2.566112627174272  0.07832129942435924 
H  5.462842580906879  1.1820808847120419  0.7781638588710512 
C  1.8373124617155332  2.273417540151122  -0.982804621633878 
C  3.1693190551369606  2.701754551100634  -0.8943245101244365 
H  3.5919570930926574  3.311740711809657  -1.4548457520202438 
C  0.8056107860354337  2.6439401877041573  -1.9827553154429294 
H  0.08940880749843916  2.00549211018384  -1.9556024472644018 
H  1.1968449244139587  2.648573605289292  -2.8597020643443254 
H  0.46286644716733827  3.518022201326116  -1.7814709602723056 
C  2.9077865909363845  0.2925812077160401  1.749063721235921 
C  2.5545876298347903  -1.0492321400296196  1.496637090928538 
C  2.8153294103850746  -2.000703576095797  2.47621316261755 
H  2.604784355312187  -2.8926671114635685  2.319848488350649 
C  3.3877191202415564  -1.6283575994453554  3.686076350636179 
H  3.585490129134096  -2.2744795901743124  4.325701046171054 
C  3.664628029986414  -0.30203009314469037  3.941450315116085 
H  4.01392308157945  -0.055021540136613506  4.767787801663472 
C  3.426752508469196  0.6746941771687913  2.973437052787781 
H  3.614205139592618  1.5690003796882257  3.1477928878497887 
C  2.6690155950946224  -1.1015564528711452  -1.4102544132477592 
C  4.0384943379823985  -1.4109724538645987  -1.3465946841579262 
H  4.407153182564102  -1.7396386926723908  -0.5578144550829074 
C  4.829568333091294  -1.2278474724269632  -2.4548130495697054 
H  5.737677563544154  -1.4252695220736944  -2.4090451582085204 
C  4.292920427120475  -0.7531040778153666  -3.631559617140579 
H  4.839051157767822  -0.6375745637854401  -4.374616442982701 
C  2.9315016934193654  -0.44312958051888507  -3.714841383173336 
H  2.5702205037113126  -0.12412062022002568  -4.509324872751453 
C  2.1267226466023104  -0.6171664735123698  -2.5965265725402005 
H  1.2209890015099236  -0.40977545281616423  -2.640258239353384 
C  1.2160531819745384  -3.177687388198198  -0.02239400764553673 
C  1.6783333818800306  -3.9753044344326267  -1.0807914397718774 
H  2.212800412620772  -3.6081324514317754  -1.7478096838630088 
C  1.3248068900117322  -5.326761739762954  -1.1211695454305317 
H  1.6425605812061645  -5.865966254515973  -1.8099377960070862 
C  0.51192515932649  -5.868376565762468  -0.1519039202401408 
H  0.2670842457177862  -6.7647413530182465  -0.1986644291249178 
C  0.05700308737525206  -5.078510061869684  0.894851903474162 
H  -0.48397728213948277  -5.450882822762285  1.5541644141585385 
C  0.40412913968726993  -3.7383486582460788  0.9658962279756778 
H  0.09661823327759156  -3.2144890437212035  1.669397753708888 
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


