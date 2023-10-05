%chk=IWITEFOS_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.523495315056794  -1.4734524169446392  1.8044587861071838e-16 
N  1.523495315056794  1.4734524169446392  -1.8044587861071838e-16 
C  1.8768067863105162  -2.741669112312963  -1.2302654786141167 
C  1.5144124884175074  -2.5167889898268303  -2.5345339607010247 
H  1.0783820425335475  -1.7362693057336247  -2.7863199273577903 
C  1.8302831957133812  -3.521051246947005  -3.479976505330956 
H  1.561504713641334  -3.416504749056406  -4.3641467872076385 
C  2.5012613551241785  -4.604332978100437  -3.129207435225085 
H  2.7342129409844027  -5.227417927734407  -3.780681121711985 
C  2.8476303365560662  -4.823948831506145  -1.8708608074492776 
H  3.2939908418929917  -5.606734411641602  -1.6419408764529393 
C  2.546090264319803  -3.8807520549982173  -0.8863566451315978 
H  2.799659011687329  -4.026486445556461  -0.004627974980112088 
C  1.2668080638029477  -2.3164947647768135  1.5800107324969077 
C  0.41906193242514145  -3.3602067349189615  1.6611540015706345 
H  0.021113639200353695  -3.680866940546594  0.8830516185603112 
C  0.12252715851614138  -3.971383173088145  2.8668759219526208 
H  -0.4687890081967956  -4.688090735446895  2.8966989171357143 
C  0.6884606814442318  -3.527582275631292  3.9801108796514595 
H  0.49325268802311006  -3.9312755867688365  4.795125216855011 
C  1.5626431958567053  -2.472172398042807  3.928070301405594 
H  1.9674098651683138  -2.1616374667901264  4.707316610274815 
C  1.840604912283971  -1.873813170516341  2.7195163063210512 
H  2.4303749087876327  -1.1540914490650234  2.6875267002710457 
C  3.0120109615610096  -0.48030306364931596  0.17724068991420555 
C  4.277980394414152  -1.057904206524889  0.33596971949096577 
H  4.367695686876289  -1.9774930642429251  0.4426006674341998 
C  5.392564375739839  -0.2403252598643659  0.3301987100772013 
H  6.238674979138956  -0.6093331244871352  0.4431726569066879 
C  5.257286437017925  1.079348665563154  0.16130353595468247 
H  6.020658830427095  1.609810247499238  0.12992824200067446 
C  4.016764469223214  1.6764098396350229  0.03429663915522731 
H  3.948722749356793  2.5987114283201627  -0.05879446604010684 
C  2.8700423304140292  0.8835299057023519  0.048157617903718974 
C  1.0962894281334612  1.8920390797422966  1.3940896620755805 
H  1.2763494835072726  1.1697375464448114  2.0158971093250444 
H  0.14149546199237562  2.062267410237087  1.4009894966786998 
C  1.8371338459582214  3.1485001664456878  1.849239248671065 
H  1.4197510478667061  3.4781950038313  2.6604313122778254 
H  2.7518779261450494  2.9078490087141917  2.0676029087559367 
C  1.8573673883539807  4.1352416881127745  0.9434888894769596 
H  2.785837913530899  4.333401589065294  0.7437400533547045 
H  1.4801763480590449  4.925290197051992  1.360865844571693 
C  1.2220541519987782  3.9472728356307143  -0.22358991772009668 
H  1.4903406698354298  4.657664612088962  -0.825785175329778 
H  0.27078959148289283  4.050884043732411  -0.06411153842607392 
C  1.4242455218938335  2.658335769757786  -0.92552198719609 
H  2.2379646904209722  2.7121749600076472  -1.4509640357273823 
H  0.6864596129179601  2.515921313351295  -1.5381872588764056 
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
