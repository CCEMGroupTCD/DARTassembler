%chk=VENIVIGU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.523495315056794  1.4734524169446392  0.0 
N  1.523495315056794  -1.4734524169446392  0.0 
C  1.8768067863105162  2.741669112312963  1.2302654786141172 
C  1.5144124884175074  2.51678898982683  2.534533960701025 
H  1.0783820425335475  1.7362693057336243  2.7863199273577903 
C  1.8302831957133812  3.5210512469470046  3.4799765053309564 
H  1.561504713641334  3.4165047490564056  4.3641467872076385 
C  2.5012613551241785  4.604332978100437  3.1292074352250854 
H  2.7342129409844027  5.227417927734406  3.7806811217119853 
C  2.8476303365560662  4.823948831506145  1.8708608074492783 
H  3.2939908418929917  5.606734411641602  1.64194087645294 
C  2.546090264319803  3.8807520549982173  0.8863566451315983 
H  2.799659011687329  4.026486445556461  0.004627974980112581 
C  1.2668080638029477  2.3164947647768135  -1.5800107324969075 
C  0.41906193242514145  3.3602067349189615  -1.661154001570634 
H  0.021113639200353695  3.680866940546594  -0.8830516185603108 
C  0.12252715851614138  3.9713831730881455  -2.8668759219526203 
H  -0.4687890081967956  4.688090735446895  -2.896698917135714 
C  0.6884606814442318  3.5275822756312922  -3.980110879651459 
H  0.49325268802311006  3.931275586768837  -4.79512521685501 
C  1.5626431958567053  2.4721723980428076  -3.9280703014055938 
H  1.9674098651683138  2.161637466790127  -4.707316610274815 
C  1.840604912283971  1.8738131705163412  -2.719516306321051 
H  2.4303749087876327  1.1540914490650236  -2.6875267002710457 
C  3.0120109615610096  0.48030306364931596  -0.1772406899142055 
C  4.277980394414152  1.057904206524889  -0.33596971949096566 
H  4.367695686876289  1.9774930642429251  -0.44260066743419957 
C  5.392564375739839  0.24032525986436593  -0.33019871007720125 
H  6.238674979138956  0.6093331244871352  -0.44317265690668783 
C  5.257286437017925  -1.079348665563154  -0.1613035359546826 
H  6.020658830427095  -1.609810247499238  -0.12992824200067465 
C  4.016764469223214  -1.6764098396350229  -0.034296639155227515 
H  3.948722749356793  -2.5987114283201627  0.05879446604010652 
C  2.8700423304140292  -0.8835299057023519  -0.048157617903719085 
C  1.0962894281334612  -1.8920390797422963  -1.3940896620755807 
H  1.2763494835072726  -1.1697375464448112  -2.0158971093250444 
H  0.14149546199237562  -2.062267410237087  -1.4009894966787 
C  1.8371338459582214  -3.1485001664456873  -1.8492392486710654 
H  1.4197510478667061  -3.4781950038312996  -2.660431312277826 
H  2.7518779261450494  -2.9078490087141913  -2.067602908755937 
C  1.8573673883539807  -4.1352416881127745  -0.9434888894769602 
H  2.785837913530899  -4.333401589065294  -0.743740053354705 
H  1.4801763480590449  -4.925290197051992  -1.3608658445716937 
C  1.2220541519987782  -3.9472728356307143  0.2235899177200962 
H  1.4903406698354298  -4.657664612088962  0.8257851753297775 
H  0.27078959148289283  -4.050884043732411  0.06411153842607342 
C  1.4242455218938335  -2.658335769757786  0.9255219871960897 
H  2.2379646904209722  -2.7121749600076472  1.450964035727382 
H  0.6864596129179601  -2.515921313351295  1.5381872588764054 
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
-N 0
6-31+g(d)
****
-C 0
6-31g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****

-Pd 0
lanl2dz

