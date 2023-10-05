%chk=EKEVOLOL_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5842275901839349  1.4079499076671722  -2.7201806104765302e-15 
N  1.5842275901839356  -1.4079499076671718  1.9984014443252818e-15 
N  2.9579048951463425  0.4567420431636664  0.15797085281134507 
H  3.7471890969928605  0.8286028118947686  0.27822191265787355 
C  1.8623293934583176  2.465167775926635  -1.4485030673245058 
C  0.7933893249913695  2.9176792832824825  -2.2200315539900455 
H  -0.08189976180137348  2.583033893958463  -2.06229082384814 
C  1.0088764880685324  3.8626050421856144  -3.2255511125816723 
H  0.2815179773127876  4.180479589825982  -3.746258889356387 
C  2.2948367341070344  4.334185335194841  -3.458043505691477 
H  2.439872885256818  4.983717118974879  -4.134922539728729 
C  3.3698425786798016  3.8692386162560695  -2.7140283225920427 
H  4.246440775058163  4.18853355533871  -2.89073744511391 
C  3.1563797872050716  2.932044300737906  -1.7064047909016358 
H  3.88868925884701  2.6101585308450677  -1.1944043720697748 
C  1.7832262493741637  2.6156912773991774  1.339336868557866 
C  1.2053362371075376  3.8869008532892115  1.2482041747100343 
H  0.7085472093169706  4.1296152873319745  0.4744408094786342 
C  1.3586695563012445  4.792234241258949  2.2882709657635045 
H  0.9710380718418887  5.656647128658137  2.2221773779199374 
C  2.0787030459600553  4.4432971645031385  3.4274594401762712 
H  2.1751213898539445  5.064083294918948  4.140201738096189 
C  2.653288771546021  3.184400959923824  3.517298852733603 
H  3.1523373297809005  2.945801260338908  4.290683427108928 
C  2.503158885691765  2.2709349814850315  2.481125516072071 
H  2.892929644960282  1.4074258513808104  2.551502411754548 
C  2.835277484972959  -0.9066194724320742  0.09693196302251372 
C  3.967583283785819  -1.7341274637713642  0.1438543696571295 
H  4.836943131742221  -1.3606060926372432  0.2298551258591092 
C  3.7940237384198916  -3.103736080170976  0.06283645746291996 
H  4.5402458705721624  -3.68853018791322  0.12181137616616408 
C  2.5022385890219705  -3.616062325637455  -0.10816642555726559 
H  2.3635082329245845  -4.551703798268176  -0.19835252753816984 
C  1.4366815097809535  -2.7485409434145778  -0.1441454787061078 
H  0.5642476062180286  -3.101035646739861  -0.2743309542461618 
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