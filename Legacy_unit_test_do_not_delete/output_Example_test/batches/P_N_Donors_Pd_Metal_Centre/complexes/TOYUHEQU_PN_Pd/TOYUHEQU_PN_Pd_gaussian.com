%chk=TOYUHEQU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.540387987651163  -1.4557832419354195  2.1272103130662645e-14 
N  1.5403879876511635  1.4557832419354195  -2.1161497194164316e-14 
C  2.5223477449493523  1.4260961609503928  0.8279473426912348 
H  3.1012765726179286  2.1538874793828127  0.8517704262173408 
C  3.147164331691868  -1.90949010293887  2.4453496671253117 
H  3.6152336224142543  -2.6375128809665793  2.0069304228448255 
H  2.476267842645001  -2.2798092544928235  3.0399104604651974 
C  2.522042227932788  -0.9870905476213159  1.4379438969701914 
C  0.6457011258341029  -5.890716354928619  0.8832669995921133 
H  0.4921694386030955  -6.784761105210018  1.0891466137122872 
C  2.6564332859245265  -1.4050305351644965  -1.4294373560801314 
C  4.815743788964142  -1.009229526881189  -2.424535061486436 
H  5.705227224472275  -0.7518577154692915  -2.3410035785160983 
C  2.147201188757418  -1.7231850011990588  -2.694088955121174 
H  1.250525392160796  -1.9473354520355297  -2.789243309659691 
C  1.592842729454723  2.5269263524488133  -1.0899201515992247 
C  1.6664083893204702  -4.209266243066045  -0.5308099479547959 
H  2.193168054717683  -3.9919547739634096  -1.2662687369478198 
C  1.1463573475637006  -3.2074892113005484  0.28404985753007045 
C  1.3017821502853326  3.8923668162072165  -0.4467311831110982 
H  0.43716790085817436  3.8715463745872865  -0.0304403846320904 
H  1.313735062124856  4.573132383808076  -1.1245508788028447 
H  1.9700100638325264  4.0847240620755825  0.21318774790870287 
C  2.7782242561254344  0.30678457797486436  1.7352179755484798 
C  3.6263583870609017  0.4394813555138097  2.985742038161107 
H  3.097346330147393  0.7510213375217448  3.737828797285778 
H  4.369709641112333  1.0457313037043448  2.845072626350492 
C  2.979007951387104  -1.705602328681172  -3.806035707115495 
H  2.6385458874768117  -1.927713828996383  -4.641772390805798 
C  4.315695730387754  -1.358494520702695  -3.67334625346258 
H  4.8745910271193775  -1.3597616770669636  -4.418256269240334 
C  2.9883863305335  2.5093738836032173  -1.7608860436723932 
H  3.154943574720677  1.6394225999792293  -2.129538927195301 
H  3.6589828503385795  2.714612727816584  -1.106486544550268 
H  3.0136474547875904  3.1647118917064234  -2.4618975633614353 
C  0.35814970526073675  -3.543078051676588  1.3963030583932161 
H  -0.00538778749290203  -2.874066613747687  1.9291114172641333 
C  0.12540145434357264  -4.88414667972623  1.696398004295424 
H  -0.37970341481272163  -5.108175625523406  2.4448146212940305 
C  1.3891065527967505  -5.554417340145967  -0.22808223616537102 
H  1.7109194118904383  -6.226013067336912  -0.7854056188982497 
C  3.995564124150055  -1.0402345891045808  -1.2978842266924249 
H  4.339234471290132  -0.8181877137948567  -0.46296163287799075 
C  4.121856128072415  -1.0169478187780956  3.210457811803941 
H  5.023569225586236  -1.1271840375895255  2.870612073513417 
H  4.115738471170554  -1.2372982258713314  4.155430833835722 
C  0.5430492367081657  2.184417223591742  -2.1499753479152623 
H  0.7530023588485201  1.3355564114854652  -2.547220534385525 
H  0.5415005182208399  2.863919141531714  -2.8285755357924516 
H  -0.32369572150774806  2.1410258496698855  -1.739847592116565 
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
