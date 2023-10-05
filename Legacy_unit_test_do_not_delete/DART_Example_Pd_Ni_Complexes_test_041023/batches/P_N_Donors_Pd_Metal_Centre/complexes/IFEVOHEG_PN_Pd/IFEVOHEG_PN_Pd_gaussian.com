%chk=IFEVOHEG_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5403879876511688  -1.455783241935419  2.01671513764636e-14 
N  1.5403879876511635  1.4557832419354195  -2.1161497194164316e-14 
C  2.492854616040767  1.4260961609503935  0.861712880852273 
H  3.070599362284704  2.1538874793828136  0.9057257767117771 
C  3.060844054312417  -1.9094901029388687  2.4999357118997567 
H  3.543928820839782  -2.6375128809665775  2.078118923445283 
H  2.3696063848330606  -2.2798092544928226  3.070720365061615 
C  2.471260712402332  -0.9870905476213154  1.471327179321812 
C  0.6154205711569993  -5.890716354928619  0.8515048160172821 
H  0.45479731636747056  -6.784761105210018  1.0519008352269688 
C  2.7056400655783093  -1.4050305351644965  -1.389617162257981 
C  4.898363584090578  -1.0092295268811868  -2.3087498312482744 
H  5.784389963608359  -0.7518577154692896  -2.19422670917598 
C  2.2408538831563427  -1.723185001199058  -2.6712703136105174 
H  1.3480451563100897  -1.9473354520355288  -2.797660236693768 
C  1.6308484401800196  2.526926352448813  -1.0874255576029268 
C  1.6848566211250344  -4.209266243066044  -0.5260885442840043 
H  2.236962540347946  -3.991954773963408  -1.2427156644486579 
C  1.1366841833170767  -3.2074892113005475  0.27012535100721674 
C  1.3175181959664302  3.892366816207216  -0.4547862701789669 
H  0.43890230683013165  3.8715463745872856  -0.06892366695341284 
H  1.3531193926497043  4.573132383808076  -1.1317759058188435 
H  1.9623082045665639  4.084724062075582  0.22805147392383096 
C  2.716911965700121  0.3067845779748648  1.7773607904116515 
C  3.520886775839292  0.43948135551381035  3.0567225218592267 
H  2.965949529407356  0.7510213375217454  3.7898888755150932 
H  4.268694492065052  1.0457313037043454  2.9420813866930904 
C  3.1109603135854944  -1.7056023286811715  -3.7535100603178124 
H  2.7998724396021637  -1.9277138289963822  -4.600619590478011 
C  4.442203023357045  -1.3584945207026933  -3.5742517067593953 
H  5.026754840814042  -1.3597616770669627  -4.299202778917013 
C  3.0489582857412456  2.5093738836032182  -1.7092789460794833 
H  3.228279867848361  1.6394225999792298  -2.0718944922321305 
H  3.6963080831036432  2.7146127278165837  -1.031874608419295 
H  3.0986689708207793  3.1647118917064234  -2.4089818279688324 
C  0.3101396188972716  -3.543078051676588  1.3541929472590346 
H  -0.07177116020171437  -2.874066613747687  1.8739894581441283 
C  0.06705998935096891  -4.88414667972623  1.6459822866130989 
H  -0.46385654682354627  -5.108175625523406  2.376315082714627 
C  1.3971586646083347  -5.554417340145966  -0.2332229405665234 
H  1.7382257895612658  -6.226013067336911  -0.7789757101382574 
C  4.039364003474021  -1.0402345891045794  -1.2114091792093757 
H  4.3536866176146205  -0.8181877137948552  -0.36500127553347916 
C  4.00824020561363  -1.0169478187780947  3.298594026502892 
H  4.9212644490947515  -1.1271840375895237  2.9904246463271087 
H  3.969147192556711  -1.2372982258713305  4.24277789335511 
C  0.6186898460930799  2.1844172235917414  -2.18347226147351 
H  0.8423787275394735  1.335556411485465  -2.5731481986155735 
H  0.6408248760623494  2.863919141531713  -2.8617131139471597 
H  -0.2618403667740572  2.1410258496698846  -1.803843307233538 
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


