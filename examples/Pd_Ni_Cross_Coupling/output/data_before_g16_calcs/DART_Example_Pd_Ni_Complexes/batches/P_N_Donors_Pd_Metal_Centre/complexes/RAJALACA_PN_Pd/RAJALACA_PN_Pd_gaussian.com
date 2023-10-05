%chk=RAJALACA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.521536413465021  -1.4754751582117547  1.8528561213240758e-16 
N  1.521536413465021  1.4754751582117542  -2.0844916858817484e-16 
C  2.7605029050923715  1.4442655341501665  -0.31656160852289245 
H  3.2790398386073503  2.179091724862893  -0.0800585894649763 
C  3.435520092672488  0.3406097077968381  -1.024436359752094 
C  3.003352850933746  -1.0032742125554337  -0.9639506276614104 
C  3.768643719421144  -1.9882390131565912  -1.587006554674831 
H  3.4965435260354507  -2.877112511094384  -1.5437578619396628 
C  4.938915523046788  -1.658323550133838  -2.276770161207759 
H  5.4249954755710394  -2.3215972247554406  -2.7115595381570365 
C  5.376356067835097  -0.34379681669110235  -2.31450511635112 
H  6.158569545008383  -0.12380704716262281  -2.7651410641307677 
C  4.637284492607517  0.6430721061071405  -1.6735580468472189 
H  4.946912119459611  1.5204553800952836  -1.675923044898431 
C  1.084863311336194  2.539392941756512  0.9562637302200812 
H  0.10437534891530631  2.555130539825017  0.9458342741082585 
C  1.5687696643247997  3.95572740077428  0.616288993311406 
H  1.2784427174833404  4.198606813555671  -0.276645861414683 
H  2.53778767121359  3.9876052502509234  0.6400428359069217 
C  0.9870092383254598  4.939586604991168  1.640594125105063 
H  1.2849639174165643  5.838664198820515  1.431182975731985 
H  0.018179443906701787  4.9247726178641384  1.5895688839795918 
C  1.425979957649839  4.575705354932909  3.0550137104586805 
H  1.027404490432784  5.1951295247763944  3.6861728304485815 
H  2.38974567685843  4.654469368242932  3.1248447829526462 
C  1.00476408835983  3.1468338464792023  3.407355566881526 
H  1.3537319153829  2.915754483263608  4.282710148521825 
H  0.03652606698536798  3.0980961299546115  3.448679881187729 
C  1.5204156054793785  2.1423330036481403  2.3741400522651435 
H  1.1768102856891183  1.259819480629459  2.5825360356740927 
H  2.489538454507501  2.105140784792684  2.4156864401739293 
C  1.244663648301877  -3.203641739538946  -0.48090117393786685 
C  0.6697357028552406  -3.441860206345253  -1.7343794544106523 
H  0.36497149902988957  -2.73032335385605  -2.2512525990181698 
C  0.5552835308616066  -4.747571650741874  -2.2066193238600964 
H  0.18342325057754705  -4.906998746949373  -3.043902778218263 
C  0.9963207997648516  -5.812175226923918  -1.4280893995676627 
H  0.9291454336158471  -6.681205713722518  -1.7512917405204207 
C  1.5360369183010671  -5.582261603429814  -0.17119251949991104 
H  1.8122788689091267  -6.298028313646724  0.3555563212862134 
C  1.6646540531324947  -4.276096217130971  0.3039035300127691 
H  2.030990753995641  -4.121397967059015  1.1440397172023171 
C  2.0693719247020534  -1.488230405153905  1.7333117568121381 
C  3.3722170764700303  -1.162987680917454  2.1145081508742876 
H  4.017697484124563  -0.9983619780392081  1.46517356740749 
C  3.7091758483966877  -1.08427130649827  3.4641985666437742 
H  4.575210086122635  -0.8544168503845961  3.7140765930708595 
C  2.750990136482332  -1.349267852067956  4.436226931620487 
H  2.9785144991161028  -1.3020886470419832  5.337573899786044 
C  1.4528203389167862  -1.6843769403153366  4.0649992849363485 
H  0.8169051256060831  -1.868960400419185  4.717841714421496 
C  1.10584508439018  -1.7432133290656657  2.7229059900775106 
H  0.23349224638869703  -1.9519324649223766  2.477130831939015 
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


