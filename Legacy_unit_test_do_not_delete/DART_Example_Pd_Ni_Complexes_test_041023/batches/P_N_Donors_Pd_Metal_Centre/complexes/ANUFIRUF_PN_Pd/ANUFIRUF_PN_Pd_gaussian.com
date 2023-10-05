%chk=ANUFIRUF_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.4855741154516657  1.5116777260712682  -1.37704196674127e-15 
N  1.4855741154516657  -1.5116777260712682  0.0 
N  3.425874038622977  -3.457928553386359  0.5080516069424422 
N  3.1209970966993366  -2.528458610310834  1.4062881766736486 
C  3.5892754518766514  -0.32298456481787285  0.3235908612944262 
C  3.204986851857543  0.9185634633307651  -0.20731587787508274 
C  4.144958132407078  1.706298279067449  -0.8661772455529855 
C  5.452104222398125  1.293600182614457  -0.9586177161557234 
C  5.832933939072718  0.07129719889444086  -0.42392225808011075 
C  4.908952673067381  -0.7233566577327792  0.19799122321668067 
C  2.55095662698902  -1.2373337783848632  0.9675994587008876 
C  1.745116289047249  -2.4445964796373394  -0.8973854018035013 
C  2.853188252214091  -3.355287039706554  -0.6789601325592306 
C  3.210940381131778  -4.340225253931327  -1.6974963510682217 
C  4.294552537362481  -5.193442978685441  -1.5429794619121 
C  4.598359175921783  -6.1229640636139315  -2.56772259999144 
C  3.8408552736853  -6.173181910110415  -3.697635652394045 
C  2.795418732468381  -5.317545948362428  -3.8586228234730577 
C  2.4457878927033585  -4.386729190443443  -2.856720084150292 
C  1.2941162369007446  -3.528132478581322  -3.005187095975988 
C  0.9318585726613708  -2.6173051367123668  -2.0676827770505817 
C  3.6360480080848805  -2.69604215714154  2.7179093199587867 
C  3.8012147983786235  -1.6018356265612135  3.5638283348883726 
C  4.311751881112897  -1.8022788955028601  4.824037378297136 
C  4.649681596755134  -3.052966857885966  5.260719251101506 
C  4.495739192898283  -4.138682147391146  4.405142494504265 
C  3.9944205557018027  -3.9752599174148924  3.1184963004209525 
C  1.3175277691477658  2.825795024844211  -1.2211933960514372 
C  1.1002443637263064  2.4486679856824773  -2.543638120061331 
C  1.1320853221732707  3.462482959977963  -3.5518549227676943 
C  1.2776389033861542  4.762801133685928  -3.219066488291838 
C  1.4452189142985243  5.122368361955305  -1.9504681719216357 
C  1.4734699638646827  4.154497184940276  -0.9407255696204433 
C  1.5287001637400541  2.374621836116652  1.60063609446278 
C  0.37401286637532993  2.5939896446652617  2.290214335039718 
C  0.40996674242060127  3.37261742339097  3.4836689512669206 
C  1.6167648579515896  3.8871732785080186  3.9224409083421468 
C  2.7504430721429145  3.646409021204382  3.2440188501268525 
C  2.7138659134961833  2.8862012532519805  2.08615151704841 
H  3.8816353536462778  2.535338982245022  -1.2504174004828434 
H  6.0937314561191185  1.8433980004072787  -1.3925268339417038 
H  6.737956162645376  -0.21401061345316652  -0.48788112809105166 
H  5.178544053988674  -1.5647281505214723  0.5495340349131882 
H  2.1770063922773866  -0.7959181555363322  1.7212833711365785 
H  4.826488502803765  -5.154712884968864  -0.7563937714048459 
H  5.334966276735427  -6.715045202711037  -2.466512700493366 
H  4.044154498098822  -6.806692104435292  -4.374392077586564 
H  2.290510676663616  -5.34606983560132  -4.66152339828928 
H  0.7701519013726292  -3.6038302268688467  -3.7932912172688615 
H  0.14457348888491617  -2.0978078447591635  -2.185061695585434 
H  3.564446283887874  -0.7280505087286795  3.2743008950830523 
H  4.430894185643694  -1.0587299659825424  5.403509543875665 
H  4.988498906486392  -3.1799850256398954  6.140713231584879 
H  4.738504631013484  -5.006682609639234  4.704667167267363 
H  3.902446554943877  -4.716563748696698  2.5341504499937337 
H  0.9362981751822397  1.54081030124948  -2.7724197453677184 
H  1.0488223095081417  3.2184278722848823  -4.4669494779329115 
H  1.2608064600118887  5.42645693987904  -3.899216427464236 
H  1.5450560560930309  6.0398280930082  -1.7290683367485884 
H  1.6034893142910835  4.426746274744353  -0.04053581650415507 
H  -0.44777155530068935  2.23100893089497  1.9759190507584574 
H  -0.38785844490413224  3.5388150531284572  3.9723552395336874 
H  1.6415638244909683  4.411567834522419  4.713850954324022 
H  3.573377403936674  3.993292590338402  3.5601056815830656 
H  3.5235122162170995  2.7164078604055635  1.618902152965288 
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


