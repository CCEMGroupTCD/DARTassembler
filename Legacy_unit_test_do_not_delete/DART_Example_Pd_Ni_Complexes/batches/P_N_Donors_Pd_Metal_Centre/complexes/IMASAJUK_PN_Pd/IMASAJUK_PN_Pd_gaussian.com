%chk=IMASAJUK_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.3785481620168318  -1.609877313648465  3.2254755843738433e-15 
N  1.3785481620168376  1.609877313648465  4.689807048571407e-16 
N  -1.4059378368740898  -3.2960483199645063  -3.188983499399123 
C  1.7861390567724087  -2.543016585979394  1.5438181405029905 
C  1.730733462880792  -1.6473519053384837  2.8157446792692715 
C  2.010223505738997  -2.5496399664752873  4.060997995134409 
C  1.0241535217894264  -3.688439401349693  4.198369512864432 
C  1.0850317793541342  -4.570687851563051  2.973629356810069 
C  0.843387434452614  -3.758754803051614  1.699460217601281 
C  2.585064791133964  -0.3848158457466625  2.7554484124233456 
C  2.1786703423208524  0.6185409388391668  3.808156237415405 
C  4.093289318868323  -0.6628609962715558  2.8475150633348347 
C  0.053049493689484084  -5.712221693564434  3.072522782476055 
C  2.9809529069540437  -1.0012386014354173  -0.6789522232075299 
C  2.8246096876160003  0.3151913227785273  -1.4556965134350794 
C  2.5873918114129673  1.507059919678937  -0.5720816329965576 
C  3.580839910472807  2.425526465212774  -0.3096450326345811 
C  3.36313254769426  3.4847623371286374  0.5144733762355671 
C  2.1300518054344852  3.5736726744484506  1.1405051445874734 
C  1.1664124647979142  2.640408920582002  0.8661140242926678 
C  0.8987475355048751  -2.901029342674376  -1.1924264198436079 
C  0.7171229029507142  -2.3558779000398844  -2.649359392716121 
C  -0.0891581634524008  -3.3166013551054334  -3.460194329584773 
C  0.4195636346814481  -4.152604209850851  -4.386982395412681 
C  -0.4383004013043279  -5.081043142636742  -4.989837845624105 
C  -1.749619159618053  -5.077151516891797  -4.703438999006958 
C  -2.172085401283686  -4.187502636286017  -3.77123125719732 
H  2.674348354061671  -2.872890471708279  1.4598104996417955 
H  0.8307609092108794  -1.356547611724948  2.8973434259312016 
H  2.8844139643427678  -2.9153743302111366  3.97961871479311 
H  1.9649132841190247  -2.0108999565086743  4.841046080497065 
H  1.2423356722636898  -4.197100417352316  4.969090720574197 
H  0.14841633706479707  -3.3305691741493852  4.290692878131283 
H  1.952525603051933  -4.9543204394404015  2.924015917454294 
H  -0.05261235157311872  -3.4440129679414127  1.7103781236885642 
H  0.9676513299657958  -4.331536158613563  0.951629893936541 
H  2.4266048965344864  0.01826094122251788  1.9099603488706132 
H  1.2555063122072905  0.8223652503057052  3.7119554392914456 
H  2.3325580284605594  0.25226169817849764  4.672515588562615 
H  2.6909569001473215  1.4124746623372388  3.7046793541863217 
H  4.293518523594737  -1.0547502257496744  3.6890152340975595 
H  4.346248244591517  -1.2632193352617567  2.1541906034303 
H  4.572913638738168  0.14976849040351847  2.752576719260274 
H  -0.820750279773887  -5.347496971983437  3.13686128854387 
H  0.11070137733884633  -6.260432369320759  2.2966496064360182 
H  0.23455660943805823  -6.238276249287508  3.841611016254001 
H  3.5841920747158014  -0.8598226232923083  0.04069682848638217 
H  3.334852048420206  -1.660756361576886  -1.2636920813731098 
H  3.6163806901500055  0.4656407843285446  -1.9565297981373047 
H  2.086085578491847  0.22919549637219758  -2.048462148473992 
H  4.434357268387527  2.316409047686164  -0.710845928608952 
H  4.0350863834430735  4.141190731417297  0.6579272217542055 
H  1.9554643330169368  4.277400902925979  1.7535135335339298 
H  0.32110329406919913  2.713952444051922  1.295245794309271 
H  1.574621621599147  -3.568863556427862  -1.2019701448241211 
H  0.07693400746922308  -3.2838774390943395  -0.9085068534075826 
H  0.2680142141192199  -1.517949988662062  -2.6192215269976487 
H  1.5693468539839017  -2.242467177090926  -3.0532831172057597 
H  1.338424488421319  -4.115743860619941  -4.624711645020335 
H  -0.09088343373282703  -5.717643948369809  -5.6052658875457 
H  -2.3545292002808904  -5.662948348740268  -5.135058130303971 
H  -3.089009251109977  -4.207895625069488  -3.5212883244541904 
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
-C 0
6-31g(d)
****
-N 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****
-P 0
6-31+g(d)
****

-Pd 0
lanl2dz

