%chk=XUDALOBE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5967438194651034  -1.3937392779856643  -2.8724277793978085e-15 
C  4.2151025571944984  -1.0763382102340586  0.9774373088574078 
C  2.919451577301593  -0.5132406155861787  0.4803164611710992 
C  2.7379595497071625  0.9440750044023191  0.38542248002790935 
C  3.8858098736931086  1.7990228406447994  0.7684169852491983 
C  4.809380309341718  2.25419833151342  -0.18855665227087143 
C  5.956054093419832  2.96547297997762  0.20495835044207242 
C  6.164643206269691  3.267914671661316  1.4992170475857778 
C  5.261631851203761  2.8375798131454277  2.4827275875232586 
C  4.1570612663222475  2.1006119625410444  2.093382714447752 
N  1.5967438194651007  1.3937392779856643  -2.817061370256243e-16 
C  1.7903904427808754  -3.1677228085161375  0.22441753203580278 
C  2.365100592554168  -3.9530114290548615  -0.8278939219548785 
C  2.595473273932848  -5.307906788486228  -0.5654218040723299 
C  2.231985392846832  -5.931711282434733  0.6439764982389019 
C  1.69188918177374  -5.152284009229065  1.6103087695450522 
C  1.4378121280656808  -3.7797270107890193  1.4513835455801734 
C  2.7413366663869536  -3.421148596666762  -2.238812432809077 
C  1.6026312638707443  -2.6253312135167066  -2.883795431776238 
C  4.014339125845844  -2.5269490755043686  -2.186942840902263 
C  3.050823327023793  -4.54986175690259  -3.2411371251111047 
C  2.4900697155063494  -7.425940485051367  0.852912786243264 
C  3.1167319060641017  -8.10674138876943  -0.25587726981877934 
C  3.2592145234966625  -7.633610068288135  2.0900509105203633 
C  1.1043947798314968  -8.048640267362392  1.2254878843809833 
C  0.8371402580933972  -3.040549783546811  2.691445068792042 
C  0.5267724040512178  -4.042206038092259  3.8508581758885936 
C  -0.5268627893677924  -2.41941377649449  2.2780802211577433 
C  1.7555936606746858  -2.0169256003624754  3.2724933102376763 
C  0.7874970354169286  3.4724590792868364  1.0041957827964398 
C  0.36582248478302315  4.76191292720981  0.8118938394252402 
C  0.44312838122424014  5.394944638330946  -0.4400633647516893 
C  1.0016442635104492  4.735515059559001  -1.515481469432254 
C  1.4497085175223978  3.4234831904113703  -1.3832933889647085 
C  1.3292819588915492  2.806525690840318  -0.12744282950334776 
C  0.6503261465336188  2.8023117068529477  2.339761302141313 
C  1.1229994276071456  3.664509360302707  3.4779588402448334 
C  -0.8050347642262965  2.3322541068785445  2.5731510237197184 
C  2.049548052413173  2.681342403150986  -2.5638314333608174 
C  2.991972488804607  3.5657629647796436  -3.402457042235693 
C  0.9667044787143162  2.1036132376475627  -3.4613583389329348 
H  4.158775768763613  -2.033017537609459  0.9885119410235785 
H  4.927377524011513  -0.79259666958411  0.39977276108298576 
H  4.384257754359801  -0.7651779115684565  1.8660378590521058 
H  4.6496862111228  2.068930617473538  -1.1221531438712364 
H  6.600684328023608  3.2443969876025225  -0.4520510243808347 
H  6.943114545493454  3.7842114255563466  1.7473491412871502 
H  5.407699413537902  3.0527588783673267  3.4056503454635165 
H  3.560861849106658  1.756203182088393  2.7691407013475864 
H  3.01303828941453  -5.835924712371968  -1.260739342971786 
H  1.4577728816570203  -5.580562328157336  2.4471181980027197 
H  1.3756662810390945  -1.897354571094467  -2.3088006190553956 
H  1.872714723027801  -2.276210784261641  -3.7402978024481035 
H  0.8459878880763776  -3.1982793521055184  -2.988002160705442 
H  3.772001428462811  -5.086804194286993  -2.8789267650432016 
H  2.276785089683453  -5.10583926807665  -3.3447062044918927 
H  3.303510399763686  -4.183863364622764  -4.097039409322775 
H  3.85472124321032  -1.7852404565910451  -1.591320354613624 
H  4.7373908736661345  -3.067539698801526  -1.8522808054437512 
H  4.225103158627226  -2.1882661077098318  -3.062143474935194 
H  2.578459707400327  -7.9668561691451565  -1.043128376821083 
H  3.9983315455609736  -7.761575956756003  -0.40812433935374953 
H  3.169581515392804  -9.048527066143086  -0.06687772831762195 
H  4.146374267235606  -7.286781699715048  1.9595652720692183 
H  2.839696070263852  -7.165068637825042  2.8070474902449494 
H  3.3057998889184788  -8.576257825836311  2.3072546372791707 
H  0.7336098910902258  -7.573797966902516  1.9718170181919996 
H  0.5084077599722763  -7.976427912735019  0.4720165458427974 
H  1.2036113523107135  -8.977781485131331  1.454356491202306 
H  -1.0826325375826724  -3.105412409884397  1.91802255639228 
H  -0.9595625427936656  -2.0449929122463293  3.0504962229345547 
H  -0.39446561421840753  -1.7447349431914427  1.609056656621098 
H  1.9676527537884587  -1.3602568930539245  2.5980767699085723 
H  1.336870027031013  -1.586085220468858  4.014842780286253 
H  2.5588593348898656  -2.4420027000451703  3.5742425537379785 
H  1.3418345501589628  -4.45081129446637  4.147966123351652 
H  0.11247307498862535  -3.569169075287528  4.577584140581511 
H  -0.06365866714773949  -4.728540928777491  3.534833449022233 
H  -0.0017957342859227232  5.241244901670297  1.5627257535500878 
H  0.11567133568808385  6.286582306601341  -0.5514421412311721 
H  1.0875448986160237  5.167915630827251  -2.366825810375087 
H  1.2032814008951433  2.0139915114230393  2.3332662154347448 
H  2.028662040862552  3.9223169256597745  3.3264643842536405 
H  0.5635291621692513  4.437691078735329  3.5220635154213396 
H  1.0684208308241856  3.180365589184174  4.309522858021214 
H  -1.0810970640134667  1.7622618080691443  1.8479430287732905 
H  -0.8741605879315584  1.8276996101319627  3.398419768014511 
H  -1.3709814091491626  3.101434682423051  2.6223308702583963 
H  2.5684171033162775  1.939600951698088  -2.230477301525812 
H  3.6858450319507985  3.896207540945729  -2.8279172560152244 
H  3.385741133005693  3.0736601727917043  -4.122846198846033 
H  2.490142288553211  4.29476194690523  -3.7625717242128704 
H  0.3952138330634143  1.5278014357731986  -2.951051548114183 
H  0.4589838322493123  2.8292486530903083  -3.8182606911364454 
H  1.354582676701798  1.6081468789767817  -4.178535165769615 
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

