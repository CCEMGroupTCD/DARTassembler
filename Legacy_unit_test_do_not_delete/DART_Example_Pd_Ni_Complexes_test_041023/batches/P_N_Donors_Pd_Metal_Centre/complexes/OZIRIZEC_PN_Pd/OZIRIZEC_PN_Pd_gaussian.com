%chk=OZIRIZEC_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.4007255994662187  -1.5906186830915816  -5.307056023443476e-15 
N  1.4007255994662187  1.590618683091582  2.7249997033822894e-17 
C  0.9062472326214394  -2.947960570214185  -1.0645124734755276 
C  2.7223876471807  -2.624998272297794  2.2307963616307513 
H  3.510343801967127  -2.4960107329327825  1.7522733707277194 
C  4.816331557549371  -5.641424563528368  -1.6586244988320902 
H  4.292500202752225  -5.690974732803893  -2.4614580758928453 
H  5.635130327013216  -6.127533753756886  -1.780987811552425 
H  4.322520140714038  -6.024500712750813  -0.9295421108264654 
C  2.790939443507363  -3.1938618469675637  3.4951281672238497 
H  3.6095550058310133  -3.4439209115492284  3.858350145278758 
C  3.0922693144403115  -1.1700646734405549  -0.5596834152430779 
C  3.9578642557054513  -2.0341462409527766  -1.237513675958751 
C  3.4682692253306273  0.19274499806837014  -0.3994846753846846 
C  0.5741231790546304  -3.783165807724595  -3.303517601713889 
H  0.5781244160462942  -3.6488519037883247  -4.224643898869686 
C  5.899908278900823  -4.124995078799927  -0.16342546298716548 
H  6.7241551965168584  -4.5982865929777805  -0.3028976874897126 
H  6.088737646188746  -3.206095709304815  0.04114422878344776 
H  5.423244106156691  -4.523699983213202  0.5690707251387308 
C  1.52893626226439  -2.2654927377038225  1.6739231602305198 
C  0.4244337200120988  -2.9999760915201383  3.6755178412599068 
H  -0.3577288973071713  -3.1604632514621054  4.152472767392406 
C  3.8116003854123663  -3.5368239382588  -1.239400012606467 
H  3.2345356936033625  -3.7785907891627257  -1.9811498837027854 
H  3.375100332962432  -3.824126061550714  -0.42119499296700524 
C  0.7071217710988038  -4.274991927533195  -0.586111682019801 
H  0.8390718493855283  -4.473240325060564  0.3123633369528889 
C  0.9457400917472808  2.5903826943561463  1.060356321971982 
H  0.09268110059348644  2.9518644379770933  0.8083469637321626 
H  0.8600047347223148  2.122131781728638  1.8944960233286106 
H  1.5815286658173022  3.3027020528826077  1.1579258897196654 
C  5.070338617073058  -1.4938361690695607  -1.8861414517449102 
H  5.655639268053183  -2.054362720453382  -2.3418082523306594 
C  2.7062469953090575  1.0976493131631817  0.5244268255047615 
H  3.2552835212142375  1.8617921985199595  0.7575566905598523 
H  2.54156813084915  0.6003971757353079  1.3418551057982797 
N  5.10071304697014  -4.2383926717726315  -1.3909043292685102 
C  4.598122002038248  0.6851114359363124  -1.0410854576686541 
H  4.885722224456848  1.558055369488772  -0.8954601132442691 
C  0.4555310649887183  -5.288192446382633  -1.4530709401180562 
H  0.1742961417708817  -6.109884640051369  -1.1141326023241143 
C  0.8196539753984222  -2.7524001739315205  -2.444563606299728 
H  0.901778007269644  -1.8867966514556538  -2.779411938036041 
C  0.37994162840618206  -2.442271083758121  2.4155617634136486 
H  -0.4387047984592163  -2.1730814347545175  2.0686050793018382 
C  1.6467925301702553  -3.40238796232497  4.201295109182455 
H  1.683955594206301  -3.7791047589838325  5.0501706075840485 
C  1.6690087250676395  2.361289632094524  -1.276333848521016 
H  0.8453845544986714  2.690189280112727  -1.6444319158727632 
H  2.2419493085055198  3.0991314369694054  -1.0596554028259118 
H  2.1032130616310933  1.794527048086768  -1.917821734204271 
C  0.39830353754027414  -5.068104351906349  -2.7820921215790704 
H  0.2812770211043214  -5.807861559279887  -3.3329597722809208 
C  5.357728796557819  -0.16564904763213395  -1.8268228612063095 
H  6.04184309649291  0.19933135349809333  -2.3413566612769188 
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


