%chk=IMIYOBOZ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.4007255994662187  1.5906186830915816  5.501850631334685e-15 
N  1.4007255994662187  -1.590618683091582  -2.220446049250313e-16 
C  0.9062472326214394  2.947960570214185  1.064512473475528 
C  2.7223876471807  2.6249982722977943  -2.230796361630751 
H  3.510343801967127  2.4960107329327825  -1.7522733707277192 
C  4.816331557549371  5.641424563528368  1.6586244988320908 
H  4.292500202752225  5.690974732803893  2.461458075892846 
H  5.635130327013216  6.127533753756886  1.7809878115524256 
H  4.322520140714038  6.024500712750813  0.9295421108264662 
C  2.790939443507363  3.193861846967564  -3.4951281672238492 
H  3.6095550058310133  3.443920911549229  -3.8583501452787576 
C  3.0922693144403115  1.1700646734405549  0.559683415243078 
C  3.9578642557054513  2.0341462409527766  1.2375136759587513 
C  3.4682692253306273  -0.1927449980683702  0.3994846753846846 
C  0.5741231790546304  3.7831658077245947  3.3035176017138896 
H  0.5781244160462942  3.6488519037883242  4.224643898869687 
C  5.899908278900823  4.124995078799927  0.16342546298716598 
H  6.7241551965168584  4.5982865929777805  0.30289768748971313 
H  6.088737646188746  3.206095709304815  -0.04114422878344737 
H  5.423244106156691  4.523699983213202  -0.5690707251387302 
C  1.52893626226439  2.2654927377038225  -1.6739231602305196 
C  0.4244337200120988  2.9999760915201388  -3.6755178412599063 
H  -0.3577288973071713  3.160463251462106  -4.152472767392406 
C  3.8116003854123663  3.5368239382588  1.2394000126064675 
H  3.2345356936033625  3.7785907891627253  1.9811498837027859 
H  3.375100332962432  3.824126061550714  0.4211949929670057 
C  0.7071217710988038  4.274991927533195  0.5861116820198016 
H  0.8390718493855283  4.473240325060564  -0.31236333695288837 
C  0.9457400917472808  -2.5903826943561463  -1.0603563219719823 
H  0.09268110059348644  -2.9518644379770933  -0.808346963732163 
H  0.8600047347223148  -2.1221317817286374  -1.8944960233286108 
H  1.5815286658173022  -3.3027020528826077  -1.1579258897196658 
C  5.070338617073058  1.4938361690695605  1.8861414517449104 
H  5.655639268053183  2.0543627204533816  2.34180825233066 
C  2.7062469953090575  -1.0976493131631817  -0.5244268255047616 
H  3.2552835212142375  -1.8617921985199595  -0.7575566905598525 
H  2.54156813084915  -0.6003971757353078  -1.3418551057982797 
N  5.10071304697014  4.2383926717726315  1.3909043292685106 
C  4.598122002038248  -0.6851114359363125  1.0410854576686541 
H  4.885722224456848  -1.558055369488772  0.8954601132442689 
C  0.4555310649887183  5.288192446382633  1.4530709401180568 
H  0.1742961417708817  6.109884640051369  1.114132602324115 
C  0.8196539753984222  2.75240017393152  2.4445636062997282 
H  0.901778007269644  1.8867966514556533  2.7794119380360414 
C  0.37994162840618206  2.4422710837581216  -2.415561763413648 
H  -0.4387047984592163  2.173081434754518  -2.068605079301838 
C  1.6467925301702553  3.4023879623249704  -4.201295109182455 
H  1.683955594206301  3.779104758983833  -5.050170607584048 
C  1.6690087250676395  -2.361289632094524  1.2763338485210158 
H  0.8453845544986714  -2.690189280112727  1.644431915872763 
H  2.2419493085055198  -3.0991314369694054  1.0596554028259113 
H  2.1032130616310933  -1.7945270480867683  1.9178217342042707 
C  0.39830353754027414  5.068104351906349  2.782092121579071 
H  0.2812770211043214  5.807861559279887  3.3329597722809217 
C  5.357728796557819  0.16564904763213373  1.8268228612063095 
H  6.04184309649291  -0.1993313534980936  2.3413566612769188 
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
