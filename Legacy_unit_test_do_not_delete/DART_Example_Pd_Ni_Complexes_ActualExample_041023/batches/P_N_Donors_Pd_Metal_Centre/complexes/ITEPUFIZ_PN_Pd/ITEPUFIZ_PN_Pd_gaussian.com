%chk=ITEPUFIZ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5512668540583225  -1.4441852884931352  -7.780543502646663e-16 
N  1.5512668540583212  1.4441852884931348  -1.0328235399107963e-17 
C  2.7513865203301053  1.0016428529121355  0.0647968312325035 
H  3.3797595219469683  1.5445323881033317  0.2883435677752536 
C  3.1046188840114235  -0.4156750849652871  -0.29796795512223895 
C  4.352168249975883  -0.9140504913170462  0.4351915867412495 
H  4.176306124957387  -0.9457281093931409  1.3983858254579353 
H  5.09861284714058  -0.3026232582357968  0.2605478006786226 
H  4.5841484756166855  -1.8106081528003526  0.1150748196901587 
C  3.3918232787037343  -0.381770318457737  -1.8178005223735103 
H  3.74414598727655  -1.2508066914217058  -2.102329012141316 
H  4.049642574336644  0.3172995225412698  -2.0104431752275946 
H  2.5607222845008826  -0.19135415289290797  -2.3036591348593456 
C  1.3067205960395962  2.845175124161038  0.23526085106593114 
C  0.7154857047167513  3.225395819784481  1.4465656582478907 
C  0.4483030855194998  4.5773061532827155  1.6337639706857665 
H  0.05977795287985965  4.859528402219532  2.453442848974366 
C  0.7285385284468422  5.526347021459733  0.6644266160681169 
C  1.2984342109114624  5.100884246775803  -0.531207758841912 
H  1.4912407327979609  5.744514705696947  -1.203961456317441 
C  1.5955117860109373  3.7594335614406917  -0.7746539812134703 
C  0.40812785778670624  2.2246461884764392  2.52031349613508 
H  0.043156802276634654  1.4111023763622041  2.1148245512666515 
H  -0.2483623887953743  2.603696383324482  3.1427262176698174 
H  1.228747716443357  2.0059834985043694  3.008523630803083 
C  0.42482645747417  6.979780307842347  0.900146873959806 
H  -0.4882894978869292  7.071702151983629  1.2417471032259768 
H  0.508286556590313  7.472670978223514  0.05592135378955046 
H  1.0574091402176506  7.343056377180659  1.5549870505977283 
C  2.1759885216895167  3.3112784423179167  -2.095925447847057 
H  3.11884121054665  3.075597305893145  -1.9758186243588436 
H  2.103122545405522  4.038345996160398  -2.7493901612399925 
H  1.6833614533383219  2.5288712581407506  -2.4221034049129457 
C  1.8362791431161853  -2.899423568335227  -1.028759658019429 
C  1.2768013968742933  -2.902127328047442  -2.3025787743556436 
H  0.6590192957411803  -2.2258167377220754  -2.5519467573312817 
C  1.6258579884730555  -3.9004311790057278  -3.214590006207175 
H  1.261861187946618  -3.8917894035358893  -4.091919464804838 
C  2.5054030185890466  -4.906222074795554  -2.8395803210727846 
H  2.7336310008681464  -5.591063859212264  -3.45850782766926 
C  3.0520343490316195  -4.915178463635032  -1.5703980718754411 
H  3.649498087612514  -5.608614778891306  -1.3184622709960023 
C  2.729232341086856  -3.917315622570055  -0.6647823343307604 
H  3.113648936850076  -3.9213921424230898  0.20383159768261325 
C  1.6213556728397034  -1.9778552380078331  1.7312229441833584 
C  1.9710910884345374  -1.0430344703702894  2.7205276351013286 
H  2.2255583950037785  -0.1630009593475518  2.467616723188878 
C  1.9490470399480688  -1.3878187677643874  4.055388642171644 
H  2.181950207487512  -0.7438672246972583  4.713725041538693 
C  1.5872693600864902  -2.6678962063934044  4.439378422591602 
H  1.578928489486337  -2.9044751708256564  5.360137209389624 
C  1.2379598588291614  -3.601417902495636  3.483098434702186 
H  0.9867731395874033  -4.479497149046572  3.74691524397161 
C  1.252939245751015  -3.259334636631964  2.134926718197362 
H  1.009457087280137  -3.9061922507091213  1.4821840387157852 
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


