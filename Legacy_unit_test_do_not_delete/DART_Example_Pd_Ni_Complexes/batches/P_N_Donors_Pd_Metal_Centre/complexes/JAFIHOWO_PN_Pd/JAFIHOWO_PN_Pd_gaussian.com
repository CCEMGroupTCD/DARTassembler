%chk=JAFIHOWO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.6009657858617716  -1.3888875233437725  1.700896659838646e-16 
N  1.6009657858617716  1.388887523343773  -1.7008966598386464e-16 
N  2.9859161010003943  -0.43537320505566074  -0.318459588357823 
C  2.782804352556466  0.9322763152070668  -0.2963040615263248 
C  1.4490577272526584  2.8470637694116983  0.25098636842041033 
H  2.367370197165414  3.239712379060861  0.30077659781664906 
C  0.7830574697553179  3.069386970398198  1.6080696051147063 
H  1.352631726871544  2.6869809072194317  2.323257838011907 
H  -0.08932893460126379  2.60198482365192  1.6289119828769163 
C  0.5693189533422545  4.561692700821974  1.8670072492414351 
H  0.09002660992253575  4.68536013761087  2.724908470022314 
H  1.446683662090288  5.015385299689127  1.9359153298754725 
C  -0.24249963344115733  5.184287497352498  0.7343834319542768 
H  -1.134566222830217  4.758938881952938  0.6984738862393194 
H  -0.3700914329580953  6.150413718840538  0.9130793939939191 
C  0.4563963778815936  5.007594193628972  -0.6090938554519807 
H  1.314585480651504  5.500499197449466  -0.5999215549932075 
H  -0.10896121607743292  5.3902675066704955  -1.3257237295422921 
C  0.7233735717648648  3.5258094951493932  -0.9086093396894462 
H  -0.1368434263819509  3.063073723391061  -1.0731033871031208 
H  1.2731053228036555  3.4510758840718685  -1.7291109541151364 
C  4.235177319494729  -1.1613056675870033  -0.7122965569237909 
H  4.049086461310891  -2.117659701067099  -0.4861847811429843 
C  4.493133227430137  -1.1685868331457892  -2.2191296910794502 
H  4.8005699859083215  -0.27408284767619395  -2.5127372003761836 
H  3.6570214623719988  -1.3841429994789918  -2.7030913274589166 
C  5.5616925655838045  -2.2162430583935158  -2.533237800723731 
H  5.208284354928032  -3.1171016081108696  -2.3249735694159632 
H  5.771584628415812  -2.188910078319935  -3.5005527908693788 
C  6.838752710645005  -1.9829484974200557  -1.7364373152089716 
H  7.27043304398572  -1.1484503213477226  -2.049719095622072 
H  7.466449919560852  -2.7304428031508396  -1.9011038473863127 
C  6.567147008478969  -1.875517566194773  -0.23698574961338337 
H  7.408046549868304  -1.6412311754461533  0.23079577888167566 
H  6.2644060599630205  -2.753938931852055  0.10334540154452687 
C  5.507058131752513  -0.8217510333394393  0.07467380450765958 
H  5.835104012784535  0.07647597139060244  -0.18365675487631583 
H  5.312762411121959  -0.8152484081461046  1.046405346040756 
C  3.943985601425166  1.8229139796753244  -0.6359773573655563 
H  3.621735916738487  2.606209611405701  -1.1266283779814803 
H  4.583210387438736  1.3301163099603586  -1.1918955075184416 
H  4.3849662964944445  2.1134626441405855  0.18959987109565166 
C  1.4575041653188077  -2.5541080453798326  -1.370531283282961 
C  0.8283539374693721  -2.098361698493448  -2.5357008452338343 
H  0.4892524119239554  -1.2122498828749713  -2.576684446760531 
C  0.7002777704548182  -2.9417054288626647  -3.6312422035995535 
H  0.28851111216669767  -2.6280490940781145  -4.42759617186477 
C  1.1765412757060285  -4.245540307183317  -3.558399471358032 
H  1.0721148698916891  -4.828561823605688  -4.30090280950357 
C  1.801696690478697  -4.699331812547712  -2.4114746364017163 
H  2.1261141688282192  -5.591622400641112  -2.369761974759562 
C  1.9563986445230963  -3.8540197118432094  -1.317460688676462 
H  2.4012274037949064  -4.1619365105846935  -0.5367724538621248 
C  1.9737878562425495  -2.3501780508773944  1.47657429393552 
C  2.9288015731846433  -1.913782572616638  2.393007448620593 
H  3.3877010869087987  -1.093619873438839  2.2503180195057206 
C  3.206213394288752  -2.683801312093274  3.5187139078569953 
H  3.8637461387820307  -2.391714121315345  4.137678157957574 
C  2.5329871883536597  -3.8719041694270735  3.7434380187964287 
H  2.749587646346116  -4.405866432414757  4.49956224357867 
C  1.540824952209395  -4.282078870707758  2.861503536196649 
H  1.0587966673350482  -5.083620400377365  3.028693005948809 
C  1.2510542626999035  -3.5232301538134028  1.7350629330714724 
H  0.5633591690735942  -3.7991521910978774  1.140082229566642 
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


