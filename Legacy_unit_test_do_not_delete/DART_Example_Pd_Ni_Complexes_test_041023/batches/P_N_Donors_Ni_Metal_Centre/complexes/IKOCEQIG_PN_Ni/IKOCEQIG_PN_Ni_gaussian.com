%chk=IKOCEQIG_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.2618195255265305  1.5773114736791842  -1.7120097971980576e-15 
N  2.7887277298102666  0.9923172823621713  0.38948443077537787 
N  1.2618195255265299  -1.5773114736791842  3.3306690738754696e-16 
C  3.372055157240896  -0.1842593678229414  -0.30998982082102877 
C  2.3732144482102937  -1.1630363061526285  -0.9156879989404944 
C  3.118238998336568  -2.369482868391386  -1.4609056822298792 
C  4.1121180803870505  -1.9920602857664123  -2.5409652242863534 
C  5.08530932339818  -0.9475710187498412  -2.021224446962071 
C  4.346658411116298  0.22660812462504243  -1.4228135270678979 
C  3.816221532633806  1.871373165196645  0.9881107416015538 
C  4.629082065762024  1.2344275987914517  2.0774925062365814 
C  6.005037165833517  1.3069904862806467  2.0370438975272775 
C  6.785858917836148  0.8295072671652428  3.080749779020144 
C  6.168129894008141  0.25199113912034043  4.179905047837295 
C  4.80714672202182  0.15907646282108678  4.224634551504765 
C  4.034675185165209  0.6562432514267846  3.1683044472977735 
C  1.6956008818533248  -2.132768800668868  1.3031309039374341 
C  0.5323888588542031  -2.7032880291259818  2.0660788245132826 
C  -0.00139521175707924  -3.9346038011908115  1.722888475834813 
C  -1.1443663881506163  -4.412511885661718  2.3822968187688076 
C  -1.7193257154727686  -3.6636307580523884  3.376570281472056 
C  -1.1799159822961212  -2.463379589154178  3.708751890404239 
C  -0.05778393427785811  -1.9766850980468753  3.0711616543308145 
C  1.3133294330362522  2.4910458432379152  -1.5724641176559486 
C  2.3954140517628533  3.2911508421930447  -1.9151646845686963 
C  2.4020272737834114  3.98138994091746  -3.1082685492371014 
C  1.3566603359240292  3.8885767537677327  -3.9632616764151054 
C  0.29871867872461166  3.083811489467004  -3.6517908923322806 
C  0.27038563865216547  2.370635053632113  -2.475353703099421 
C  0.9600711397771649  2.8521513442785076  1.2566734285517462 
C  0.6603207992040804  2.4336581132235877  2.536778345692353 
C  0.45397951468426  3.3516706073712794  3.5466614007494583 
C  0.47561211918345403  4.685102237773394  3.253614531689689 
C  0.7364784190477132  5.112532809509114  1.996630931979309 
C  0.9775081915260524  4.200783960915812  0.9976481313927082 
H  0.9405824431084753  -2.192989965977732  -0.40460702045383357 
H  3.891243150232164  -0.6820966486032329  0.3546551970637012 
H  1.9600272695225316  -0.7129917861922164  -1.680724348025732 
H  3.5893668071742826  -2.8074588335243886  -0.7353625396769421 
H  2.47798995161161  -3.0004438527038753  -1.8261148338720847 
H  3.6385517853206166  -1.6386783765485444  -3.310679392282891 
H  4.601395133951148  -2.779628594708204  -2.823519402464978 
H  5.658989300661414  -1.3456402653618817  -1.3483269743936308 
H  5.648961676470639  -0.6378233892214884  -2.749176005018981 
H  3.851697481121869  0.6790229690663128  -2.1231925860839835 
H  4.989262661527455  0.8556988156890326  -1.0615394762835093 
H  3.3798449368508225  2.659749672152128  1.3477531055320724 
H  4.418869510810984  2.166292337955174  0.28583481581567816 
H  6.418487203074386  1.6852153010716298  1.2941729086053293 
H  7.7122813569762965  0.8970220363005996  3.043719248288547 
H  6.6803921681529825  -0.07191142770008851  4.885367215281507 
H  4.391186136879084  -0.2358882828209723  4.957049830774385 
H  3.1066848904490083  0.5956839728099923  3.207827563140992 
H  2.112226810540905  -1.4321511208758606  1.8286494310196781 
H  2.3539565813788865  -2.8279361936013188  1.152859952966482 
H  0.39830882013882185  -4.4448977338929785  1.0557629335764438 
H  -1.5102738522568329  -5.2345048126362075  2.1449664883951525 
H  -2.4748893747027894  -3.9766769656275134  3.820396743925326 
H  -1.5765289476874198  -1.9573701431325263  4.381344468582441 
H  0.298920841646407  -1.1553295721290424  3.3224072371318867 
H  3.1211586039772894  3.361062826503146  -1.3361832518949788 
H  3.1296664047008758  4.516609448343489  -3.3274952219489675 
H  1.3597600494210607  4.370796348750975  -4.759272175083569 
H  -0.41520962654291105  3.0149458362256483  -4.243843725060478 
H  -0.4455763015849519  1.8093857634490775  -2.2864577677250746 
H  0.5947019105968516  1.523838246350619  2.7208051679254397 
H  0.3031185201213167  3.0662878548132046  4.417495767814666 
H  0.3070022552531745  5.305707480513865  3.925708012143454 
H  0.754268156417177  6.023126889351667  1.8084300366572692 
H  1.155159953802665  4.5010940346738515  0.1353419439166634 
Br  -1.5839191898578664  1.5839191898578664  0.0 
C  -1.2937473829822488  -1.2398794889499034  -0.4598005675376139 
C  -1.601562203225019  -1.4938654694550033  -1.7919791134336869 
H  -1.1309417375171225  -1.0428217585005106  -2.480879663633267 
C  -2.60218398542082  -2.411889969420291  -2.121542410022959 
H  -2.8231586547736818  -2.5626138747145775  -3.03270050271447 
C  -3.2683208079706407  -3.09759503004543  -1.143064365570095 
H  -3.9596020908062854  -3.707191723557817  -1.3683280760515548 
C  -2.9347449854640226  -2.896087320001541  0.14198826503765638 
H  -3.3668350933519986  -3.405606940434718  0.8171022443470677 
C  -1.9675509241539992  -1.9516175659140573  0.5112421903797603 
H  -1.7755298364628027  -1.8022912010877177  1.4301324877763764 

-Ni 0
lanl2dz
F 1 1.0
3.130 1.0
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

-Ni 0
lanl2dz


