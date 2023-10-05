%chk=OYIJOFED_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.2279884975031314  1.6037905879509329  -2.3285434136825415e-15 
C  4.75592708231504  -0.30639598691486625  -0.3626422493286809 
C  3.2911887516630696  -0.38712956957382727  -0.11572250254854127 
C  2.865507888807899  0.9360187117184444  0.4959143821405725 
C  4.17161503322756  1.7762150129834482  0.35220188220890647 
C  5.061627460318958  1.259265845037985  1.5048756640671503 
C  5.391132721882465  -0.18800897839489372  1.0665870855491963 
C  4.898805006056908  1.1365408761556886  -0.9009180969192315 
C  4.162486035580703  1.29124512611453  -2.2337789781673796 
C  6.306662887454297  1.6433169391318505  -1.066960223446551 
C  5.360481559219685  -1.4331112535936135  -1.1824438547721035 
N  2.651982728636685  -1.4632379548992565  -0.3782575538447843 
N  1.227988497503132  -1.6037905879509333  1.1102230246251565e-16 
C  0.7472498946365274  -2.634779317042155  -0.9785483565980581 
C  1.230888803847746  -2.1936979967534493  1.3675829727161166 
C  1.3841272775374815  2.6736311350214663  -1.4348351979280616 
C  2.138504021845273  3.845835457782785  -1.378018973596765 
C  2.3310926010486153  4.608698389634228  -2.530290398378062 
C  1.7673332647074516  4.200265315565058  -3.738671333635909 
C  1.011101184511613  3.030018916751912  -3.795368132024791 
C  0.8195556327723412  2.266041444128482  -2.644542353319403 
C  0.7616071324323718  2.6620513688478002  1.3965289799822784 
C  0.9243947983139935  2.186816150655262  2.6975881151588315 
C  0.4775731451523848  2.9474506556616724  3.7793567457169517 
C  -0.13072678557905704  4.182328441377042  3.5587390025074885 
C  -0.2932632435222575  4.658573616862115  2.257082472286668 
C  0.15454011705080983  3.8980928725784896  1.1759866033493358 
H  2.599915716559736  0.9035368774993287  1.4174417330565854 
H  4.043115370636483  2.727640480537497  0.32458021673529774 
H  5.857866512481066  1.7875882122882922  1.5999628244783115 
H  4.57701505198516  1.2577272838202433  2.332782824634286 
H  6.341460729716343  -0.31827518744337335  1.0352621385152183 
H  4.999221957922114  -0.8300999831286768  1.6626824426289102 
H  3.218675967163451  1.1191763335469291  -2.2311587914279407 
H  4.585986940592703  0.7402596254894109  -2.8959957791887003 
H  4.319677588069497  2.2152149829854793  -2.438169925667615 
H  6.816528198912473  1.5969447877058838  -0.25461950198553124 
H  6.23692375599939  2.5583570890275906  -1.3491312134910523 
H  6.743588874697088  1.130185015186941  -1.7483218383769232 
H  4.896756810575687  -1.543105875912346  -2.0149603174057553 
H  5.219660656747513  -2.206458122766002  -0.6323600126460337 
H  6.2992589343072956  -1.323849183273146  -1.3494048083298407 
H  -0.2039747433536172  -2.7612028860063873  -0.9416508053122016 
H  1.187847411598138  -3.4725294261237947  -0.8171312794935212 
H  0.9939733615249139  -2.307219712972405  -1.8470522061394732 
H  0.311017706770105  -2.315986574497614  1.614278754022774 
H  1.6542515419215105  -1.6040919988730526  1.9972315804920435 
H  1.6792622951658398  -3.04296241458508  1.3725950469696602 
H  2.5279549700568387  4.125944490330457  -0.5464564510993003 
H  2.851211565598452  5.4150067628814975  -2.49072690577493 
H  1.8990443085097957  4.725350115492669  -4.5321160096061375 
H  0.6235055721880958  2.747951960256067  -4.627050080464631 
H  0.300500184584901  1.4598996449975625  -2.6833770208333307 
H  1.34305779973613  1.3358597725263275  2.848546631851016 
H  0.5890756137777252  2.6187248083092562  4.674595572384387 
H  -0.4376643557495612  4.705503652844034  4.302299683369713 
H  -0.7107269108956893  5.508459967511843  2.1046184194662656 
H  0.04393221089693333  4.224791198367661  0.2799058061983082 
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


