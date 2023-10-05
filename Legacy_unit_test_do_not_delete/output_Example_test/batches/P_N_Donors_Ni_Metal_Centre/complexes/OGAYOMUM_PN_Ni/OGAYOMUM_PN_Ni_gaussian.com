%chk=OGAYOMUM_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.2199979651622344  1.6098773136484659  -2.751401578964528e-15 
N  1.2199979651622364  -1.6098773136484659  -8.881784197001252e-16 
N  -1.4514978809708483  3.296048319964506  3.28421801675119 
C  1.5734620904537095  2.543016585979395  -1.5571024052918605 
C  1.4737006521073859  1.6473519053384846  -2.8263204934356914 
C  1.7095617231922384  2.5496399664752887  -4.080569296455545 
C  0.7192982296041398  3.688439401349694  -4.18344378501533 
C  0.8228822168168449  4.570687851563051  -2.961574328121796 
C  0.6258529367934713  3.758754803051614  -1.6797481123076623 
C  2.329615854164598  0.38481584574666394  -2.795876690819597 
C  1.8867299966022975  -0.6185409388391658  -3.8337602729216282 
C  3.833708532487641  0.6628609962715584  -2.940523534147217 
C  -0.21192274390553067  5.712221693564434  -3.0243918482092917 
C  2.8451216592993362  1.001238601435419  0.6226155047458155 
C  2.7159816648217427  -0.31519132277852613  1.4043429230064306 
C  2.4480705807333427  -1.5070599196789365  0.529545101366924 
C  3.4317545927843423  -2.425526465212773  0.23259753163384111 
C  3.1854185337561174  -3.484762337128636  -0.583420969179153 
C  1.9312407573325494  -3.573672674448451  -1.1660374785941832 
C  0.9777645517435525  -2.6404089205820025  -0.858182981965114 
C  0.7821047021354688  2.901029342674376  1.2084448262697294 
C  0.6514369378826574  2.3558779000398844  2.670828883205084 
C  -0.12605523267426766  3.316601355105434  3.509308684756455 
C  0.41470110293263895  4.15260420985085  4.41777804161917 
C  -0.42160099366438897  5.081043142636741  5.050205271674025 
C  -1.7421161076243281  5.077151516891794  4.809745256177707 
C  -2.196858575281504  4.187502636286014  3.892849249344882 
H  2.4640621390228974  2.87289047170828  -1.5041439970575463 
H  0.5714285822106726  1.356547611724948  -2.8764609431767334 
H  2.5860597445143894  2.915374330211139  -4.02974839720045 
H  1.6370558177521308  2.010899956508675  -4.8585608936707745 
H  0.9104496871603913  4.197100417352317  -4.961309937430588 
H  -0.1591275186708796  3.3305691741493857  -4.245148122392006 
H  1.6915790707793463  4.954320439440403  -2.942266209772347 
H  -0.2699820597992011  3.4440129679414127  -1.6593894259236959 
H  0.7761400360137742  4.331536158613563  -0.9367100940820985 
H  2.2007595969374423  -0.01826094122251698  -1.9453735051818084 
H  0.967485692509729  -0.8223652503057055  -3.7054001176998015 
H  2.010358232196373  -0.2522616981784971  -4.702963682501315 
H  2.402315774444301  -1.412474662337238  -3.7482249680495534 
H  4.004447830557852  1.0547502257496775  -3.7884989842125814 
H  4.110710037170201  1.2632193352617596  -2.2564495679648044 
H  4.316353978606801  -0.14976849040351548  -2.862381671316447 
H  -1.0874356036365447  5.3474969719834355  -3.05819598868495 
H  -0.1272283968772152  6.260432369320757  -2.2510033348680114 
H  -0.057366989666387314  6.238276249287507  -3.7993460812028417 
H  3.4228779603525137  0.8598226232923105  -0.11764789953618959 
H  3.21921234172242  1.6607563615768886  1.1946482532683271 
H  3.524749171421898  -0.4656407843285426  1.8772387041116 
H  1.9985946669560133  -0.2291954963721965  2.0225215809475587 
H  4.298753720331976  -2.316409047686161  0.6037667006444 
H  3.8519565663751205  -4.1411907314172955  -0.7502382771660407 
H  1.7353659546373565  -4.27740090292598  -1.7725794244750959 
H  0.11799383772952177  -2.713952444051924  -1.2575524718857338 
H  1.4579001351964864  3.5688635564278637  1.1943950720319005 
H  -0.04911684927617288  3.2838774390943395  0.9533790944740819 
H  0.20155003758463175  1.5179499886620622  2.656383043865059 
H  1.5172384718020806  2.242467177090927  3.0447643613065476 
H  1.3412988426299661  4.115743860619942  4.623294691656049 
H  -0.05291753420297329  5.717643948369808  5.65313373414664 
H  -2.3315943629850704  5.662948348740265  5.262212512731309 
H  -3.12194674243003  4.207895625069487  3.675058755956612 
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

-Ni 0
lanl2dz
