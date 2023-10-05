%chk=UZIZEQOH_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3926231417723909  1.4631133876087656  3.560052170000538e-15 
C  2.8223683940623627  0.9812711900247522  1.0056805955082198 
C  3.5255478324341745  1.915772819617939  1.7236923166791158 
H  3.2319435885121504  2.796584837274497  1.7405777394680098 
C  4.673961867790661  1.5639511601624891  2.4306008986996384 
H  5.12790038654763  2.2047144550238773  2.928070527066134 
C  5.131579478038601  0.2822584452425425  2.388546204406409 
H  5.909152201875116  0.052458740142786064  2.841606437357303 
C  4.446580929760058  -0.6691639440194911  1.6796461860207565 
H  4.768262440084683  -1.541621971408776  1.6480989775980481 
C  3.277382764225765  -0.3445936311361375  1.0053660278475445 
C  2.6228243037801735  -1.4166849764033838  0.2687809212494025 
H  3.15300276469178  -2.1221043675052744  -0.024563624026019304 
C  1.0492974989838966  3.174787876551153  0.4679079791023959 
C  1.264810994915265  4.227736484668154  -0.3763948941225647 
H  1.6300085525358157  4.083794397346952  -1.2179418918628078 
C  0.9347593539245489  5.529457418349662  0.038897370547363336 
H  1.0850826026392126  6.251231560356477  -0.5310448009281324 
C  0.3958621161669068  5.734775090612211  1.2650372812402302 
H  0.16323798266344847  6.594768024622912  1.5333062991372532 
C  0.19153634162094235  4.680460118810119  2.1122569918762446 
H  -0.16053245803493565  4.833015099777655  2.959401647818093 
C  0.5019164838144714  3.399573995511843  1.7255220872272705 
H  0.3455303450471483  2.686967710974999  2.301212454386606 
C  2.015766116055553  1.4646944544075462  -1.6927133643441798 
C  3.365575378753984  1.6482290712555778  -1.9542483118215654 
H  3.96601320171848  1.7990596677554629  -1.2609603111226075 
C  3.8041934177295302  1.6010935104821822  -3.268329349960326 
H  4.70632580036915  1.7408480925478853  -3.4523472050304287 
C  2.9549102764096773  1.3608471440850536  -4.282794790132082 
H  3.277742169596705  1.3019180578432048  -5.152614214803674 
C  1.6088549798375669  1.2006911288475954  -4.041118561773944 
H  1.013796178871167  1.0606355174811337  -4.741563783841446 
C  1.1584871597614874  1.25334876676718  -2.719080808981432 
H  0.25320588719615045  1.1406208323522284  -2.5401859226793215 
C  0.9476285054095877  -2.5540917311319467  -0.8889275659871105 
H  -0.00019221769047916837  -2.49127402966165  -1.0240650265377633 
H  1.1593544058651903  -3.399882180785669  -0.48675490611225064 
H  1.3978570774421961  -2.481372572178798  -1.733688743389724 
N  1.3926231417723893  -1.463113387608766  2.220446049250313e-16 
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
