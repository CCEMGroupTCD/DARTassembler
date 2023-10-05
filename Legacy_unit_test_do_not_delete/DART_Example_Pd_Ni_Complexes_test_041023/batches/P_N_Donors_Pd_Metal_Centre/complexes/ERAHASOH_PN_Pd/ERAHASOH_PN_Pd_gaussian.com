%chk=ERAHASOH_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5334272773757505  1.4631133876087656  4.8074036150194605e-15 
C  2.9631725296657216  0.9812711900247522  1.0056805955082175 
C  3.6663519680375334  1.915772819617939  1.7236923166791158 
H  3.3727477241155093  2.796584837274497  1.7405777394680089 
C  4.81476600339402  1.5639511601624891  2.4306008986996384 
H  5.268704522150992  2.2047144550238773  2.928070527066139 
C  5.272383613641962  0.2822584452425423  2.388546204406408 
H  6.049956337478477  0.052458740142786064  2.841606437357304 
C  4.587385065363417  -0.6691639440194911  1.6796461860207568 
H  4.909066575688044  -1.541621971408776  1.648098977598045 
C  3.418186899829124  -0.3445936311361375  1.0053660278475445 
C  2.763628439383531  -1.4166849764033838  0.2687809212494017 
H  3.293806900295141  -2.1221043675052744  -0.02456362402601836 
C  1.1901016345872535  3.174787876551153  0.46790797910239756 
C  1.4056151305186242  4.227736484668154  -0.37639489412256544 
H  1.7708126881391744  4.083794397346952  -1.2179418918628073 
C  1.0755634895279078  5.529457418349662  0.03889737054736342 
H  1.2258867382425724  6.251231560356477  -0.5310448009281338 
C  0.5366662517702663  5.734775090612211  1.2650372812402304 
H  0.3040421182668074  6.594768024622912  1.5333062991372532 
C  0.3323404772243017  4.680460118810119  2.1122569918762446 
H  -0.019728322431575407  4.833015099777655  2.9594016478180984 
C  0.6427206194178303  3.399573995511843  1.7255220872272705 
H  0.48633448065050966  2.686967710974999  2.301212454386601 
C  2.156570251658912  1.4646944544075462  -1.6927133643441794 
C  3.506379514357346  1.6482290712555778  -1.9542483118215668 
H  4.106817337321837  1.7990596677554629  -1.2609603111226082 
C  3.944997553332889  1.6010935104821822  -3.268329349960326 
H  4.847129935972511  1.7408480925478853  -3.4523472050304287 
C  3.095714412013038  1.3608471440850536  -4.282794790132083 
H  3.418546305200062  1.3019180578432048  -5.152614214803663 
C  1.7496591154409258  1.2006911288475954  -4.041118561773944 
H  1.1546003144745272  1.0606355174811337  -4.741563783841447 
C  1.2992912953648452  1.25334876676718  -2.7190808089814302 
H  0.3940100227995078  1.1406208323522284  -2.5401859226793215 
C  1.088432641012947  -2.5540917311319467  -0.8889275659871111 
H  0.14061191791287975  -2.49127402966165  -1.0240650265377638 
H  1.3001585414685484  -3.399882180785669  -0.48675490611225136 
H  1.538661213045557  -2.481372572178798  -1.733688743389725 
N  1.5334272773757482  -1.463113387608766  2.220446049250313e-16 
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


