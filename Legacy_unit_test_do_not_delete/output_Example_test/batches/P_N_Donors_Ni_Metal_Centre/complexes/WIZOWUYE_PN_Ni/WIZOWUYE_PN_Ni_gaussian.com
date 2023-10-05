%chk=WIZOWUYE_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
C  2.9023021966446416  0.4779613527704973  0.4710341145537845 
H  3.6592356481868604  0.7730962753971785  -0.07482244845753026 
C  3.1600464746921055  0.8643134374329425  1.925681310203228 
C  4.163081553056004  0.38121387885700375  2.6739561514716526 
H  4.692628273549268  -0.2972159191131098  2.3236400941162065 
C  4.462912540010926  0.8920791178795551  4.053327963863124 
H  5.424329036057394  0.933131708467559  4.178358263933874 
H  4.101050199655224  0.2785385416777745  4.709698086573086 
C  3.868066260695441  2.2744274790997485  4.273053616021379 
H  4.352181474285809  2.92581480944314  3.7395031911893315 
H  3.9599627616359836  2.5227577513035375  5.206325606934065 
C  2.3905127949262317  2.293985097628476  3.882083546385713 
H  1.8914475173413814  1.705122085490106  4.4701019913723306 
H  2.038333923489108  3.191477561449388  3.9876129730777587 
C  2.211884207217044  1.8465120398227632  2.4493690702282627 
C  1.2618448330079863  2.284555714325719  1.585390328197891 
C  0.2140117778508388  3.2860531293218798  1.8330810956400017 
C  -0.47534383499386146  3.359921269306674  3.0526557046470546 
H  -0.3244853216977348  2.728853745666233  3.7199502195146716 
C  -1.3870937940789212  4.3912934927666445  3.2465341428078687 
H  -1.8419159832188565  4.475697405044666  4.053260401124839 
C  -1.6064388454815977  5.288088526236737  2.2197595247107045 
H  -2.2048341732290693  5.993479672058642  2.3235196959040634 
C  -0.9193159669050222  5.1134629513322185  1.0376368075929006 
H  -1.1021444198955932  5.701291561519234  0.33891524514997945 
C  2.686322287674809  -0.9926521644580752  0.2504513594118877 
C  3.7439982001465495  -1.877110549988679  0.2940270228224239 
H  4.611742211226295  -1.557173960028622  0.3998930213297895 
C  3.521373252630987  -3.2427488304728453  0.1814425229126509 
H  4.231440535866046  -3.843037976899829  0.207292938018264 
C  2.2133071627342904  -3.694262076768027  0.028659039242489033 
H  2.028324316341521  -4.6050250367700585  0.002572314097932191 
C  1.1926504675237104  -2.7642542025568724  -0.08413065189588702 
H  0.3242323373213569  -3.0659296896181973  -0.22362814879881843 
C  1.8496542853777682  2.579207416757804  -1.3470299694331997 
H  1.038936102216495  3.0752645522257764  -1.5862500267925561 
C  2.927023036267933  3.6111383944479694  -0.9645396402890514 
H  2.6297938064264255  4.124870903381496  -0.19725327641133011 
H  3.7467834808964504  3.15329927582208  -0.7238980709824798 
C  3.185939565282818  4.550386776999513  -2.1472532010793164 
H  3.884466832549916  5.179437512417846  -1.9132472841583565 
H  2.379482816076567  5.053046981013204  -2.343004882963769 
C  3.609091538452623  3.7489751077792253  -3.3917774652221016 
H  3.7363986547175942  4.356757358825836  -4.138022715432126 
H  4.452567273131519  3.3063566572872984  -3.2185173939649387 
C  2.5619037529309017  2.717663686010704  -3.7539812106250454 
H  2.870752662641693  2.200032871914345  -4.514262843487748 
H  1.7419665902367683  3.168021499067834  -4.012477998712464 
C  2.274252006449056  1.7745770441273347  -2.5778929832172617 
H  1.5685679869904037  1.1548007614243505  -2.8214647469329206 
H  3.0701418990161335  1.2595407556756046  -2.3715345075312593 
N  1.4144073378981024  -1.4420651450263953  1.243791779212181e-15 
N  -0.011553666914833283  4.156975091018008  0.8273291049568835 
P  1.4144073378981037  1.4420651450263953  -3.986466513268852e-16 
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
