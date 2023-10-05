%chk=UCUBUGUJ_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.425530962483805  1.4310700454554977  -1.6122970180068557e-16 
N  1.425530962483804  -1.4310700454554977  -1.7763568394002505e-15 
C  2.324915910448582  -0.9457075978200882  1.0303216923853398 
C  3.0827618312573075  -1.766337649248638  1.8502588174174819 
H  2.9849736602039005  -2.6897995330304094  1.8044927513395603 
C  3.9859265012620395  -1.2087400372142008  2.737149499468422 
H  4.4982241191806605  -1.7592441075787297  3.2818888859524464 
C  4.128093717346493  0.16611107244881906  2.817290805479487 
H  4.753322668392806  0.5338574428457277  3.3999054086550733 
C  3.3525353163863345  0.9839917930444759  2.043690266714103 
H  3.432094832010885  1.9082021202748423  2.1224527864448035 
C  2.438430577388539  0.43738010028766405  1.1334814700727354 
C  1.606573228820258  -2.6054947917054636  -0.48846832363570347 
H  2.2318779722655435  -3.1500315860256642  -0.06778312368648493 
C  0.9112994920353867  -3.158170485331804  -1.6522042595816053 
C  0.5215028819513715  -2.372395845179873  -2.72802672061369 
H  0.6617599061699129  -1.4529514167408752  -2.7060845199527783 
C  -0.07632928560840413  -2.951226988093176  -3.8303901622355916 
H  -0.3378988512807257  -2.4228953020230164  -4.551035117680566 
C  -0.284053195966145  -4.309599313174765  -3.8644608757317123 
H  -0.7206889155269174  -4.690730416248336  -4.590906428181145 
C  0.14421004434055962  -5.103958218673464  -2.842154475352177 
H  0.030018746128693152  -6.024456033705793  -2.8932282448762576 
C  0.7513688454619427  -4.542443135919964  -1.7259209952678223 
H  1.0478178958931466  -5.085681513711988  -1.033656879924301 
C  2.5451905173153033  2.0085842192896264  -1.2906780505979254 
C  2.2178399780516127  1.8318098847255375  -2.6242357430444265 
H  1.4265419186126678  1.401752266683432  -2.8581274194506667 
C  3.065517287073384  2.291798857376053  -3.604610740078508 
H  2.851748447344521  2.162087690833597  -4.499929290336006 
C  4.237131941123234  2.9508473316244284  -3.2600246390724106 
H  4.791271095875939  3.28702252201972  -3.9272847221979315 
C  4.578222798014974  3.1069171974631113  -1.9583948129229145 
H  5.375943000473411  3.530819383052435  -1.7340888931564338 
C  3.74432263219211  2.6400504482626235  -0.9657761842786591 
H  3.9837772816411157  2.7445365183014983  -0.07312103136399184 
C  0.8506367872822924  2.8591304688357666  0.9622055157264366 
C  -0.00597679394688333  2.643829749644254  2.0150743170352095 
H  -0.2635297769488385  1.7762655699993455  2.2273194939977286 
C  -0.4920331160691216  3.705302068783254  2.7659754606131663 
H  -1.079301463340569  3.5504366932526934  3.4708162254351658 
C  -0.10308299214960681  4.970581634063839  2.463857465805898 
H  -0.4092076713783921  5.68191423670803  2.97988493521282 
C  0.7316883590828316  5.2094174591063585  1.4135077479319569 
H  0.9746077690952882  6.0817177800893685  1.203138072712029 
C  1.2192163007513033  4.159497289390444  0.6625599179755621 
H  1.7982901390466357  4.325711071696459  -0.045294762289719126 
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


