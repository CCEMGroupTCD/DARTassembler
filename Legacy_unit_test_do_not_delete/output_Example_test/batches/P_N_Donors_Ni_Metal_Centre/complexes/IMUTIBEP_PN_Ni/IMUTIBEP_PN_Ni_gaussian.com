%chk=IMUTIBEP_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3886098030044296  -1.4669229069722782  7.638818861401961e-17 
N  1.3886098030044294  1.4669229069722782  1.4642785047449287e-17 
C  2.8480983295993347  -1.0417979075357642  -0.9819165369785756 
C  3.5817207175296963  -2.022019296925502  -1.6461218589312383 
H  3.290237440544148  -2.92579760981339  -1.6144619411914023 
C  4.739376076067403  -1.6959510572609906  -2.356747487466434 
H  5.224194897909445  -2.3742036840834593  -2.811646832592459 
C  5.180033028854346  -0.3848137433241872  -2.39754790367887 
H  5.970860896469775  -0.166357822105106  -2.8759511585665605 
C  4.47129631764522  0.6114981371798505  -1.7427770499056134 
H  4.787010684023086  1.5072900590008889  -1.7612733106843792 
C  3.2907907156243152  0.30069876399328466  -1.0547918296979593 
C  2.6231549470452187  1.413902274715243  -0.36535039329446956 
H  3.157538012985304  2.174226641250667  -0.16840482596942072 
C  0.8840565237517153  2.624729139956611  0.8148388400524281 
H  0.14462182286325276  3.044583968095606  0.28861601043518265 
C  0.25524356176690044  2.0661756494455488  2.0841423880967223 
H  -0.08226731979057056  2.8042553261987737  2.6330315987705637 
H  -0.48569279837384016  1.470317672966  1.8462275249899882 
H  0.9288874177412068  1.5646111831038256  2.588665223468357 
C  1.91097288680668  3.7114701421067178  1.064413593454619 
C  2.6567626263381108  3.785861225974401  2.2353403721434857 
H  2.592090295727776  3.0944889458228264  2.8833541913568266 
C  3.497287700274233  4.871755769190571  2.4582222267167775 
H  4.010891601381454  4.91232701169904  3.256294965913368 
C  3.591413016488719  5.88946321509946  1.5308371290420637 
H  4.13494665912984  6.647916852786793  1.7089726639477996 
C  2.890713861759478  5.8003328154182245  0.33572953705947634 
H  2.988796040924998  6.471967681379602  -0.3287301479141099 
C  2.041361132010449  4.71903339779997  0.12109251805689938 
H  1.5426735802893963  4.670110898424647  -0.6861627605343019 
C  1.9277633348310501  -1.5063759482499917  1.7262289062151943 
C  3.2382132767477434  -1.7838634958677169  2.081701268921032 
H  3.8937507352496263  -1.9243701741807668  1.4090456335099075 
C  3.5926381778117897  -1.8560215073739068  3.425771271146089 
H  4.489264259216506  -2.0478553138149134  3.674116657354892 
C  2.6199498715720715  -1.6443903415338703  4.410404689505277 
H  2.8662121439771924  -1.6721777015963912  5.32728106929193 
C  1.3128836921787066  -1.4009860847028397  4.066080844204705 
H  0.6521685297243394  -1.2834046090982336  4.738656089584403 
C  0.9667864489763494  -1.3273058186941213  2.712356161468736 
H  0.06639392967010993  -1.1528624967010357  2.4660244773730833 
C  1.0280882939940303  -3.186817874191357  -0.43494941461664277 
C  1.1712547744437076  -4.191728369482886  0.4782627761301785 
H  1.4574773129248795  -3.995932855189081  1.3624840578508934 
C  0.8945814617779437  -5.513733817610127  0.10303446805567278 
H  0.9909320002112166  -6.211334908394196  0.7410719370812061 
C  0.4950072091423755  -5.808615990997698  -1.1484790718027549 
H  0.31931093441753045  -6.709763613491661  -1.3922255620775825 
C  0.343355438340073  -4.8164664526119685  -2.05781210974515 
H  0.0607160760464418  -5.03021358854508  -2.939521394465832 
C  0.5971192868541657  -3.4828637692659274  -1.7210281847137254 
H  0.47588828397004546  -2.7914084718000693  -2.3615878974229454 
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
