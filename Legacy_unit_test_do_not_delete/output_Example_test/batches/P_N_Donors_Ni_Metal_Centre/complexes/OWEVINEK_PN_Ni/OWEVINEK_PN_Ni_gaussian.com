%chk=OWEVINEK_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.566647818592296  1.2750351416725738  2.8416164763587913e-15 
N  3.723113394413012  -0.8306112575930589  0.20412486756878884 
N  1.5666478185922954  -1.2750351416725738  1.3322676295501878e-15 
C  4.992451982856414  -0.14038283859784428  0.3912333053794368 
H  4.849279994964789  0.8270546066022073  0.3385363272121923 
H  5.620400756446095  -0.41675335458044693  -0.30766711359506804 
H  5.361311518124024  -0.3678757723373023  1.269647705500917 
C  3.5702999094506107  -2.2015201077656545  0.23213510377086571 
H  4.264253542946546  -2.8446954550147048  0.32185114051382835 
C  2.248696634030126  -2.459697011476857  0.10516243431796896 
H  1.8551677664910493  -3.3242202383065296  0.09121644258464112 
C  2.483490369151496  -0.30313204308470265  0.06420026821918362 
C  2.0056188698353203  2.2909084636430403  1.5104086797550342 
C  3.295292412522359  3.102451366407771  1.3957536743462415 
H  3.441098653290157  3.6022320659138884  2.225667556959587 
H  3.2204812536822303  3.7282154994287344  0.6459018820982464 
H  4.050572207569587  2.495825388252083  1.2423028228657749 
C  0.8295424975776304  3.2485339400658315  1.783466775264488 
H  1.0304464588752567  3.7976816240124185  2.5702278823046316 
H  0.014560902075057225  2.72817546654983  1.9482295988388467 
H  0.6943358692706214  3.829956597191601  1.0071720107306295 
C  2.09996047675234  1.2956790628789137  2.6738639264486856 
H  2.3271797442040656  1.7760486385727632  3.4973889626476202 
H  2.7946414737883387  0.6325092196612598  2.4803242111348864 
H  1.2380828240350508  0.84230174322191  2.786523536634324 
C  1.9432960066399163  2.0844406719752184  -1.6376903278568884 
C  3.430777205912781  2.283093835893303  -1.9291114188491865 
H  3.5387802525318257  2.706682941281814  -2.8064772671229705 
H  3.88362083119069  1.4139376426402512  -1.9295389408135453 
H  3.8250735424336098  2.8562721006439658  -1.2397343987624327 
C  1.3645590889979862  1.141283259527354  -2.709512745752222 
H  1.5366898894660683  1.5137887572008977  -3.600074283060688 
H  0.3983327000769217  1.048963490328218  -2.5760965339651114 
H  1.7904512775373354  0.2612946197411956  -2.634337898971709 
C  1.2196298952708817  3.4253333886710573  -1.689207136454558 
H  1.4075026600642948  3.8645684510594913  -2.5446973784897984 
H  1.5315387000351885  3.9938362764025577  -0.9548398345194443 
H  0.25508240103208646  3.278896565953313  -1.6025745311029282 
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
