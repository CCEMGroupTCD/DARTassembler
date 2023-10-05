%chk=OYASACEG_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.3218704843894507  -1.5273370363151677  8.196243455592256e-17 
N  1.3218704843894507  1.5273370363151673  -1.8704484127425744e-16 
C  2.694895753248839  -0.4594862959606674  -0.4873722882951074 
C  3.920144058921787  -1.0266015720294808  -0.8797102970260283 
H  3.99351543488043  -1.949750023534098  -0.9625401688994097 
C  5.002109225049302  -0.24013505780286418  -1.1409321400769916 
H  5.812514067633567  -0.6223994081839561  -1.3917724230136133 
C  4.8848961544374685  1.114451949391976  -1.0299743732615394 
H  5.628387533155105  1.6502881437230372  -1.186967343223207 
C  3.690564449374337  1.7044119409167686  -0.6920408935926904 
H  3.6199162542845267  2.6304849041000775  -0.6517320931709625 
C  2.5844854984805914  0.9032105388003815  -0.41024945008079305 
C  1.5759469472549248  2.2411766074620694  1.2894811079200008 
H  1.8312610257972222  1.6076254612562972  1.9625191502619392 
H  0.7764392472365323  2.6957165017311855  1.566249137376326 
H  2.2809058362806374  2.8806231747739277  1.1673052005702615 
C  0.9233296109948611  2.5559826021886884  -1.0368043750113347 
H  1.5882052310430173  3.246935996827794  -1.073362770795082 
H  0.07651188293734057  2.9391761973705357  -0.7974130341923057 
H  0.8532751017275042  2.1336106995493256  -1.8955090457754673 
C  1.8102204790951824  -2.268757202968109  1.5859706993356777 
C  0.9291519438002962  -3.112217918821638  2.2412509957929836 
H  0.10632520896927944  -3.3011933860154046  1.8499634196469674 
C  1.2402945165036183  -3.6767042335646107  3.460276737899986 
H  0.6403906175530265  -4.248285514737613  3.882744856266948 
C  2.444859437775329  -3.3824921544257944  4.040390643034607 
H  2.6597464808042757  -3.7553289752848245  4.86505426157857 
C  3.3354912769911413  -2.5486329702914055  3.425080007656177 
H  4.149894705483733  -2.3630900722525463  3.834948094541289 
C  3.038202343436068  -1.9755362082880672  2.193449628294319 
H  3.6448315845353054  -1.4064008345416341  1.7785977629328766 
C  1.4205072811830215  -2.8474322098495755  -1.2464786541010235 
C  1.6413133557176205  -4.163871592378305  -0.8983015733274596 
H  1.6909831379124582  -4.404543046196276  -0.001948259028971006 
C  1.7888010637988636  -5.126077385898226  -1.8851820761653988 
H  1.9101533055175777  -6.016924260326558  -1.6498083738299152 
C  1.7552605929877936  -4.764710872236267  -3.214315944374908 
H  1.8699413133460672  -5.40909870648565  -3.87547833853349 
C  1.5532862987772147  -3.4561873031216934  -3.5605762706466297 
H  1.541139332780587  -3.211610447379868  -4.457647986203706 
C  1.3689835466747495  -2.5010704622969304  -2.5904960440981175 
H  1.2081650598446017  -1.6180619036530965  -2.8355251233017933 
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


