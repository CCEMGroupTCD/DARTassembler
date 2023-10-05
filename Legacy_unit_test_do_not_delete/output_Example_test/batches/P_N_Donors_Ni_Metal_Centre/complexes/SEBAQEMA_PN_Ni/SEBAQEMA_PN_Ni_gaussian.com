%chk=SEBAQEMA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.2753892239234257  -1.566359577970525  7.908059252415185e-17 
N  1.2753892239234257  1.566359577970525  -2.4733487557879807e-16 
C  1.4926420454301421  2.5204897094032974  0.8119605679092236 
H  2.183379367305586  3.1303689851681478  0.5835853696917087 
C  0.7697387603875849  2.7714419846843685  2.0835305329547826 
H  0.046233818554626716  2.0906171281200057  2.193373256097576 
C  0.14380174730145967  4.179909932292191  2.030118584757432 
H  -0.33592225212366467  4.3551467737303735  2.8669999047654495 
H  0.8510415752989597  4.8476953467723405  1.9134064349972075 
H  -0.48200222130358816  4.23059898614461  1.2782983118878106 
C  1.7394071944006604  2.6735337315879097  3.2598851232269 
H  1.2559807597315706  2.8408416848099587  4.096697850814021 
H  2.134181304999411  1.7770274441750051  3.285435431054813 
H  2.449873087067389  3.3405873473060637  3.1553908692479427 
C  2.1611665364931962  1.4762380490538385  -1.1716981362801837 
H  2.6786848999115045  2.316265679685634  -1.2557221405883174 
H  1.6130630013983613  1.3720520452171951  -1.9894703548234136 
C  3.1300557114803564  0.29075458400814946  -1.0523489817126288 
H  3.8626593651700327  0.41220860559780736  -1.706254286202682 
H  3.5306022866883544  0.2975935053839671  -0.14718775284587224 
C  2.4818850377784374  -1.0659951207777023  -1.2890122579434757 
H  3.1941175759093667  -1.7514982262595302  -1.3410202765786867 
H  2.0206838259230473  -1.0470600209340073  -2.1654453667299025 
C  0.7901664948020956  -3.2460865493527695  -0.5338326935989355 
C  0.700102141502966  -3.5756519041155914  -1.8833361093800536 
H  0.9081552095428754  -2.919854404328905  -2.5384825435096405 
C  0.3126577875360206  -4.842046720413212  -2.2867251794274033 
H  0.2596655273672326  -5.05232581573228  -3.2114650630887454 
C  0.0013397802925441127  -5.8064123113314565  -1.3344081285357945 
H  -0.2632919986019433  -6.678104513558152  -1.6043788420888114 
C  0.07758458935075185  -5.486313639831243  0.008479993980568595 
H  -0.13914756501370507  -6.143533140690563  0.660107710514155 
C  0.46445338852440576  -4.22866547417035  0.41539962196100344 
H  0.5118961386486564  -4.024721740884091  1.3426044855774162 
C  2.3516134307833747  -1.842488399910884  1.4417646618857678 
C  3.39665840629721  -2.7629828223142763  1.3715824171400102 
H  3.506001376825329  -3.3010709068472277  0.5970800103782483 
C  4.282914984931564  -2.8912957862254154  2.444755730319067 
H  4.988260683504597  -3.5262205978632166  2.404388056727327 
C  4.133390435891615  -2.095498140186151  3.5649142127637345 
H  4.749921961063174  -2.172237116147048  4.2833130351724105 
C  3.0869924941112927  -1.1846615129740776  3.646247084143215 
H  2.9821668315952046  -0.6465043798002393  4.421441705185015 
C  2.196352841604105  -1.065899721226985  2.5868977998167746 
H  1.4757760355665435  -0.4496331700620867  2.64434390078136 
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