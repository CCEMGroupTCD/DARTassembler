%chk=SARARATA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.1237551901993599  -1.6784737926163749  1.8763724849470497e-16 
N  2.4136518658434745  0.3517457764888316  1.4227040582582875 
H  2.6351157885274  0.2016549343748742  2.2601396797021 
N  1.1237551901993599  1.6784737926163757  -2.0555375575803629e-16 
C  2.6351719558922504  -0.7013795763024261  0.45567521638588115 
H  3.314746311850334  -1.3231134214125782  0.8180936079841852 
H  3.0098375205341785  -0.2988224598479832  -0.36783041362732644 
C  1.8755971843929435  1.5711523284431368  1.1014635396120724 
C  2.1289630974090805  2.6762755082984397  1.9345815153322068 
H  2.6499893587501306  2.578412678552404  2.723115092009254 
C  1.6096302585458127  3.8938963099996355  1.5886273626721343 
H  1.7645312538834421  4.651408471102887  2.1422778096367274 
C  0.8567636034641016  4.022215963215232  0.42925323571191826 
H  0.50333176853495  4.863662985631029  0.16462807340454177 
C  0.6403903459587525  2.901839429323288  -0.31610268312203804 
H  0.11982264811650678  2.9850655277077314  -1.1059034671878363 
C  1.7582600437255074  -2.766955482373427  -1.3003989433690868 
C  1.2591151849053326  -2.6366164529751877  -2.5928132605789185 
H  0.591727828470697  -1.9887989834482631  -2.784091938551477 
C  1.7453196087964302  -3.4638166488362923  -3.607735806469044 
H  1.416642204854845  -3.3698129707961084  -4.494769809567693 
C  2.696850737486288  -4.41075250329635  -3.3277002218155576 
H  3.010634868733952  -4.984389562293352  -4.016789101077475 
C  3.1987691330283186  -4.535332848798229  -2.053654221181678 
H  3.8622283737609604  -5.1905534306374586  -1.8704436182659767 
C  2.745131144666939  -3.7166778801924405  -1.0345824141471935 
H  3.1029257662113556  -3.8002832844694163  -0.15892353925949274 
C  0.7519849687933391  -2.8073096840880787  1.3755344679800072 
C  0.007831657515379398  -3.957539817851779  1.10628724200357 
H  -0.195446427115068  -4.184686664991823  0.20612004078851803 
C  -0.434996397551759  -4.767706503791092  2.1272624741417414 
H  -0.9407722288386657  -5.545814467552728  1.9270718019895285 
C  -0.14743968729183665  -4.453402201084251  3.4307944621207427 
H  -0.47046293365365477  -5.0029020378152715  4.136326036001388 
C  0.6112577860715498  -3.3384750730991706  3.7146723111949025 
H  0.8271520384588418  -3.1331741762273437  4.616858613202323 
C  1.06153941180076  -2.513919607241909  2.6894382885173456 
H  1.5836830797488524  -1.7460045640222923  2.893445804111095 
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


