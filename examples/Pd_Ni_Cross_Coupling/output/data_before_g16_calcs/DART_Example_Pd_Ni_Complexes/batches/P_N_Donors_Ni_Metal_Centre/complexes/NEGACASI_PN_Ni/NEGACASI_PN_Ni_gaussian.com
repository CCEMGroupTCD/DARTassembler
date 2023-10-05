%chk=NEGACASI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.3300473957720467  1.5202216696916258  -1.7597714260673189e-15 
N  1.3300473957720464  -1.5202216696916262  4.440892098500626e-16 
C  1.3770823827656873  2.924049813670721  -1.1273266395196873 
C  0.550557079383757  4.02595453151348  -0.923406347088311 
C  0.5506455131676784  5.064101047554219  -1.833949181562316 
C  1.3346268560308698  5.015566934752783  -2.9577539394882706 
C  2.190978761734631  3.933354986455978  -3.157014980741096 
C  2.197768522417962  2.8927190966117653  -2.23226181696279 
C  1.460010635233429  2.0157356924280645  1.7699313403802752 
C  2.782823306799574  2.7554309667548003  1.730475199125053 
C  3.6540946013005406  2.1680333340259255  0.947775133682911 
C  3.048651164551056  0.9307166633712083  0.3095572700838972 
C  2.7336603144585556  -0.08233746926222163  1.459258839497606 
C  1.7077306108850543  0.6106939293848657  2.3737962820415826 
C  2.9836739485261683  3.942332539053801  2.6649920301089454 
C  5.118875092554527  2.464886446961451  0.7216067853617295 
C  2.353510190733549  -1.418377949735046  0.8618653862611183 
C  3.176724336689829  -2.530193123599033  1.098744446254356 
C  2.963209714255812  -3.714238267071516  0.43995915936180935 
C  1.9720776728170069  -3.8143784196789188  -0.47695466225347316 
C  1.1417709410441266  -2.6970794567292957  -0.6564968817015739 
H  -0.0006887806616882042  4.063739715177946  -0.17612615085786182 
H  0.011874960485646735  5.806745349045468  -1.682720602830416 
H  1.2947551685365133  5.701580014813371  -3.584217380179067 
H  2.7520320380200802  3.907815954370033  -3.8996275910632736 
H  2.763990755114131  2.1688177919926632  -2.3628265303095906 
H  0.6991934122470251  2.5139207302414164  2.133721473982287 
H  3.526246061413773  0.5680377948928761  -0.4654968201897985 
H  3.5545001448640745  -0.21143062647251742  1.9785179458995126 
H  2.054953353126597  0.6873527892889222  3.276756849546328 
H  0.8810324538537504  0.10314644444701959  2.4018164122039187 
H  2.156854187851053  4.145075097027284  3.1074454723416083 
H  3.270244555049299  4.704950516493348  2.1568808533346338 
H  3.652856241431696  3.723739858377562  3.3187217898847656 
H  5.475868814901872  1.8444470506520132  0.08074634992218166 
H  5.595082761415217  2.3763883037746805  1.5504367044469114 
H  5.216070576004578  3.3593024108358573  0.38907728056023266 
H  3.874409800640347  -2.4640015792958847  1.7104346479536372 
H  3.50215756050134  -4.449379445677854  0.6225341481871949 
H  1.8454800036587395  -4.594063148314865  -0.9684649714999721 
H  0.4303998212192085  -2.7642947990949445  -1.2513328655069684 
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
-H 0
6-31g(d,p)
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
-Br 0
6-31+g(d)
****

-Ni 0
lanl2dz


