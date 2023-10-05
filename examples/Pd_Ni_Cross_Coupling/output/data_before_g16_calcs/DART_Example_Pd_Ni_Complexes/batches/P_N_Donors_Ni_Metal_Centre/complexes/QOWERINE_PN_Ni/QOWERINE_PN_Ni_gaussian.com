%chk=QOWERINE_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.4109171476383713  -1.445480128711564  -5.49862016810709e-15 
N  2.729691373469656  -0.4631437362592668  -0.4932524945167039 
C  4.055383527036209  -0.6682204987217341  -0.9170603937706131 
H  4.458659661359137  -1.5149516060465775  -1.0665219987867613 
C  4.670345146055349  0.5332372076169903  -1.0803353673968963 
H  5.571917346047121  0.6651142125052666  -1.3477926233905984 
C  3.722650652781111  1.561086591699715  -0.7783393391804992 
C  3.6686334914080483  2.9571923630040855  -0.7571008927254784 
H  4.419595690279644  3.4816841702181875  -1.00837128499493 
C  2.4841439618078014  3.550909950643956  -0.35821402857161316 
H  2.425122566361024  4.499689510292723  -0.3460209128004625 
C  1.390462376221661  2.8054969189474113  0.021131554672456674 
H  0.6016834860061988  3.251435989430907  0.3043958442798722 
N  1.4109171476383717  1.445480128711564  2.670689485642905e-16 
C  2.547800106527677  0.8942529902864789  -0.40194891423017376 
C  2.0049679299634864  -2.4405052324340386  1.390903375860829 
C  1.3947464814873765  -2.3119998062366056  2.644038930881741 
H  0.6095468837501156  -1.7867185001028507  2.7378515160698713 
C  1.9327262855885556  -2.944729351827268  3.7397505318588977 
H  1.5242388170884915  -2.8429328988167675  4.591397142979311 
C  3.0552159102940752  -3.724513845305364  3.609172945116077 
H  3.4100314337869815  -4.171751242180709  4.368903104530998 
C  3.6765107883958876  -3.8626790687719885  2.3728114019302216 
H  4.460592711689539  -4.391161598435621  2.2860813990218865 
C  3.1418805739917692  -3.2259710769296124  1.2733779446576659 
H  3.5562729140563136  -3.3257557554838346  0.4244323065737656 
C  1.1356227884420695  -2.5734537601812146  -1.3840877257254705 
C  0.8230253342615339  -3.905557662374827  -1.1703682942261355 
H  0.8173643771163788  -4.257708567039095  -0.28805462896397677 
C  0.5168705304391775  -4.725269273392368  -2.253333686498468 
H  0.2929113470356741  -5.6360965333793605  -2.108212197277116 
C  0.5368537522208526  -4.219049315424331  -3.5390994797135837 
H  0.32150375620501315  -4.778709862761698  -4.276028398246445 
C  0.8734737337125043  -2.887504410884129  -3.7473265799036932 
H  0.9092314070446229  -2.5425999392852088  -4.631098628719979 
C  1.15653459240059  -2.064200156220473  -2.6817885399766745 
H  1.3655609121890668  -1.1500921973806864  -2.829058656502212 
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


