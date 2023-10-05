%chk=VOBOHORI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4642750578357884  1.3914016512136245  8.001339902573492e-16 
N  1.4642750578357877  -1.391401651213624  2.220446049250313e-16 
N  1.8215774118790935  -2.624512527935286  0.21733638329982896 
N  2.471247377377032  -3.010062882523473  -0.8867006405389866 
C  1.0258363787385432  2.9162986054103857  -0.8883496855140606 
C  1.0548490351401825  2.6064283715897387  -2.2055682376004313 
C  0.7948116580145483  3.370577202427054  -3.4762630659631126 
H  1.58818080094803  3.9017067503583847  -3.740286921618714 
H  0.01888486715954385  3.9773047560875363  -3.375692539154107 
C  0.5031892911395945  2.261461618401441  -4.505286732556735 
H  -0.4640286976604713  2.049550266461604  -4.523476594035491 
H  0.7844902489200021  2.5433673964044146  -5.41155208193316 
C  1.317251804611222  1.0409297802232036  -4.039557244384108 
H  0.8597584165858803  0.19524757985089547  -4.277688082132082 
H  2.2219208873067418  1.0429696081891429  -4.4413196989723165 
C  1.384663263513846  1.2172638422002011  -2.558475881394455 
C  1.664181450627511  0.4371958768838664  -1.5071095263307146 
C  1.8854670147067727  -0.9650313292562321  -1.2411737606077529 
C  2.532111049905301  -2.0191856805326247  -1.8146635298514802 
H  2.9371853907426395  -2.053758303727403  -2.6736881000010744 
C  0.9152577625550874  4.217669175454358  -0.23404619395499007 
C  0.699799746064368  4.309306076193112  1.1407047840599076 
H  0.5552178001909781  3.5189157679447622  1.6474050520236696 
C  0.6910044339586778  5.539735972913963  1.7759053994094949 
H  0.5425646508702912  5.585324890351766  2.7132841336485303 
C  0.8974136245043756  6.697814546843015  1.060598225766773 
H  0.9024005006086319  7.538177061132195  1.5025964900668554 
C  1.0956859918121031  6.62908170114272  -0.30634456904366625 
H  1.2187867311392313  7.426886105091619  -0.8064813813429255 
C  1.1165874096891957  5.409989871317393  -0.9474731218190654 
H  1.2687628898711423  5.375703538569123  -1.8851278027214122 
C  3.0671508047024827  -4.307699624168242  -0.9210350552760765 
C  3.092914468696624  -5.018971376028521  -2.0940297349840744 
H  2.728324373784364  -4.648214345123721  -2.889058072042909 
C  3.6534294631257556  -6.2804290316731475  -2.1030088237848474 
H  3.670343395022634  -6.782502199683744  -2.9085600591954943 
C  4.189496961885494  -6.817062630137492  -0.9478613060310632 
H  4.569233616022262  -7.687752799969697  -0.9573133778982997 
C  4.170808267806969  -6.081798571691346  0.21725149411792066 
H  4.541328370561456  -6.4474528560284625  1.0121962581114599 
C  3.616505569686197  -4.814482625600615  0.24145483356120812 
H  3.613707158516336  -4.302702910966174  1.04204270043382 
C  3.0932694518866093  1.6365131613965596  0.7293192934009838 
C  4.194665058405944  1.8612586256032784  -0.08908158285617174 
H  4.093061331362682  1.8830105510617545  -1.0332414614730334 
C  5.434674905568666  2.04984434401976  0.4800401227120086 
H  6.188256028259459  2.2159485417533276  -0.07337272229559333 
C  5.586326475569925  2.0010408084333458  1.85312516240565 
H  6.445129571284984  2.1319841955143395  2.238297243649855 
C  4.503588170222733  1.7635478496519803  2.665432066637792 
H  4.614912502872794  1.7209503294116515  3.6073601287685495 
C  3.250525289401998  1.5835024725450229  2.106296896136379 
H  2.4982488386051305  1.4248867846827142  2.664360553025511 
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
