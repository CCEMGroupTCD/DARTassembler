%chk=EQAXEYOX_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3795191399542088  -1.4754751582117542  1.9318691422499145e-16 
N  1.3795191399542088  1.4754751582117542  -2.0844916858817484e-16 
C  2.5446926001918375  1.4442655341501665  -0.5268965998408409 
H  3.0964201108127054  2.179091724862893  -0.38402958661968595 
C  3.086533599369852  0.34060970779683825  -1.341232647530084 
C  2.671435186260038  -1.0032742125554333  -1.2066207756452942 
C  3.3169070405454186  -1.9882390131565904  -1.9531024478262924 
H  3.056450717183  -2.8771125110943836  -1.8632610971888135 
C  4.3496235925921205  -1.6583235501338374  -2.8356025613601634 
H  4.752818515445654  -2.3215972247554397  -3.348193408663471 
C  4.773865826386562  -0.343796816691102  -2.9487249911850593 
H  5.465943612094332  -0.12380704716262246  -3.528344711202149 
C  4.1573216996717814  0.6430721061071407  -2.189176915663442 
H  4.461834709540728  1.5204553800952836  -2.245272257238345 
C  1.1155335375675488  2.5393929417565118  1.017563423865946 
H  0.1481303343924374  2.555130539825017  1.1775523625260809 
C  1.533052272205603  3.9557274007742795  0.5987242107720759 
H  1.0920795337492992  4.198606813555671  -0.23023021185038253 
H  2.491473529677416  3.987605250250923  0.453848968101287 
C  1.1379988137966266  4.939586604991167  1.708489483826183 
H  1.3950630273400348  5.838664198820514  1.4505204733048702 
H  0.17502728076734742  4.924772617864138  1.826474959036755 
C  1.8159119649856852  4.575705354932908  3.0251943920349147 
H  1.532991385713978  5.195129524776394  3.7159766903303244 
H  2.7771619558327996  4.654469368242931  2.9266084127896996 
C  1.4622788325007237  3.146833846479202  3.4453267520574538 
H  1.9579487820191097  2.9157544832636075  4.246785103464688 
H  0.5159264141625606  3.098096129954611  3.654156025133626 
C  1.790680453162231  2.1423330036481403  2.3382661564717417 
H  1.4884828530086105  1.259819480629459  2.6031625742477282 
H  2.7522946030554056  2.105140784792684  2.2108949447248563 
C  1.0233450817314176  -3.203641739538945  -0.42551675341049794 
C  0.23948736448301267  -3.441860206345252  -1.5601166922359082 
H  -0.1504008659308611  -2.730323353856049  -2.0162156237568056 
C  0.04477038540302347  -4.747571650741873  -2.005307765854479 
H  -0.4668332476979724  -4.906998746949371  -2.7652981431571986 
C  0.6142976097988322  -5.812175226923916  -1.3151907784875252 
H  0.4920192908797234  -6.681205713722516  -1.6218180697336195 
C  1.364072080497972  -5.582261603429812  -0.17110970670123954 
H  1.7275862714546781  -6.298028313646722  0.29966772427688426 
C  1.573234795264812  -4.27609621713097  0.274434435204105 
H  2.0798941363852634  -4.121397967059014  1.038193365417091 
C  2.2200182266949393  -1.4882304051539046  1.611848218308257 
C  3.569264292292502  -1.1629876809174535  1.7610166962142546 
H  4.092182635053059  -0.9983619780392078  1.0094604676082868 
C  4.135475184444264  -1.0842713064982692  3.0316900051363866 
H  5.031743280058902  -0.8544168503845954  3.1273865556909195 
C  3.3606374207442644  -1.349267852067955  4.1553382778536765 
H  3.741222435532716  -1.3020886470419823  5.003482569310984 
C  2.017726735029126  -1.6843769403153357  4.01517523291612 
H  1.504837400886072  -1.8689604004191842  4.768525036910544 
C  1.4429707593566843  -1.7432133290656648  2.753722971517499 
H  0.5411925128040732  -1.9519324649223757  2.66316417088645 
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