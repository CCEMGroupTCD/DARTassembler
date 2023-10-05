%chk=TAXOLAKA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.46320914089545  1.392522534826635  0.0 
N  1.46320914089545  -1.392522534826635  0.0 
C  2.6782370118806766  -1.0064872308687556  0.06661851715014255 
C  2.9704561271158063  0.39486579660159915  0.2442533422573582 
C  4.158960464606878  0.8973306418740017  0.6299750152682331 
C  5.431146147029731  0.21068785076180063  0.7925922248494599 
C  5.906903698020687  -0.7517503909768517  -0.091024187126992 
C  7.152083587284637  -1.3177095336761946  0.10230644202323745 
C  7.915795860160932  -0.9627639743488303  1.1763724707998686 
C  7.473728818236916  -0.02358411384997039  2.0497625784138287 
C  6.250277005029873  0.5758302648221241  1.8604329362508312 
C  1.5265123352321373  2.79649786815465  1.1298630355017745 
C  1.6972224103482727  2.548815046962144  2.4836679218646522 
C  1.7162394956405116  3.598723068868602  3.3854740795301757 
C  1.5170589652067072  4.888670272648415  2.933643472159751 
C  1.336183435455697  5.134377283804682  1.6142493356690604 
C  1.3450299601645364  4.099562611897489  0.6930669032452772 
C  1.5637259710725653  2.015500560893553  -1.702363496991451 
C  2.7544176577286343  2.0734494399830004  -2.3922215148041412 
C  2.770365069709473  2.499341707553259  -3.7063304853522383 
C  1.6271252749724836  2.8753695109849566  -4.324724819718333 
C  0.44318309272547785  2.8383057804461718  -3.6551146353761435 
C  0.3928836977910901  2.39442985417475  -2.3428601962729827 
H  1.3090538846964948  -2.2275207132545045  -0.1329926276467058 
H  3.370247993680824  -1.6238392007649018  0.0003795079987981356 
H  4.167950075332783  1.8092442578056058  0.8165859312874849 
H  5.385832884464599  -1.014241827099648  -0.8152706805408034 
H  7.473650173337953  -1.9436181909522574  -0.5060263069562364 
H  8.742848238525802  -1.3667493002093112  1.310473738245149 
H  8.000574552639993  0.21573170450516077  2.77871237176472 
H  5.96481047560534  1.232652755751415  2.453283411486989 
H  1.799798175697201  1.6737514917775922  2.7870239944111823 
H  1.8618585084211328  3.4377916552682475  4.289030487970418 
H  1.506268072132328  5.593097608987446  3.541154216463879 
H  1.2040844268802302  6.007474587146609  1.3230539359278646 
H  1.2311917912103374  4.278009971839809  -0.2128705951610916 
H  3.5476695991514755  1.8251515407730277  -1.9744030366622336 
H  3.575143018970934  2.530489935977008  -4.170570391272198 
H  1.6501293459656226  3.159160659865142  -5.210252332862356 
H  -0.3366144194031986  3.1121710935468165  -4.079917072231492 
H  -0.4209061878782956  2.351008273203691  -1.895321183126419 
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


