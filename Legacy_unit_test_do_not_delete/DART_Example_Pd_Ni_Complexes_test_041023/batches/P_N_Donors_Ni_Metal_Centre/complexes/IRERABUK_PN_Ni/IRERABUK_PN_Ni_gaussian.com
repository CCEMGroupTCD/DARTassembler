%chk=IRERABUK_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
C  1.5928047820010465  3.1574853717474083  0.1234880714240848 
C  1.4128922704919362  3.751156732151612  1.383265200141044 
C  1.538520095814157  5.129159538098731  1.495185423097711 
H  1.3976001895542598  5.5442192936733115  2.3383796292697383 
C  1.8706085222195852  5.920840506430368  0.38309656809147824 
C  2.0257227693357143  5.316609641855731  -0.8410492127030371 
H  2.234613741489887  5.855751453240702  -1.5947646291362052 
C  1.8888070626126074  3.9317550276998015  -1.0192389343582033 
C  1.1118058294528177  2.9349679082036064  2.594349317238118 
H  0.9007959811026093  3.5300103694637577  3.344991291282136 
H  1.8929998545645645  2.387626898466016  2.821118349932541 
H  0.3455873558878202  2.3519736660973876  2.413406164732318 
C  2.02092483235037  7.427820090517111  0.5508146352590135 
H  1.8753223894359856  7.668444202932282  1.4898650892505914 
H  1.362018635998991  7.886640564199139  -0.009359888742183587 
H  2.9241293626831295  7.697705844143776  0.28066107452775385 
C  2.016491169387588  3.345700409846386  -2.395694634251879 
H  2.1716565910914998  4.064282072393226  -3.043940410295875 
H  1.1909179558701397  2.8701476670970445  -2.625819074470837 
H  2.770357089282838  2.7189277019633775  -2.4155701838920702 
C  2.912288846310666  0.5625626703139157  -0.26110676423823165 
C  4.24627504889699  1.1975690961835377  -0.31842046603496343 
C  4.759749138677261  1.79105592507787  0.8353670070523965 
H  4.236961752201704  1.8013960335247163  1.6283638549019914 
C  6.024571205532457  2.366218250918115  0.842249457945963 
H  6.371181377486595  2.7433250998289824  1.6424867122664732 
C  6.776864236062822  2.3876727783494225  -0.31554618441420906 
H  7.634355732120158  2.794445823203058  -0.31811604235718993 
C  6.2689298780147045  1.8063711916428973  -1.4888713453166527 
H  6.781080127777605  1.8256162630050687  -2.288331818946272 
C  5.0244606514500045  1.205398629408686  -1.4826445022018533 
H  4.69389804852946  0.7950022235739863  -2.2738152837834074 
C  2.7264803579660937  -0.8936258733968194  -0.29659761974657806 
C  3.7792773115637743  -1.7681845765020259  -0.563963285676987 
H  4.632048587351036  -1.4273192736169178  -0.8071307310335082 
C  3.5781513391840365  -3.1440073244711626  -0.4741915228283293 
H  4.292005286263043  -3.748364725974614  -0.6437889682328355 
C  2.3347484419690505  -3.610963684515544  -0.13784644257366432 
H  2.1804783761879856  -4.543988136111363  -0.052365687264502137 
C  1.3100197453160707  -2.716400459964525  0.07711233379238608 
H  0.44866464577374776  -3.055140268164091  0.2889893614130678 
N  1.476435892614373  -1.3784908614133058  -2.577521681806882e-15 
P  1.4764358926143726  1.3784908614133058  -6.129056519584308e-16 
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


