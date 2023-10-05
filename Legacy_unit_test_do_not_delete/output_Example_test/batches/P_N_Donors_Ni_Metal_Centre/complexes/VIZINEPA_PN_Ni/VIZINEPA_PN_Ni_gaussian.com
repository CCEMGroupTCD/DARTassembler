%chk=VIZINEPA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4478761160748537  -1.4084582892297535  -1.5920830166417833e-15 
C  2.8735433123291583  -0.9702478610476525  -1.0507021491159707 
H  3.6459507412330066  -1.5183620515091842  -0.7636665799109108 
H  2.6588881605124355  -1.2297752235511605  -1.9821964131404732 
N  3.2891879635529158  0.4242680788415748  -1.056928488646756 
C  2.767937435362362  1.4281124495351585  -0.31150461843068583 
N  1.4478761160748521  1.408458289229753  -5.055533009513243e-16 
C  0.9376769550357636  2.367122905232268  0.8150978040086836 
H  0.030723240392191098  2.2894721037592554  1.0873266462074043 
C  1.6511384391293553  3.4282782453860436  1.2609048864761063 
H  1.2535576420916987  4.091372212840874  1.812149405861045 
C  2.9861647252199335  3.511569572432834  0.8816063097769081 
H  3.5084043634647077  4.257230330036343  1.14954134820255 
C  3.54845795005195  2.5265262661388803  0.12434744780597942 
H  4.4669104884958255  2.576995623693569  -0.11172567343843974 
C  4.595113103385696  0.6054745561307069  -1.6918022727497624 
H  4.677511040860991  1.5282571830204044  -2.012807215718515 
H  4.679326521681789  -0.012393080883884878  -2.4489327920850554 
H  5.303487891814625  0.419931784294053  -1.040084036901281 
C  1.2491610395079351  -3.1739466278016146  -0.3055340156378305 
C  1.1814565847523026  -4.079003699599685  0.7634283099820084 
H  1.2731220466371784  -3.7725425936509613  1.6567760158827598 
C  0.9810129687129687  -5.4192461792206625  0.5047790636879255 
H  0.9418158822883661  -6.0346699987460735  1.2273552332206517 
C  0.8366599216911051  -5.875244837406866  -0.7855503909661649 
H  0.6875098671271034  -6.799313373871486  -0.9482858239231318 
C  0.9102527941855045  -4.9880192343800385  -1.8487851291363075 
H  0.8168312902536226  -5.305359057943029  -2.7393887833979376 
C  1.1188964523427012  -3.641290357490601  -1.6108048441064282 
H  1.173538369983006  -3.0333788761585962  -2.339460648791567 
C  1.954617594594641  -1.1790842300400584  1.7177792480171885 
C  3.275512111920141  -1.4195820442744036  2.0995790980014366 
H  3.9073961345462127  -1.7377237907687555  1.4654950727433953 
C  3.65804207583638  -1.190770806153718  3.4135746132205633 
H  4.560818182357042  -1.3336706928873672  3.6730978951567166 
C  2.7382278537242426  -0.7579017608563835  4.346103586724151 
H  3.0073917463742954  -0.6156131417808421  5.24613796816469 
C  1.4261603478524199  -0.5299422412488948  3.975856010829759 
H  0.7946454430132703  -0.23828686306609626  4.622778885168762 
C  1.0266025512038746  -0.7272253252598221  2.654894261988953 
H  0.12833953490484507  -0.5546601404342805  2.394994010862631 
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
