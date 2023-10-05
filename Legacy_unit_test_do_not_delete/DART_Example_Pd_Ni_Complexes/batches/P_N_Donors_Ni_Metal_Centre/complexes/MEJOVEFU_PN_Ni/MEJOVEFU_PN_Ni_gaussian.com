%chk=MEJOVEFU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.4099740662508653  -1.4464000596308066  -1.7219015812851847e-15 
N  1.4099740662508653  1.4464000596308058  -1.398378247419013e-15 
N  1.4481426084277977  -2.4003919761910195  1.3414872476195636 
N  1.5487998439039223  -2.7009041907099745  -1.1036938744798535 
C  1.241770217959501  2.740755585637814  0.10604870444799497 
C  2.2976633421815036  3.658916360837851  0.1638500115680488 
C  3.562790891792705  3.212320624271801  0.10422223871165953 
C  3.8186562556948145  1.8197311283102664  0.01702096983992274 
C  5.115244207715403  1.2612378541837481  -0.016256581610501187 
C  5.289489216769365  -0.09331240800311182  -0.06740366392509754 
C  4.192500368022222  -0.9544149824304915  -0.07646799175487051 
C  2.904922914095818  -0.4389187267193868  -0.06817218026544875 
C  2.7061908453287042  0.949391562421887  -0.024560646148662053 
C  0.5426150083023205  -2.2874634410237205  2.477281032316415 
C  -0.1943190177411671  -3.6168309885408267  2.5144957453704846 
C  0.8691431407821831  -4.593115762312503  2.0387730027036426 
C  1.6594072229210066  -3.8425890579115842  0.9758003048122883 
C  1.1611111240311074  -3.9779108754797905  -0.4501532685978652 
C  1.334317029957216  -2.5244372939475754  -2.5011448707114186 
C  0.5866522800079332  -3.426972388790707  -3.2262822606312858 
C  0.39450460228695095  -3.2300864663214512  -4.603014784936546 
C  0.9434869686868406  -2.1455628177256405  -5.229608398866989 
C  1.6802354964243724  -1.2538900495917664  -4.5116219495205625 
C  1.8927961727358376  -1.438097152755193  -3.158947002598659 
H  0.33291694006658723  3.0458258924743817  0.15712743925976172 
H  2.1274228215511677  4.600595729509658  0.24069785964060011 
H  4.293960269101531  3.834289130847796  0.11898497942901373 
H  5.877761168242891  1.8444920751600211  -0.004524206688254154 
H  6.18087522258043  -0.4485080964844417  -0.09133004117834825 
H  4.327353627737512  -1.9048649086776597  -0.09085778122207398 
H  -0.08979484859434606  -1.574116263574515  2.3636884827139153 
H  1.0181590444923327  -2.165127121399143  3.302235245921543 
H  -0.9461629493217218  -3.6415932766515877  1.9180443461491894 
H  -0.4778779463700724  -3.8655036483662184  3.3972727569935905 
H  0.4749495677969293  -5.37506589600421  1.6453369908330375 
H  1.4553163920458052  -4.854866138689448  2.7525659068219115 
H  2.594111150080052  -4.055463986911927  1.0275531130121651 
H  1.5877329565948561  -4.732071574476466  -0.863616195676953 
H  0.20989439615954963  -4.107804602840364  -0.45180741784201883 
H  0.19678728275632684  -4.1898042814287875  -2.7930354516502898 
H  -0.12518361779539222  -3.866021630922305  -5.100299351690561 
H  0.8131049306921709  -2.0064783585813655  -6.170517006221023 
H  2.057347709500763  -0.4925523286288113  -4.958564522327845 
H  2.4259300438435036  -0.808505340354105  -2.667887092572108 
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
-C 0
6-31g(d)
****
-N 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****
-P 0
6-31+g(d)
****

-Ni 0
lanl2dz

