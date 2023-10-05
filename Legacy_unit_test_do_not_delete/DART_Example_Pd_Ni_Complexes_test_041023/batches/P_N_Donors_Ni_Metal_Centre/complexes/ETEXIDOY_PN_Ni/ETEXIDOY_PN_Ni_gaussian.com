%chk=ETEXIDOY_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.1832764163541842  -1.6370573974360227  -2.9909897787488293e-15 
N  1.1832764163541845  1.6370573974360227  6.321855582898153e-16 
C  2.330443988997211  1.551926448320692  -0.7254889535885958 
C  3.038965598223384  2.7032424546818308  -1.0760113819098358 
H  3.820313687319449  2.630847809151614  -1.6125007005029615 
C  2.6173876824121125  3.929303155512233  -0.6558644685812893 
H  3.1310894150369615  4.705679965507176  -0.8505477981201125 
C  1.4186219394842623  4.034166968026607  0.06430537555802206 
H  1.080146189297586  4.878779788056969  0.3374053769786262 
C  0.7523912058083323  2.870873104303404  0.3605049470121146 
H  -0.060914446326389315  2.9350304577135513  0.8467784196683351 
C  2.9188743168015843  0.201473441320123  -1.0904622611100843 
H  3.4476945688685676  -0.11459882090501465  -0.3139371964433209 
H  3.5563141390988706  0.3440774302612757  -1.8351834406717442 
C  1.975439283722247  -0.9325006959993674  -1.4997535684318137 
H  2.483568445947027  -1.6371915619438433  -1.9744201755901365 
H  1.2809866504334901  -0.584390800223846  -2.1133656369095317 
C  2.6431240857316998  -2.3161543684653756  0.9101793538374027 
C  3.149649133212955  -1.6190514038488881  2.0310759962612677 
C  4.304525683593687  -2.0886891521673894  2.6464338098421254 
H  4.632232327804815  -1.6327984783468525  3.4138527168696813 
C  4.997018856457576  -3.203897774576947  2.181630427711255 
C  4.477108462012476  -3.855717898093235  1.1190987301477493 
H  4.93703551865411  -4.6242226872169585  0.8030825570443986 
C  3.3162574846599933  -3.480901498276699  0.4620934905884527 
C  2.525291549804991  -0.37245118800046495  2.5937248209180526 
H  3.0101274407178673  -0.10037269173173072  3.401268448687393 
H  1.588692510209897  -0.5496742441097969  2.8213302852919235 
H  2.570822002532325  0.3445568606106535  1.9285183220762714 
C  6.331888084486745  -3.617105555843172  2.807430214833237 
H  6.209873881064768  -4.4408669414317945  3.3241123354610336 
H  6.650048233198894  -2.9033675941626518  3.398736372274235 
H  6.99108497365665  -3.770320565526302  2.0976903173239285 
C  2.855951181867474  -4.292239287703175  -0.7330341591034487 
H  3.47310512454922  -5.038168047997436  -0.8795539038420377 
H  2.8386516742722208  -3.720584856731049  -1.5284696175706498 
H  1.955701591264913  -4.641094407807369  -0.5626087801335303 
C  0.022252585093617494  -3.0205972437159203  -0.30166843885440364 
C  -0.5025893368850842  -3.634119837103945  0.8722602788862357 
C  -1.504948183936474  -4.59692933142505  0.7306396257913509 
H  -1.8358410935400535  -5.031463820633638  1.507933083928178 
C  -2.031929344330467  -4.942195441755242  -0.5011145622052989 
C  -1.5086173019902998  -4.356902139577397  -1.6475604382795095 
H  -1.853533461544462  -4.608412255010654  -2.4962667978133366 
C  -0.48010346648482205  -3.3990054723890446  -1.5800901102582792 
C  -0.04555612947854937  -3.282443267362594  2.2603946123074397 
H  -0.6530744480861872  -3.6821390239223897  2.9162904765604867 
H  -0.0448794963389767  -2.307841987347088  2.3663279269904556 
H  0.8613279755018646  -3.6256823565234617  2.4022415643576602 
C  -3.201712066551349  -5.921860246124252  -0.5945914249798289 
H  -4.039985087119514  -5.422144535922973  -0.6981490030007791 
H  -3.242664835390358  -6.460579454370258  0.22298253827439604 
H  -3.074975025782555  -6.510765831542297  -1.3675484250520271 
C  -0.0376746177615328  -2.836091571848174  -2.897450852194716 
H  -0.29266933872259204  -1.8916241802823102  -2.9524060935349814 
H  -0.46397306865695076  -3.3338735990857664  -3.6252765759249157 
H  0.9374010989102128  -2.917371428136903  -2.975438939226041 
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

