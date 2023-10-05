%chk=CAKEBUPA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.146417727096018  -1.6630773869546769  6.1421150650580466e-15 
N  1.1464177270960179  1.6630773869546778  1.8376365058192276e-17 
C  0.6027598252866173  2.738004641707086  -0.8699289124926688 
H  1.3149113546271325  3.3768519336815714  -1.0846695858520528 
H  -0.12200745126535284  3.202846573335771  -0.4023713080852858 
H  0.25538678625594535  2.3452232834107374  -1.698613554759849 
C  1.6378732846405089  2.344808786764096  1.2206032580132526 
H  2.2990844058225095  3.0256788796907195  0.9708484586933182 
H  2.054128468928152  1.6892773012471745  1.8177208871699062 
H  0.8855893862619033  2.7738902658326623  1.6805460376724524 
C  2.3288981248220164  1.0822119638675884  -0.7132003289789421 
H  2.8752714637939754  1.8245248800816567  -1.0701217893626027 
H  2.001666883704061  0.5536888765508141  -1.4837652587607 
C  3.2071130652345987  0.20197325309982908  0.13472020415023847 
C  4.477376779636506  0.6396413710776941  0.4781455444785851 
H  4.767316532921624  1.4968369088093672  0.190863436511736 
C  5.326759127827245  -0.14688989334845468  1.2270401409754779 
H  6.1982301381185785  0.16580376064455232  1.4374264926580198 
C  4.913663885090749  -1.3854732185489054  1.6705337113503995 
H  5.490703826954327  -1.9199698342642524  2.2039262903041306 
C  3.6404724500247596  -1.8472882735935618  1.3332209712875105 
H  3.357001251577125  -2.7050734381265436  1.6289120113238762 
C  2.789301489446419  -1.0714138704409903  0.5723826662575601 
C  1.5710885247004538  -2.185570790939419  -1.7231893820629596 
C  2.8900915864391585  -2.2453793246734626  -2.1746650651474 
H  3.59223628131792  -1.9431603576292686  -1.6117285470457765 
C  3.18907182270842  -2.739640753184126  -3.429314258590344 
H  4.090937790261535  -2.7457579880643546  -3.725214657385148 
C  2.2046648787881233  -3.2188227451052973  -4.2505467514813535 
H  2.4138157530228215  -3.5824454662988967  -5.1027781099977085 
C  0.8793459537473218  -3.1597329121421125  -3.8022212017196755 
H  0.1839472819049055  -3.488286423413499  -4.358989822810242 
C  0.5670765932839682  -2.638040108304088  -2.5767843143446325 
H  -0.342660588902131  -2.583893815590138  -2.3056082675900043 
C  0.9343294253381372  -3.306551055973114  0.8143898215845743 
C  1.1435633784290526  -4.515602845007389  0.16614214954326273 
H  1.5181649806964348  -4.517276774929013  -0.7055509429711556 
C  0.8140060685861057  -5.721892907476184  0.7709548103880256 
H  0.9546273023080842  -6.537775367693278  0.3028606107480255 
C  0.2800445832664711  -5.743917943885988  2.053931353418517 
H  0.04205603120010193  -6.567437213490081  2.4628813768640674 
C  0.10204900953468798  -4.543939677297733  2.728918786966562 
H  -0.24400435526577202  -4.547432199681574  3.614708542494163 
C  0.42366844374804535  -3.342956854329176  2.1206487659472857 
H  0.296050636297148  -2.529557492633168  2.5962685969194177 
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

