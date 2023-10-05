%chk=ODELUMUN_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.4338262830622126  1.4227586548673676  2.055286847288396e-15 
N  1.4338262830622124  -1.4227586548673676  -7.771561172376096e-16 
N  2.703631584638671  -1.2706913177544756  -0.4896859616854054 
C  1.5184671355858153  -2.233176033848266  1.0733667117671524 
C  2.837504798782566  -2.612295238749735  1.2492492906184305 
H  3.1679717862045207  -3.1998779796032304  1.9187819207981307 
C  3.579043854627657  -1.983130111595089  0.28567943771263676 
C  0.31473261997965785  -2.6028624934600924  1.888195655483634 
H  -0.48929662225835324  -2.544925752244003  1.3298619276438854 
H  0.23119918776191706  -1.986795596238012  2.645626757201207 
H  0.4152732997800668  -3.519473631644714  2.2234997608285556 
C  5.058411931700457  -2.0275131626374927  0.0773029819453678 
H  5.421929339957911  -1.1188113982328167  0.12574363905916683 
H  5.2530596218996894  -2.409423965387556  -0.8040362917982468 
H  5.469393905736998  -2.581731017719531  0.7726197415055736 
C  2.844851597213702  -0.36438705062618526  -1.6326070715759249 
H  2.0310386301349537  -0.4917940492621542  -2.200027428377722 
C  4.049259622880897  -0.6434620776110069  -2.528921596183168 
H  4.064498458856959  -1.6007483505431985  -2.781270916637362 
H  4.885815464613099  -0.43755693059413003  -2.041329569899224 
C  3.9497190715969364  0.2300302787688986  -3.7915826718467893 
H  4.735157659608674  0.06749352019155741  -4.372271544284814 
H  3.135262296811513  -0.01063040449499586  -4.300966071200468 
C  3.8955849940836202  1.6934135608717173  -3.407854823868873 
H  3.781969956652568  2.2392983933715787  -4.2268393065766645 
H  4.754264327136393  1.9514137626235872  -2.9881595309669238 
C  2.7520296079460858  2.003439121771306  -2.437509895286869 
H  1.8849445795586932  1.8522560808077972  -2.8909242991951074 
H  2.7950360222347816  2.9564116694072187  -2.1712504132353425 
C  2.836101928934135  1.1132474257537763  -1.1772143866123677 
H  3.6967360147039363  1.3086430320633473  -0.7066184692415385 
C  1.1498334711686387  3.2047327945617416  -0.12918474466159777 
C  0.24524308437426767  3.73691664556814  -1.0583854231198346 
H  -0.2970006009106827  3.1595200033550723  -1.5815864870805916 
C  0.14571572718371617  5.107207034124412  -1.2072598876897254 
H  -0.4706266736640117  5.47176464073293  -1.8331923172546882 
C  0.9356714396439301  5.953776876366757  -0.4522714211835986 
H  0.8530451853001229  6.894578006214589  -0.5545114630955462 
C  1.8362566353920464  5.440848476667069  0.4404347913974927 
H  2.3880767237472136  6.024878562976727  0.9471137512654259 
C  1.9448589748066432  4.075001596486848  0.6032942454541647 
H  2.571910879388411  3.7248223174534365  1.2260833923538565 
C  2.0022651123077013  1.1990367685790184  1.7111859703146477 
C  3.327379143160618  1.0601568824338734  2.0903988161691442 
H  4.009064427629888  1.028062275352213  1.4309940170249984 
C  3.6562115560156663  0.9711112037682874  3.4257524777011623 
H  4.565907305315544  0.8749084363264961  3.679052993378295 
C  2.6803780807412227  1.0199621032249337  4.399556717862472 
H  2.9188770102410295  0.9714031858262262  5.318573383744298 
C  1.3439672461985885  1.139978965829026  4.029044364481682 
H  0.6658665789623552  1.1601958591872898  4.693560957439888 
C  1.0061437089353213  1.23109924090624  2.7036639537539404 
H  0.09309601244881782  1.315385208585824  2.4546674284910046 
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

