%chk=YUJEBIKI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
N  1.2580152065853571  -1.5803473478953922  1.2373387775461442e-15 
P  1.2580152065853574  1.5803473478953922  -3.0455903457662574e-16 
C  2.9082963831975857  1.0417394455671105  -0.5034160363645166 
C  3.359191107732588  -0.1930004188042529  -0.03404689089508668 
C  4.645972400424913  -0.6035053314514327  -0.37157544001305576 
H  4.974142913313849  -1.4034987327019066  -0.027261045271984197 
C  5.448792299735348  0.15570516844141688  -1.2115128996632754 
H  6.296620959468298  -0.1431107283479131  -1.4481059755582777 
C  4.977569797185427  1.3598029590531795  -1.6928668169990335 
H  5.511708262437631  1.8740157054359392  -2.25321649478563 
C  3.718683337615685  1.802445458457876  -1.346787606943513 
H  3.4082803610274635  2.6122649843389136  -1.677583794386772 
C  1.5041741250460765  2.286163086116148  1.666356978724011 
C  0.3949951032424669  2.6823763495065505  2.418823441434741 
H  -0.4637114651500185  2.493955677053402  2.114029708572787 
C  0.571929651342852  3.34736346224448  3.6077571803710184 
H  -0.1692075129646795  3.632556360248971  4.092896989714694 
C  1.839552855894879  3.597220002027236  4.088594398989852 
H  1.9513007112147964  4.043775343762598  4.897389325778495 
C  2.937713663020061  3.1890001439946687  3.3718660200330253 
H  3.7906764460702327  3.3493478551086016  3.7052094949531127 
C  2.786347604956874  2.5404153051204705  2.1565929999881446 
H  3.533633919502269  2.2779004468821906  1.6706629604790344 
C  0.9180075749681409  3.042284992854778  -1.0352780702268758 
C  0.9385269445068767  4.3282248229221345  -0.5057724810967196 
H  1.104287934323049  4.453986555154407  0.40075934457799106 
C  0.7135128529968807  5.4221717215595335  -1.3228055597315878 
H  0.7100039432151386  6.281095382745598  -0.9628115433828809 
C  0.49324989844944833  5.236584445198683  -2.677822970581751 
H  0.34781227582809904  5.9705358778028375  -3.2302164829482085 
C  0.4893638428131013  3.979245276916339  -3.2024341737357274 
H  0.33714030178083343  3.861185294565959  -4.11193425296773 
C  0.709606996293909  2.8762701035651568  -2.3929441983093622 
H  0.7174076690959934  2.0219467728125755  -2.763014807836322 
C  2.4737193951759107  -1.0976462815659578  0.7831157362790375 
H  2.173902162354243  -0.6215096391992774  1.57200551609556 
H  2.9882042167060687  -1.866832949053227  1.0757755711982095 
C  1.6954012917456387  -2.0304416331468413  -1.376527241609183 
H  2.34146970053408  -1.3804315000150236  -1.7260110437968479 
C  2.413495342604647  -3.3532548253124776  -1.1418892379664025 
H  2.266244897509717  -3.959428184147646  -1.8833423891889565 
H  3.36725554915217  -3.213579822665775  -1.0377802467980315 
C  1.8084014419429888  -3.909604194437498  0.13666683857318243 
H  2.4981905433848994  -4.0483256598747985  0.8046948773252344 
H  1.3713935956419365  -4.7571146621733895  -0.03740015815339997 
C  0.5386301685499888  -2.1417192774154867  -2.3613900109465074 
H  -0.13541192675642155  -2.7198942534527686  -1.998408687956927 
H  0.8585579545614707  -2.501690343749478  -3.1911914589689694 
H  0.16337546565398453  -1.2710287600222927  -2.514362364252653 
C  0.5452676575928177  -2.779201303972275  2.121218124899967 
H  -0.04016784612671764  -2.0421222361640092  2.305186733047001 
H  1.3817418061283764  -2.644099523794432  2.5722496956624976 
H  0.14165890493664834  -3.5927330580592263  2.4320942254151876 
C  0.7855210461081249  -2.872679443570342  0.631706793380698 
H  -0.0737785563379949  -3.098189879914992  0.21648979428792833 
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


