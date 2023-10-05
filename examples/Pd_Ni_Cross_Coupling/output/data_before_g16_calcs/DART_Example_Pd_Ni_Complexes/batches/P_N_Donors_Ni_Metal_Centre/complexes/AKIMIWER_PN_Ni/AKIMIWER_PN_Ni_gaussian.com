%chk=AKIMIWER_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.2932781748332416  1.5516222357584324  0.0 
N  1.2932781748332416  -1.5516222357584328  0.0 
N  2.715831284347233  -0.32193046154422933  1.404337595110974 
C  0.9040803646750701  2.5725040995876953  1.4460195623409193 
C  -0.23084989316232596  2.309520561716093  2.2149327100069702 
H  -0.8079317867137807  1.6208762542453972  1.9713965198354066 
C  -0.4986385937983293  3.069906435809536  3.3327110844134755 
H  -1.239841592339432  2.8764711149672215  3.8586664357545173 
C  0.33372847156747665  4.120017959381226  3.6730455325028286 
H  0.1333285941901965  4.651753704065701  4.409409794924635 
C  1.4541700147245058  4.377202524543726  2.9308522101062287 
H  2.017599913256587  5.078533863397663  3.167713651842729 
C  1.7470958764746298  3.599409233633524  1.830997980613506 
H  2.5212586119173084  3.7671633615400295  1.342809705844174 
C  1.7117609321332965  2.6777035361762116  -1.3496427329823195 
C  1.421890196062348  4.024896832860467  -1.2856601749004952 
H  0.9319074258769076  4.360245023939761  -0.5702541323894825 
C  1.8572887127358952  4.878589596777665  -2.2797353161820277 
H  1.6788502968306385  5.789454043953606  -2.216849621464607 
C  2.5462619636304358  4.394676197822678  -3.3525448667844007 
H  2.8566261433489153  4.979380840545552  -4.006330417807054 
C  2.782680422946303  3.0539482119394736  -3.4705750322688234 
H  3.215523352351508  2.7199348625196182  -4.222662906704549 
C  2.368649946640655  2.184051913153682  -2.453643769080931 
H  2.536955660068213  1.2733203720390538  -2.526623059613381 
C  2.878390419061641  0.735591557787969  0.3952926078341573 
H  3.249250477927242  0.3499803753679529  -0.4145253054839295 
H  3.5057383287066384  1.3980282458633928  0.7250375464696408 
C  2.2518388959910904  -1.535728808588034  0.9716747052517949 
C  2.745248380859179  -2.7562646847539893  1.4761115746755267 
H  3.340731309115645  -2.754381652316928  2.18952125104946 
C  2.3606703781213505  -3.9200334735644686  0.9284933079962874 
H  2.6640402422122405  -4.722432067312022  1.2862161508696879 
C  1.510388643990041  -3.930376870707236  -0.16794576771117994 
H  1.3039713569963003  -4.724244467221736  -0.607008807415551 
C  0.9899414678444471  -2.7401086711172384  -0.5777511360104169 
H  0.39443048728450636  -2.7425750607524257  -1.2919725435268268 
C  3.64205244616259  -0.24587679507797144  2.5594625620398705 
H  3.4346890010501125  -1.0063441709681058  3.14394320486551 
C  5.098765214903539  -0.3974391209974959  2.131058184437459 
H  5.321509214945159  0.2983047882616572  1.4952371929828947 
H  5.218381573255746  -1.2552532553440585  1.693914103235597 
C  6.047922427384682  -0.30473931837466606  3.3454958290239682 
H  5.9173831683170555  -1.0789857214708327  3.917636178182033 
H  6.966837219835858  -0.3116671161029241  3.0379601403095746 
C  5.787371647487237  0.9655105656356637  4.1464995139509435 
H  6.367012386235798  0.9803202657485048  4.923481403550442 
H  6.000784247603518  1.7377260394712297  3.5994413231882842 
C  4.353747812445057  1.0541919223659892  4.592439707886747 
H  4.213605718788114  1.8777649580897098  5.083728053376716 
H  4.149213797403961  0.31095221271056483  5.183413565753991 
C  3.4330382529409715  1.014244753801663  3.381486405417373 
H  2.510625497741448  1.0541470048606891  3.6788737108049396 
H  3.6026012317783525  1.7903869920645392  2.825853179096185 
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


