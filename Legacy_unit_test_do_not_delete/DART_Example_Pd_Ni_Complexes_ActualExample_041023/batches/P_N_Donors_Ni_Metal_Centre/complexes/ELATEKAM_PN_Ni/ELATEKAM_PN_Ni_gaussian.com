%chk=ELATEKAM_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.368839441826541  -1.4853883608336238  2.0037617829716163e-15 
N  1.3688394418265428  1.4853883608336238  -1.5374156464789625e-17 
C  1.289638534486775  2.446163079843188  -1.2092103727593018 
C  0.6743484987117164  3.7580244161649574  -0.722242092334382 
C  0.5544468573751127  3.566517711588144  0.7613187990708097 
C  0.10696903573055905  4.452866251616506  1.774568583046045 
C  0.18683535889444025  4.0220012523507105  3.121918508949496 
C  0.7083307238721326  2.772695273570487  3.4847143822595203 
C  1.0676045967123753  1.8553658694035569  2.4515069062194255 
C  0.9948945649482728  2.3121192565239737  1.1676398824841323 
C  1.2529274952402538  -2.5698957007833787  1.4427948815523792 
C  1.9698338876714274  -2.291999400472838  2.6196572994947016 
C  1.8261976141500338  -3.11740666925736  3.7184209690645593 
C  0.9843720085551657  -4.219553557685397  3.673881030882532 
C  0.28035349197909665  -4.487256035693799  2.543171496304265 
C  0.4032486293162022  -3.6671613102198855  1.4304406454383356 
C  1.628723081776299  -2.554327078724821  -1.4338776241644995 
C  2.381280884923372  -3.716938517951946  -1.325162793986181 
C  2.6332797274429143  -4.496199905425312  -2.4924020596733505 
C  2.121216324756265  -4.090412081921887  -3.7042327260312273 
C  1.4108469098405227  -2.946341982270455  -3.811562411746913 
C  1.139858140843623  -2.158243252761747  -2.6903035122842103 
C  2.8500489907214375  -0.4909548685912193  0.19133424097425666 
C  4.130190542166848  -1.1031232778614277  0.33106615198450823 
C  5.23547672998307  -0.3174260549849643  0.42728442980583314 
C  5.098588398445027  1.0648714566554824  0.42942545791337233 
C  3.885407944751689  1.716979902191159  0.31447145825478756 
C  2.7359498834297393  0.8952794720659267  0.16251347592405413 
C  3.902524243671026  3.2236257760985216  0.32432148105799596 
H  0.7449158182645341  2.0666554205239547  -1.8888701324549293 
H  2.160212676881891  2.6053929528087068  -1.5538276824644477 
H  -0.17767668507242895  3.8990594818244366  -1.1163033140398815 
H  1.2431724462350988  4.49256610827009  -0.9238518125446299 
H  -0.23101820048478028  5.315485179289331  1.556805882438027 
H  -0.13268766312582803  4.596102051498899  3.8044459961735835 
H  0.823009781032598  2.5473904483649226  4.402887306924544 
H  1.340836074659705  0.967419151375024  2.65219550928666 
H  2.5495380662148883  -1.5395836218997732  2.6604411386169895 
H  2.3078864791356692  -2.9301515665959257  4.514441624497643 
H  0.8975747884390028  -4.783303885263057  4.435916030561823 
H  -0.2981730532572686  -5.241429920139281  2.5147039976910457 
H  -0.09883717692408145  -3.86212056305628  0.6469160236774496 
H  2.726715200949757  -3.988709558202621  -0.4808582073632858 
H  3.1586108684112513  -5.2864284180152765  -2.4357937544631563 
H  2.2664487812768357  -4.624753100693948  -4.475474691454955 
H  1.0908406745184651  -2.6734400189263368  -4.665079692039437 
H  0.6269569182810724  -1.3635962862690496  -2.7771715992770933 
H  4.209509794560972  -2.0511631354039292  0.3630770049945482 
H  6.098742298000803  -0.7107610207121139  0.491061101179148 
H  5.88401708450088  1.5941307606389858  0.5108365644497836 
H  4.42898670274306  3.5365126528606154  -0.40117741731326634 
H  4.277772467454442  3.5248699333833375  1.1430609684709134 
H  3.016449821455118  3.547745159752374  0.24067625743693094 
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


