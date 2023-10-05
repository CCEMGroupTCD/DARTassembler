%chk=UQOKEWAH_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.368839441826541  1.4853883608336238  -1.821854172813053e-15 
N  1.3688394418265428  -1.4853883608336238  -1.6653345369377348e-16 
C  1.289638534486775  -2.446163079843188  1.2092103727593015 
C  0.6743484987117164  -3.7580244161649574  0.7222420923343815 
C  0.5544468573751127  -3.566517711588144  -0.7613187990708101 
C  0.10696903573055905  -4.452866251616506  -1.7745685830460454 
C  0.18683535889444025  -4.0220012523507105  -3.1219185089494963 
C  0.7083307238721326  -2.7726952735704864  -3.484714382259521 
C  1.0676045967123753  -1.8553658694035566  -2.451506906219426 
C  0.9948945649482728  -2.3121192565239737  -1.1676398824841325 
C  1.2529274952402538  2.5698957007833787  -1.442794881552379 
C  1.9698338876714274  2.2919994004728386  -2.619657299494701 
C  1.8261976141500338  3.1174066692573605  -3.718420969064559 
C  0.9843720085551657  4.219553557685398  -3.6738810308825314 
C  0.28035349197909665  4.487256035693799  -2.5431714963042644 
C  0.4032486293162022  3.6671613102198855  -1.4304406454383352 
C  1.628723081776299  2.554327078724821  1.4338776241644997 
C  2.381280884923372  3.716938517951946  1.3251627939861814 
C  2.6332797274429143  4.496199905425312  2.492402059673351 
C  2.121216324756265  4.090412081921886  3.7042327260312278 
C  1.4108469098405227  2.9463419822704546  3.8115624117469133 
C  1.139858140843623  2.1582432527617468  2.6903035122842107 
C  2.8500489907214375  0.4909548685912193  -0.1913342409742566 
C  4.130190542166848  1.1031232778614277  -0.3310661519845081 
C  5.23547672998307  0.31742605498496435  -0.4272844298058331 
C  5.098588398445027  -1.0648714566554824  -0.42942545791337244 
C  3.885407944751689  -1.716979902191159  -0.3144714582547878 
C  2.7359498834297393  -0.8952794720659267  -0.16251347592405424 
C  3.902524243671026  -3.2236257760985216  -0.32432148105799635 
H  0.7449158182645341  -2.066655420523955  1.888870132454929 
H  2.160212676881891  -2.6053929528087068  1.5538276824644475 
H  -0.17767668507242895  -3.8990594818244366  1.116303314039881 
H  1.2431724462350988  -4.49256610827009  0.9238518125446293 
H  -0.23101820048478028  -5.315485179289331  -1.5568058824380278 
H  -0.13268766312582803  -4.596102051498898  -3.804445996173584 
H  0.823009781032598  -2.547390448364922  -4.402887306924544 
H  1.340836074659705  -0.9674191513750237  -2.65219550928666 
H  2.5495380662148883  1.5395836218997734  -2.6604411386169895 
H  2.3078864791356692  2.930151566595926  -4.514441624497643 
H  0.8975747884390028  4.783303885263058  -4.435916030561822 
H  -0.2981730532572686  5.241429920139281  -2.5147039976910452 
H  -0.09883717692408145  3.86212056305628  -0.6469160236774492 
H  2.726715200949757  3.988709558202621  0.48085820736328627 
H  3.1586108684112513  5.2864284180152765  2.4357937544631567 
H  2.2664487812768357  4.6247531006939475  4.475474691454956 
H  1.0908406745184651  2.6734400189263363  4.665079692039437 
H  0.6269569182810724  1.3635962862690492  2.7771715992770933 
H  4.209509794560972  2.0511631354039292  -0.3630770049945479 
H  6.098742298000803  0.710761020712114  -0.49106110117914786 
H  5.88401708450088  -1.5941307606389858  -0.5108365644497839 
H  4.42898670274306  -3.5365126528606154  0.4011774173132659 
H  4.277772467454442  -3.5248699333833375  -1.1430609684709139 
H  3.016449821455118  -3.547745159752374  -0.24067625743693138 
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


