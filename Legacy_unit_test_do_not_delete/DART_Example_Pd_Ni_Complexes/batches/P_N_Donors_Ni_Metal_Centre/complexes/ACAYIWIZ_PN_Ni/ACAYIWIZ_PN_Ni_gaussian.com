%chk=ACAYIWIZ_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.2841575205168565  -1.5591790989171186  1.0618727371285693e-15 
C  1.8945650650070116  -1.9274494890209821  1.6796055666144691 
N  1.2841575205168563  1.559179098917119  5.307005967277214e-16 
C  0.9767899593053488  2.8836958727978708  0.5270962986242883 
C  2.751785096687721  -1.0486719252686143  -0.9570841934314308 
C  1.0323167653368912  -1.741298071136032  2.7697845699599917 
H  0.16201756968859593  -1.4485496560777427  2.6218145542986733 
C  0.7625737506496232  3.7573244642445043  -1.8942927160027094 
H  0.6381366448158884  2.8027410741556458  -2.078008932102924 
C  0.8324650290747467  -3.153279931673419  -0.7520289144577478 
C  -0.108604746392855  -5.498664199933412  -1.919427848798517 
H  -0.4511934815401879  -6.2749534734624275  -2.3021692005674876 
C  3.200217553826647  -2.3709047472913185  1.9218696835437434 
H  3.7844879690099247  -2.508282271880227  1.2113446187646075 
C  0.9135922794861151  3.042969724129302  1.9228906866154298 
C  3.152316085971094  0.3023959709487246  -1.0611085965036862 
C  1.073385842550871  1.8469601949513295  2.852386003074075 
H  0.6641614918189218  1.076355152955529  2.4069134449757463 
C  0.768108969392557  3.9446259296135757  -0.3743176068564653 
C  4.3052253327711085  0.6178860005044298  -1.7876335856774979 
H  4.576138154788803  1.5063254401827022  -1.8480015812594128 
C  0.42649364773353937  -5.539209907635564  -0.6408674223069297 
H  0.47380224246515656  -6.344316758972805  -0.17907877498614988 
C  0.899673705634236  -4.356636724554424  -0.04947523478460366 
H  1.2556199978799376  -4.376597436586553  0.8092072542715819 
C  3.5426683544978963  -2.023125117371573  -1.572364730191654 
H  3.3059266010312096  -2.918207620061139  -1.49309547042155 
C  -0.13852815744927938  -4.302874196530404  -2.635175567783483 
H  -0.4633292206585067  -4.29048107888701  -3.505744444997254 
C  0.3175252937555979  -3.132693845686062  -2.049624360678636 
H  0.28079992599704195  -2.3320255062585864  -2.520981441458648 
C  0.7190369490731238  4.334318983748275  2.41561582610947 
H  0.6903642432734769  4.476486009431179  3.335167626037362 
C  0.36936992558031667  2.016004059240593  4.20579945453609 
H  0.5075048132042903  1.2308556206486303  4.73951310641103 
H  -0.571568645263707  2.144847400675532  4.063189990383536 
H  0.7311587978228309  2.7804270724965945  4.661304853764338 
C  2.5487969200765868  1.509062412479834  3.077496378035562 
H  2.616978942083427  0.750838789663786  3.662254626954529 
H  2.991937197773964  2.26179914157986  3.4751814556360165 
H  2.9602776072780173  1.3026631186312267  2.2354295127813435 
C  1.4619242381791235  -1.9878919295754065  4.062856547174764 
H  0.8765489780001996  -1.8658425190027756  4.775946419595744 
C  4.671110328644721  -1.6819061352998528  -2.2987361897476895 
H  5.174966066066936  -2.34685898178174  -2.7093675911557917 
C  0.564093005431088  5.208785926607463  0.18220855465734134 
H  0.4239283134572369  5.935740470051001  -0.3818953993180588 
C  5.049947619269675  -0.3605055589126261  -2.417132912487125 
H  5.8016589045286295  -0.12976731001622857  -2.916310059510685 
C  -0.3896624420497008  4.514879206236842  -2.555709543063139 
H  -0.3625742536406362  4.372958893712332  -3.5049395915454356 
H  -0.30554899814900915  5.45234571410317  -2.371085234847576 
H  -1.2246077857211604  4.193973201694128  -2.2076986533918443 
C  3.626200821804576  -2.610457673096574  3.235652778708636 
H  4.495595107726409  -2.9036798727255477  3.3946910673412436 
C  2.7603617496410546  -2.41373325239787  4.299395860197671 
H  3.0499642577007613  -2.5671424471248394  5.169204658899314 
C  0.5644380545500756  5.4116606317774405  1.5513825154612042 
H  0.46116146666634394  6.26990540681924  1.8931699969589098 
C  2.091923541805489  4.184676992547874  -2.529278447579892 
H  2.051149719503287  4.049173104542492  -3.4793956254849716 
H  2.8064468932896496  3.6568934713529555  -2.161760999097566 
H  2.2513874039025605  5.113021985759552  -2.345069407997402 
C  2.487832349552769  1.4486425735624004  -0.43412314210130887 
H  3.015579972753986  2.2068702690343764  -0.3360924427019606 
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

