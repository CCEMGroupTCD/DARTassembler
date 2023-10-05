%chk=FUSEPAKO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.2841575205168547  1.559179098917119  -2.980939406545539e-15 
C  1.8945650650070154  1.9274494890209828  -1.6796055666144694 
N  1.2841575205168563  -1.559179098917119  -7.216449660063518e-16 
C  0.9767899593053483  -2.8836958727978708  -0.5270962986242897 
C  2.7517850966877244  1.0486719252686145  0.9570841934314354 
C  1.0323167653368952  1.7412980711360329  -2.7697845699600014 
H  0.1620175696885986  1.4485496560777433  -2.6218145542986724 
C  0.7625737506496202  -3.7573244642445056  1.8942927160027034 
H  0.6381366448158902  -2.802741074155646  2.0780089321029163 
C  0.832465029074746  3.1532799316734197  0.7520289144577457 
C  -0.108604746392861  5.498664199933413  1.9194278487985212 
H  -0.451193481540189  6.274953473462428  2.302169200567487 
C  3.2002175538266444  2.370904747291319  -1.9218696835437437 
H  3.784487969009931  2.508282271880228  -1.2113446187646042 
C  0.9135922794861148  -3.0429697241293017  -1.9228906866154312 
C  3.1523160859710955  -0.3023959709487245  1.0611085965036877 
C  1.0733858425508698  -1.846960194951329  -2.8523860030740824 
H  0.6641614918189287  -1.0763551529555289  -2.406913444975738 
C  0.7681089693925555  -3.9446259296135766  0.37431760685646254 
C  4.305225332771112  -0.6178860005044299  1.7876335856775045 
H  4.576138154788794  -1.5063254401827024  1.8480015812594148 
C  0.4264936477335374  5.5392099076355645  0.6408674223069273 
H  0.4738022424651548  6.344316758972806  0.1790787749861506 
C  0.8996737056342351  4.356636724554425  0.04947523478460222 
H  1.2556199978799394  4.376597436586554  -0.8092072542715838 
C  3.542668354497896  2.023125117371573  1.5723647301916541 
H  3.3059266010312105  2.91820762006114  1.4930954704215555 
C  -0.13852815744928737  4.3028741965304045  2.6351755677834863 
H  -0.46332922065851  4.2904810788870105  3.505744444997259 
C  0.3175252937555999  3.1326938456860622  2.0496243606786386 
H  0.28079992599703774  2.332025506258586  2.5209814414586518 
C  0.7190369490731248  -4.334318983748276  -2.41561582610948 
H  0.6903642432734789  -4.47648600943118  -3.3351676260373635 
C  0.36936992558032467  -2.0160040592405926  -4.205799454536089 
H  0.5075048132042927  -1.2308556206486296  -4.73951310641102 
H  -0.5715686452637119  -2.1448474006755314  -4.06318999038355 
H  0.7311587978228395  -2.780427072496594  -4.661304853764321 
C  2.548796920076592  -1.5090624124798335  -3.0774963780355713 
H  2.6169789420834197  -0.7508387896637855  -3.662254626954518 
H  2.9919371977739653  -2.2617991415798597  -3.4751814556360237 
H  2.960277607278025  -1.3026631186312265  -2.2354295127813573 
C  1.461924238179128  1.9878919295754074  -4.062856547174757 
H  0.8765489780002036  1.8658425190027768  -4.775946419595722 
C  4.671110328644721  1.681906135299853  2.298736189747691 
H  5.174966066066903  2.3468589817817396  2.7093675911557678 
C  0.564093005431088  -5.208785926607464  -0.18220855465734292 
H  0.4239283134572356  -5.935740470051002  0.38189539931805966 
C  5.049947619269678  0.36050555891262603  2.4171329124871344 
H  5.801658904528626  0.12976731001622843  2.9163100595106903 
C  -0.3896624420497057  -4.514879206236843  2.5557095430631467 
H  -0.36257425364063667  -4.372958893712333  3.504939591545428 
H  -0.3055489981490087  -5.452345714103171  2.3710852348475706 
H  -1.2246077857211621  -4.193973201694129  2.207698653391843 
C  3.626200821804582  2.610457673096575  -3.2356527787086344 
H  4.4955951077264125  2.903679872725549  -3.394691067341233 
C  2.7603617496410555  2.4137332523978703  -4.29939586019766 
H  3.049964257700777  2.5671424471248407  -5.169204658899296 
C  0.5644380545500756  -5.411660631777441  -1.5513825154612044 
H  0.46116146666634406  -6.269905406819241  -1.8931699969589184 
C  2.09192354180549  -4.184676992547875  2.529278447579893 
H  2.051149719503286  -4.049173104542493  3.479395625484969 
H  2.8064468932896496  -3.656893471352957  2.161760999097571 
H  2.251387403902556  -5.113021985759553  2.345069407997394 
C  2.4878323495527654  -1.4486425735624004  0.4341231421013091 
H  3.0155799727539887  -2.2068702690343764  0.3360924427019661 
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
