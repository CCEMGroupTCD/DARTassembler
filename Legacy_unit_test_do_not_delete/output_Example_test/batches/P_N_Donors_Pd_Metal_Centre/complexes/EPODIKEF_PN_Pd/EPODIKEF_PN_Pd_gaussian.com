%chk=EPODIKEF_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5374627507682908  -1.458872266512734  9.482120580096076e-17 
N  1.5374627507682908  1.458872266512734  -6.763802269245076e-17 
H  1.840751666755772  1.1875866720953105  0.7956001129340912 
C  2.5549019479874193  -0.4419790236614062  -1.126206443155724 
C  3.388846350591275  -0.9789124193948318  -2.087172281212375 
H  3.4638071349252  -1.9293330023526332  -2.1703624627375793 
C  4.096945703124574  -0.14901763162560258  -2.9428662324445884 
H  4.664800503846504  -0.5048706572238995  -3.62458717203029 
C  3.9656980013557095  1.1960045440926677  -2.845440324152772 
H  4.4755809135151114  1.7588405845383197  -3.426346433143997 
C  3.123393924638438  1.7612964815227519  -1.8863343144635338 
H  3.0362149664945175  2.7115478104576565  -1.8231835336734783 
C  2.412465470974051  0.9401555855733497  -1.0134140587873224 
C  1.0924881597704146  2.855016313358404  -0.1537928276080449 
H  0.5155557734542386  3.0898373993873305  0.5770866010091111 
H  0.6158286736320883  2.95063598699402  -0.9813732420018295 
H  1.856695927405354  3.4366198364880916  -0.15715879991519557 
C  2.6455853618561136  -1.975835113303144  1.356105590960501 
C  3.996101531181071  -2.274271727661939  1.1355840418956842 
H  4.375498940438323  -2.1633638204283834  0.26336259526531 
C  4.8040400714035085  -2.658735996998369  2.1863159919318638 
H  5.71159646005809  -2.9163507578152097  2.0271106640639096 
C  4.29457771088018  -2.7553828338986492  3.4541297290901434 
H  4.877724571861357  -2.980230880934395  4.178658716965742 
C  2.975295407446681  -2.4797993100847147  3.6841923406756223 
H  2.6159890204580822  -2.5600028362506424  4.564988273147929 
C  2.1535306465195907  -2.064089968066932  2.635964597931878 
H  1.2296107655162838  -1.8713160107743887  2.8010322049194696 
C  1.2077597273817833  -2.988682463902294  -0.9293902981701145 
C  1.4511483962380725  -4.240973752136235  -0.4055487712884146 
H  1.8431170368758272  -4.323361424284743  0.4640452697643969 
C  1.1007304971180283  -5.382555361594752  -1.1120657030900702 
H  1.3027497391414418  -6.246065388426157  -0.7538660528753786 
C  0.5168999718629201  -5.291917796514379  -2.352297298227311 
H  0.29082682860202724  -6.088423005827616  -2.835577584844735 
C  0.27446183530736246  -4.048312299552822  -2.887850935009357 
H  -0.12734143561977218  -3.991978142076204  -3.7547702012468043 
C  0.592111335546919  -2.8941477927280546  -2.1894230531321255 
H  0.4018987466666981  -2.0364463621627635  -2.5674487758518167 
Br  -1.654629867976521  1.654629867976521  0.0 
C  -1.352040782177733  -1.3199206589279684  -0.4819701001757076 
C  -1.6380700204479186  -1.599217458033982  -1.8139382277462257 
H  -1.1579680913617836  -1.159339363685303  -2.503515938957912 
C  -2.63068216563725  -2.52632649780156  -2.142289956238523 
H  -2.8368222554424674  -2.694399528648736  -3.0538898079700485 
C  -3.3101544432602585  -3.1961573187854975  -1.1619874135498987 
H  -3.9960123560788676  -3.811990145619169  -1.3868663427944008 
C  -2.9974710358883443  -2.9701089012107795  0.12423560353005347 
H  -3.438659168087169  -3.468581175155742  0.8017084178296768 
C  -2.0389854250767776  -2.015940632791841  0.4912562833419403 
H  -1.8619130127668426  -1.8492165108145995  1.4101568590026907 

-Pd 0
lanl2dz
F 1 1.0
1.472 1.0
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

-Pd 0
lanl2dz
