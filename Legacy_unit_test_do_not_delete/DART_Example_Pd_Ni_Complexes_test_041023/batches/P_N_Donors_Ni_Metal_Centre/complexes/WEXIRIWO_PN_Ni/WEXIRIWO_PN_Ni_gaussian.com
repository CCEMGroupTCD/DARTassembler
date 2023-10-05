%chk=WEXIRIWO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.391624933306384  -1.4640628555495834  1.3192029502533475e-15 
N  1.3916249333063846  1.4640628555495834  -1.792959889799331e-16 
C  2.5386271037396213  1.0120762707580475  0.38719000269888715 
C  3.777237326619521  1.6976292893292562  0.8351165307203703 
C  4.168768460316614  3.017547906672085  0.8911828747930776 
C  5.447016617965514  3.2999088882604197  1.3629361827238178 
C  6.285071433171835  2.300148284095274  1.7659996033773395 
C  5.904244446000941  0.9664804585653113  1.7213661295204952 
C  4.644547759402419  0.6718012421883277  1.2358268136992157 
C  3.9926313392207433  -0.6313964093112622  1.0213909906111256 
C  2.7765035865967915  -0.4226094717550503  0.4949024399247893 
C  1.172071765574698  2.896296626833173  0.00853937549072483 
C  0.649518950817198  3.4574840662150566  1.1562253533039306 
C  0.44044591504038855  4.8436006610782005  1.1478892177442 
C  0.7431038345615535  5.581693230604088  0.02321570398039922 
C  1.2481771159188098  4.978415781379786  -1.0926428173906906 
C  1.4714700638680918  3.614168964770104  -1.140056452206733 
C  0.312844324736542  2.6282664199929786  2.354956574782473 
C  2.0132372536943324  2.944299432678293  -2.359675040947876 
C  1.9765975750484892  -2.551527030304804  -1.3157469578293781 
C  3.314727538435024  -2.7151792185232253  -1.5770607966870696 
C  3.742413405589189  -3.57451957760072  -2.5676870765207984 
C  2.8430028317036107  -4.2528728882085405  -3.3224519045104657 
C  1.5021162420296734  -4.114402091526864  -3.0739219825714987 
C  1.0738349074399953  -3.2615672386603225  -2.0837244304748306 
C  0.9865534812187831  -2.5002103127778104  1.4170778438685616 
C  1.3757875211573638  -3.812421781512075  1.509639301697181 
C  1.0536158714182595  -4.563788801602932  2.631786297070009 
C  0.37398409582284686  -3.997491311807815  3.6538026387460008 
C  -0.011314015317922532  -2.690577955767686  3.5820543610294786 
C  0.27896384427682874  -1.944647061571174  2.467549661482494 
H  3.597005342407355  3.699148150963102  0.6204555692231606 
H  5.735605945756057  4.182073759190893  1.4026476468323588 
H  7.133105474115329  2.5190332273101634  2.0803156020711064 
H  6.478894637892766  0.29258333983154095  2.006946941762842 
H  4.363981526586326  -1.460443277346424  1.2183752258804343 
H  0.09589973299797472  5.266691164093166  1.9030261223506255 
H  0.6002771621958684  6.500065505334928  0.023407189353226093 
H  1.4475382438205089  5.498471421145461  -1.8383101173507275 
H  -0.2935648399496902  1.927937625287569  2.099037125949394 
H  -0.09981783277962952  3.1816124337444798  3.022422014201969 
H  1.1153469787466723  2.2398721493203957  2.7111899825581047 
H  1.9279757940201416  3.533823035936438  -3.1130646909270774 
H  1.5193833107124342  2.1383556027891952  -2.52860152219024 
H  2.9384615121819015  2.7297278459399115  -2.222441387383078 
H  3.9396970987221183  -2.2408651633802  -1.079033893991683 
H  4.65261557921536  -3.688252723939915  -2.719051354175654 
H  3.135936177709176  -4.810525209558858  -4.005919762782231 
H  0.8837986954682377  -4.594976192108006  -3.5743660208663703 
H  0.16204689977567477  -3.1609545740867735  -1.9291887038386875 
H  1.8581088094142593  -4.202130766913404  0.8164898505744868 
H  1.3065961055556876  -5.4561072222554365  2.68130848704577 
H  0.16917140758461424  -4.50031041360359  4.40815568845542 
H  -0.4740454716139699  -2.3042315582370474  4.289502465464891 
H  -0.0014830220910304615  -1.0596145200457616  2.4180363745231097 
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

