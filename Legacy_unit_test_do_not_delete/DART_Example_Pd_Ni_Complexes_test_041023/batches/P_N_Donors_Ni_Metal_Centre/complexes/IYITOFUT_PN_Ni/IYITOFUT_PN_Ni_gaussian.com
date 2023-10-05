%chk=IYITOFUT_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.391624933306384  1.4640628555495834  -1.1399069612734144e-15 
N  1.3916249333063846  -1.4640628555495834  0.0 
C  2.5386271037396213  -1.0120762707580475  -0.38719000269888726 
C  3.777237326619521  -1.6976292893292562  -0.8351165307203705 
C  4.168768460316614  -3.017547906672085  -0.8911828747930779 
C  5.447016617965514  -3.2999088882604197  -1.3629361827238182 
C  6.285071433171835  -2.300148284095274  -1.7659996033773397 
C  5.904244446000941  -0.9664804585653111  -1.7213661295204954 
C  4.644547759402419  -0.6718012421883276  -1.2358268136992157 
C  3.9926313392207433  0.6313964093112623  -1.0213909906111256 
C  2.7765035865967915  0.42260947175505037  -0.4949024399247892 
C  1.172071765574698  -2.896296626833173  -0.008539375490725185 
C  0.649518950817198  -3.4574840662150566  -1.156225353303931 
C  0.44044591504038855  -4.8436006610782005  -1.1478892177442006 
C  0.7431038345615535  -5.581693230604088  -0.023215703980399904 
C  1.2481771159188098  -4.978415781379786  1.09264281739069 
C  1.4714700638680918  -3.614168964770104  1.1400564522067325 
C  0.312844324736542  -2.628266419992978  -2.3549565747824737 
C  2.0132372536943324  -2.9442994326782936  2.3596750409478755 
C  1.9765975750484892  2.551527030304804  1.3157469578293783 
C  3.314727538435024  2.7151792185232253  1.5770607966870698 
C  3.742413405589189  3.5745195776007197  2.567687076520799 
C  2.8430028317036107  4.2528728882085405  3.322451904510466 
C  1.5021162420296734  4.114402091526864  3.073921982571499 
C  1.0738349074399953  3.261567238660322  2.083724430474831 
C  0.9865534812187831  2.5002103127778104  -1.4170778438685614 
C  1.3757875211573638  3.812421781512075  -1.5096393016971805 
C  1.0536158714182595  4.563788801602932  -2.6317862970700086 
C  0.37398409582284686  3.9974913118078153  -3.6538026387460003 
C  -0.011314015317922532  2.6905779557676865  -3.582054361029478 
C  0.27896384427682874  1.9446470615711742  -2.4675496614824937 
H  3.597005342407355  -3.699148150963102  -0.620455569223161 
H  5.735605945756057  -4.182073759190893  -1.4026476468323592 
H  7.133105474115329  -2.519033227310163  -2.080315602071107 
H  6.478894637892766  -0.2925833398315407  -2.006946941762842 
H  4.363981526586326  1.4604432773464242  -1.2183752258804341 
H  0.09589973299797472  -5.266691164093166  -1.9030261223506262 
H  0.6002771621958684  -6.500065505334928  -0.023407189353226887 
H  1.4475382438205089  -5.498471421145461  1.8383101173507268 
H  -0.2935648399496902  -1.9279376252875688  -2.0990371259493945 
H  -0.09981783277962952  -3.1816124337444793  -3.0224220142019695 
H  1.1153469787466723  -2.2398721493203952  -2.711189982558105 
H  1.9279757940201416  -3.5338230359364387  3.113064690927077 
H  1.5193833107124342  -2.1383556027891957  2.5286015221902396 
H  2.9384615121819015  -2.729727845939912  2.2224413873830775 
H  3.9396970987221183  2.2408651633802  1.0790338939916833 
H  4.65261557921536  3.6882527239399145  2.7190513541756545 
H  3.135936177709176  4.8105252095588575  4.005919762782232 
H  0.8837986954682377  4.594976192108006  3.574366020866371 
H  0.16204689977567477  3.160954574086773  1.929188703838688 
H  1.8581088094142593  4.202130766913404  -0.8164898505744862 
H  1.3065961055556876  5.4561072222554365  -2.681308487045769 
H  0.16917140758461424  4.500310413603591  -4.408155688455419 
H  -0.4740454716139699  2.304231558237048  -4.289502465464891 
H  -0.0014830220910304615  1.0596145200457618  -2.4180363745231097 
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


