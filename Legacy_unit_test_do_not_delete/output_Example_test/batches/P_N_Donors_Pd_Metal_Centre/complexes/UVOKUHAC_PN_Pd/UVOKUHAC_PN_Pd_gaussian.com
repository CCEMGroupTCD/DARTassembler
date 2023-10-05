%chk=UVOKUHAC_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  3.114063430460039  1.0189966883176338  -0.36724486404533685 
H  3.7356537023036482  1.781392396802081  -0.18549234650047552 
H  3.0676209658321794  0.8029923386616297  -1.2938229872540572 
C  3.5582827951483873  -0.2171750517881596  0.4130326026289888 
H  4.124453897639425  -0.7512810554380751  -0.166510976033755 
H  4.105845259581552  0.07756466876402027  1.1571418047443316 
C  2.4678732561107832  -1.1205847283103665  0.9637471870934574 
H  2.8902454109352926  -1.9091268299944726  1.340863623622708 
H  2.0272493551959414  -0.6536976597570214  1.6903550202195996 
N  1.4223457201046448  -1.5713155801747787  1.733840950213971e-16 
C  2.051973004796302  -2.2792438054369937  -1.1298832993937669 
H  2.6083215215799385  -1.6697121288698542  -1.62045675274597 
H  1.3699913217855384  -2.627159689997365  -1.7093668982198182 
H  2.5885854797695176  -3.002893627501077  -0.7955722515210578 
C  0.5787250231904559  -2.529236709444416  0.7606952476759882 
H  -0.13199249679610503  -2.8456329453582105  0.1985157935123242 
H  0.2068344757563747  -2.0902278854351284  1.5288216260612029 
H  1.1163105927283834  -3.271596649559216  1.0461532038400183 
P  1.4223457201046443  1.5713155801747787  2.9613945353890404e-17 
C  1.2832209760081126  3.0810869844245374  -1.003077615386668 
C  1.0880591116128338  4.348580968528272  -0.4483997837162383 
H  1.0051248527847574  4.448918591338875  0.47274126963084506 
C  1.0198964582191985  5.455275468309516  -1.2786189576915805 
H  0.8747388344052248  6.2978865241056585  -0.9135850374885536 
C  1.1637852667989153  5.316760177208162  -2.637312600459017 
H  1.1333854922928805  6.067477763721746  -3.184497550004618 
C  1.3522867881396756  4.064760843133478  -3.196099647489261 
H  1.4430986730363176  3.9709837770500305  -4.116905883044986 
C  1.4045678388490237  2.9600690759921338  -2.3808131142097806 
H  1.5220775949244933  2.1190754578590965  -2.7579184595532067 
C  1.3809339384142922  2.0304940094497845  1.7468020306908882 
C  2.4717280100052976  2.6456946427258052  2.353500773230158 
H  3.2186854188408764  2.878126395842754  1.849576364683967 
C  2.436847656315841  2.9074221882225895  3.71189923938609 
H  3.165350666667371  3.3126370145970445  4.123618219268099 
C  1.3317409070928696  2.5716625059178306  4.456448021240208 
H  1.308817596538106  2.7634910404293214  5.367060581652399 
C  0.2673815097969221  1.9595762661231357  3.8626379771673602 
H  -0.47711515025443596  1.7284924631103689  4.370988190266133 
C  0.28985925952215164  1.681400652952304  2.518460850072729 
H  -0.43526256196336144  1.2557081673437807  2.1238828454881062 
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
