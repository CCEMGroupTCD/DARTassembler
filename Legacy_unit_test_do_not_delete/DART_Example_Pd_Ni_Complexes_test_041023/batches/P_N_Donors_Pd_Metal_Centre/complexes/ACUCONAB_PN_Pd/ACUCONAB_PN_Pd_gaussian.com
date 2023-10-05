%chk=ACUCONAB_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5133013967812217  -1.4839201065084333  9.015945847919875e-16 
N  1.5133013967812228  1.4839201065084329  -1.236439674256494e-15 
C  2.0211904484530514  -1.506746234311025  1.7426458673756682 
C  1.0936581109081405  -1.2327803023063892  2.718959157723256 
H  0.22740654376187797  -0.9992705847143045  2.4729308174721116 
C  1.424718920016291  -1.2964833174386678  4.0439005928575735 
H  0.7847339837139771  -1.1127198292266076  4.692554500350803 
C  2.6877645935083345  -1.6290488193383454  4.411067072290554 
H  2.9091455235251367  -1.688298471401402  5.312350447570125 
C  3.6404902733804003  -1.8774630691239065  3.455355651716459 
H  4.511480370428134  -2.0857293926574383  3.7059180178582674 
C  3.3022813543857255  -1.817662743883489  2.1298912092610522 
H  3.947329745936689  -1.9875511488550868  1.4834852977675552 
C  1.1737827083832244  -3.2150437254506503  -0.42669471073462734 
C  1.0682162628281635  -4.171538784567859  0.5164740214958577 
H  1.2514087566237047  -3.9698179077584443  1.4061860356196532 
C  0.6791609162062965  -5.478376256136938  0.15252454674987082 
H  0.5778055572012155  -6.134260377344825  0.8031685436854606 
C  0.45669441630422414  -5.770150825469213  -1.139654975264121 
H  0.22961688666020796  -6.63948254308798  -1.3809379424888701 
C  0.5575423767086308  -4.830033784207147  -2.089159039294897 
H  0.40394572533400375  -5.053699393658954  -2.9793124383419767 
C  0.8896821367175491  -3.5317227386528556  -1.7535172677393516 
H  0.9223887873014839  -2.8734844172527354  -2.4100487576527727 
C  3.007115162030858  -1.0327935114418452  -0.9278878649785179 
C  3.7869177413671395  -2.0017073253957705  -1.555573875966302 
H  3.5473422209855134  -2.8973974061288645  -1.4865638092353117 
C  4.909679770273957  -1.652151825254664  -2.2774243374652228 
H  5.398065801209728  -2.310779918312616  -2.714427426495061 
C  5.307976644779242  -0.3626699410160908  -2.3601722074066585 
H  6.064940145209841  -0.1327368344401402  -2.850329575520398 
C  4.571105813722927  0.6154842187892458  -1.6999877195596025 
H  4.859194734752325  1.5004432164066395  -1.71882454827716 
C  3.420818798907992  0.2914491862487206  -1.0212855011790893 
C  2.7341452892915346  1.4002507469145817  -0.33466211077063196 
H  3.2642548773122604  2.1330193779418334  -0.11705581933825525 
C  1.1224663725979263  2.7059316029809324  0.7981528527548482 
H  0.1509290838636148  2.663702216961326  0.9277179543648525 
C  1.4111266495065125  3.9815054978948234  0.022131736322259028 
H  1.1748767379231477  4.741125732876521  0.5577229310624415 
H  0.8954074364720019  3.9871619547102446  -0.7879865438093416 
H  2.3455807218929525  4.018258915252533  -0.1943424490046432 
C  1.7351636069051937  2.688867887832832  2.14040410296757 
H  1.403219454173191  1.9364405367221647  2.6345537774262544 
H  1.51230293169261  3.4997914070453278  2.6037428833775595 
H  2.689118950782505  2.62055689888174  2.0561203520041937 
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

-Pd 0
lanl2dz


