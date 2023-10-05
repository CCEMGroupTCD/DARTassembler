%chk=YEKEGEDE_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
N  1.4411345096832562  -1.415355547203598  -1.0551462493087116e-15 
P  1.4411345096832562  1.4153555472035984  -2.843533665163494e-16 
C  2.2758723606851223  0.428748268376065  -1.236709993276258 
C  2.4838130898666497  1.1426354623015158  -2.359727347325912 
C  2.973774514135397  0.5839987709329254  -3.6703838201613253 
H  3.8214182819232114  0.13753776242627613  -3.5193118478643144 
H  2.3405936646307244  -0.08193504468126216  -3.979297586513789 
C  3.153809406203804  1.5986226496023488  -4.733819291862238 
H  4.044016984733347  1.9677030520723113  -4.633399654135495 
H  3.1352943082161895  1.131520626152337  -5.583298921303377 
C  2.2830834434414817  2.651466970553486  -4.8300297050854475 
H  1.3896872727946845  2.3246695937997406  -5.016099495273185 
H  2.5531876670865015  3.235883113645067  -5.554236581649787 
C  2.286597937858599  3.493454573848314  -3.3998322836687676 
H  3.1374110586668684  3.946686248580093  -3.292555527191181 
H  1.5887349290390154  4.166815572703328  -3.41863457454626 
C  2.0623650340396424  2.5658502997180768  -2.227060385361066 
C  1.4921537837245902  2.8867769572617865  -1.053436190615591 
C  1.0674544704140811  4.207492881746806  -0.5291583348127009 
C  1.474331121935316  5.39371154483619  -1.1646314528793178 
H  1.9940340229176656  5.362996184949488  -1.9346836597962105 
C  1.088029156105298  6.6181210428057575  -0.6193125439115597 
H  1.337686834342808  7.399015240863406  -1.0601684672436185 
C  0.3516936675189444  6.725385634589745  0.5447712570227873 
H  0.12548790794389042  7.553854949114095  0.9044937601090338 
C  -0.024142046433863085  5.557539531703609  1.1404020706072149 
H  -0.48826990128778913  5.586088139210631  1.9443743222167484 
C  0.27448178427946557  4.318993522212587  0.5643024263580902 
H  -0.08220363428557853  3.5492350490483338  0.9441377922403696 
C  2.331200917929954  -1.0026438535559439  -0.9335146205803884 
C  3.3195239958705827  -1.867884765932738  -1.4251314095118943 
H  3.9065338396417584  -1.587829665443574  -2.0924114586709854 
C  3.402083173988134  -3.136849806942069  -0.9091221588150259 
H  4.046027393912407  -3.7196244684959603  -1.2399371153645284 
C  2.540029265125474  -3.5806494982924453  0.11090573402514106 
C  2.722170860327054  -4.845255125858143  0.8360998018625736 
C  3.5181960741455995  -5.946711740826808  0.3594122896421503 
C  3.8750025756643893  -6.146901043479282  -1.0021445375147133 
H  3.5892217418500927  -5.5374516781206875  -1.6444143609444108 
C  4.634210625142276  -7.219449272084647  -1.3825700031758474 
H  4.86460804388543  -7.324231985723644  -2.2793154421712 
C  5.0732660508906635  -8.17597441707371  -0.43610275592184744 
H  5.614413344764215  -8.882675826466015  -0.7043368594513254 
C  4.70280467258046  -8.060542170863354  0.8709483309585975 
H  4.99193077165513  -8.684639141935554  1.4965551160340091 
C  3.871851406144625  -6.975349900246528  1.2721892246784638 
C  3.3921758631412473  -6.943282964776021  2.6177045962346313 
H  3.6891855687938158  -7.570830771936471  3.23661507803012 
C  2.517528850144913  -6.007359621968472  2.97977287576839 
H  2.1210986985756923  -6.063208832547521  3.81964330661147 
C  2.17262054874785  -4.932264213605254  2.129217106554082 
C  1.2517435517811928  -3.9374611566759423  2.5907487353529834 
H  0.8896868259022994  -4.014152561330032  3.442745194088616 
C  0.8922671653265946  -2.869804956800517  1.7992333681938306 
H  0.19056925392623225  -2.3117643594675683  2.049967695450087 
C  1.589120610787226  -2.638694649820566  0.6223727853050446 
C  2.433553332534296  1.7430845838519955  1.4771314797906603 
C  2.164204180204564  1.0097918383997453  2.630275236486394 
H  1.5068821366926217  0.3532621615733681  2.629741795261119 
C  2.9027910038048788  1.283082454855963  3.795632672194873 
H  2.7364316014390297  0.8064987714608484  4.575313621023326 
C  3.86933934698596  2.252402639251048  3.7792119408409817 
H  4.35560090644354  2.4271722019158664  4.550006437292421 
C  4.128616419696822  2.970851323052845  2.638022557649041 
H  4.78805284002889  3.6279694173985138  2.6491287551604272 
C  3.423455528291289  2.7275060055596603  1.4710538718282318 
H  3.6068651002083927  3.2127740086116177  0.6986916872650928 
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


