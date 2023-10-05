%chk=IJUMOYIC_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.2813461261501509  1.561490347392515  -1.6019706750515173e-15 
N  1.2813461261501513  -1.5614903473925155  -1.1102230246251565e-16 
C  1.2331026502414182  2.7401838071845805  1.3880151853908842 
C  0.4948200900362898  2.443479587543613  2.5365375740904783 
H  0.0238325867659257  1.6201913987010055  2.5963893900647537 
C  0.4461758869087905  3.3496497904928813  3.592895219392961 
H  -0.06199216240264871  3.146765683313471  4.368575110013498 
C  1.1414248308349584  4.552412995616875  3.5135526811841284 
H  1.110058601253845  5.168301967717101  4.235582973960851 
C  1.8820721460009655  4.851745047704487  2.3752597807712155 
H  2.3526944781721935  5.67517648751267  2.319464580286387 
C  1.937176087668496  3.951312792488909  1.3182037560591322 
H  2.451807829901747  4.156413781768028  0.5473265842070133 
C  1.070369629691227  2.6027525349654788  -1.4868637041345516 
C  0.1977329614247605  3.6901349692119236  -1.4751668289495892 
H  -0.26918716529912445  3.9107107163339903  -0.6772192335263548 
C  0.0042970050805895  4.4546766821937105  -2.622272974824555 
H  -0.6009956397874197  5.18569838649443  -2.60949753315007 
C  0.6953558418987954  4.148303651583009  -3.781632806874095 
H  0.5627669068664478  4.672267079449361  -4.563421448749616 
C  1.5733623463085609  3.0876033771523934  -3.8064764385927563 
H  2.053765733074028  2.8894969780212305  -4.6020094741152375 
C  1.7604806928397663  2.300950003309139  -2.6654639172712398 
H  2.3564199054358097  1.5619314260076662  -2.691393192111154 
C  3.0024577649621325  1.047391770772747  -0.15646162553001447 
C  4.211787719846178  1.9539771828208632  -0.24192372755050529 
H  4.471042201694161  2.2869282515662865  0.653366124655482 
H  4.03576253596929  2.728061330710992  -0.8345131947420253 
C  5.301405149236235  1.0484685561674236  -0.8292799181416137 
H  6.186867738468678  1.2669589653907083  -0.44342466788933754 
H  5.347871427155635  1.1452584948643725  -1.8137797113276717 
C  4.87034694293267  -0.3801659555509107  -0.43802191540467417 
H  5.083127027436514  -1.0265745111244367  -1.1579669583905123 
H  5.313390367569035  -0.6701735279679866  0.3978079846511102 
C  3.370421721229899  -0.24849678427276212  -0.2455390625712885 
C  2.5697047092083034  -1.4630109329940948  -0.13001625436912892 
H  3.049896320727584  -2.282782276529801  -0.15300213142826152 
C  0.8370207077224687  -2.9494818541086243  0.11522569081160573 
C  0.5448493769290192  -3.654957130878354  -1.0566311252410587 
C  0.18779187135735564  -5.0038324880121134  -0.9229135922387762 
H  -0.022717657764845578  -5.512990470530303  -1.6966311485774868 
C  0.13815931854852592  -5.607795014309391  0.3346997690652199 
H  -0.08072091425260997  -6.529352424128124  0.4084743859541407 
C  0.40493223069889117  -4.873440955115357  1.4711572029746804 
H  0.35621770339438075  -5.2939589372447315  2.320991375937474 
C  0.7478164195335777  -3.512593673191147  1.3901164211011947 
C  0.6040681022540414  -2.991198592933668  -2.394457966007988 
H  0.8699182721733859  -2.054485709986639  -2.2828536557197157 
H  -0.2783664627794149  -3.03128338763286  -2.817982803140548 
H  1.2590958531602683  -3.453011543749642  -2.960111193702658 
C  1.0103473641543435  -2.6957889570237974  2.6292480722034264 
H  1.2316272617634083  -1.775626185989503  2.373722453690095 
H  1.760178671258398  -3.0864705683566465  3.1256451909027634 
H  0.20996027565365183  -2.6944459548057544  3.1944777289859103 
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
