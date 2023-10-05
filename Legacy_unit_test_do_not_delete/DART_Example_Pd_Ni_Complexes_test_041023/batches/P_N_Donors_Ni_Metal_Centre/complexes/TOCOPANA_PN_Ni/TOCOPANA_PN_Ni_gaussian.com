%chk=TOCOPANA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.2813461261501509  -1.561490347392515  1.7931980906348906e-15 
N  1.2813461261501513  1.5614903473925155  -8.020511312085761e-17 
C  1.2331026502414182  -2.7401838071845805  -1.3880151853908838 
C  0.4948200900362898  -2.4434795875436133  -2.536537574090478 
H  0.0238325867659257  -1.6201913987010057  -2.5963893900647537 
C  0.4461758869087905  -3.3496497904928817  -3.5928952193929606 
H  -0.06199216240264871  -3.1467656833134714  -4.368575110013498 
C  1.1414248308349584  -4.552412995616875  -3.513552681184128 
H  1.110058601253845  -5.168301967717102  -4.23558297396085 
C  1.8820721460009655  -4.851745047704487  -2.375259780771215 
H  2.3526944781721935  -5.67517648751267  -2.319464580286386 
C  1.937176087668496  -3.951312792488909  -1.3182037560591318 
H  2.451807829901747  -4.156413781768028  -0.5473265842070127 
C  1.070369629691227  -2.6027525349654788  1.4868637041345518 
C  0.1977329614247605  -3.6901349692119236  1.4751668289495896 
H  -0.26918716529912445  -3.9107107163339903  0.6772192335263553 
C  0.0042970050805895  -4.4546766821937105  2.6222729748245555 
H  -0.6009956397874197  -5.18569838649443  2.6094975331500705 
C  0.6953558418987954  -4.148303651583008  3.7816328068740956 
H  0.5627669068664478  -4.67226707944936  4.563421448749617 
C  1.5733623463085609  -3.087603377152393  3.8064764385927568 
H  2.053765733074028  -2.88949697802123  4.6020094741152375 
C  1.7604806928397663  -2.3009500033091386  2.66546391727124 
H  2.3564199054358097  -1.561931426007666  2.691393192111154 
C  3.0024577649621325  -1.047391770772747  0.1564616255300146 
C  4.211787719846178  -1.9539771828208632  0.24192372755050554 
H  4.471042201694161  -2.2869282515662865  -0.6533661246554817 
H  4.03576253596929  -2.728061330710992  0.8345131947420257 
C  5.301405149236235  -1.0484685561674236  0.8292799181416138 
H  6.186867738468678  -1.2669589653907083  0.4434246678893377 
H  5.347871427155635  -1.1452584948643723  1.813779711327672 
C  4.87034694293267  0.38016595555091076  0.4380219154046741 
H  5.083127027436514  1.0265745111244369  1.1579669583905121 
H  5.313390367569035  0.6701735279679866  -0.3978079846511103 
C  3.370421721229899  0.24849678427276214  0.24553906257128846 
C  2.5697047092083034  1.4630109329940948  0.13001625436912875 
H  3.049896320727584  2.282782276529801  0.15300213142826125 
C  0.8370207077224687  2.9494818541086243  -0.1152256908116061 
C  0.5448493769290192  3.654957130878354  1.0566311252410583 
C  0.18779187135735564  5.0038324880121134  0.9229135922387756 
H  -0.022717657764845578  5.512990470530303  1.696631148577486 
C  0.13815931854852592  5.607795014309391  -0.3346997690652206 
H  -0.08072091425260997  6.529352424128124  -0.4084743859541415 
C  0.40493223069889117  4.873440955115357  -1.471157202974681 
H  0.35621770339438075  5.2939589372447315  -2.3209913759374743 
C  0.7478164195335777  3.512593673191147  -1.3901164211011952 
C  0.6040681022540414  2.9911985929336686  2.3944579660079874 
H  0.8699182721733859  2.0544857099866394  2.282853655719715 
H  -0.2783664627794149  3.0312833876328606  2.8179828031405476 
H  1.2590958531602683  3.4530115437496423  2.9601111937026574 
C  1.0103473641543435  2.695788957023797  -2.629248072203427 
H  1.2316272617634083  1.7756261859895028  -2.373722453690095 
H  1.760178671258398  3.086470568356646  -3.125645190902764 
H  0.20996027565365183  2.694445954805754  -3.1944777289859108 
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


