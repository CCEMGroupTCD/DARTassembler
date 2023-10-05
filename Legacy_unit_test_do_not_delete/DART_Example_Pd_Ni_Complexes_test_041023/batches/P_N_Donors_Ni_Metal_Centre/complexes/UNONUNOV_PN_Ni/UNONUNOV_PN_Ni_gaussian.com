%chk=UNONUNOV_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.3704309969859845  1.4839201065084335  3.9352402058811553e-16 
N  1.3704309969859847  -1.4839201065084335  1.0547118733938987e-15 
C  1.8783200486578133  1.5067462343110263  -1.7426458673756693 
C  0.9507877111129037  1.2327803023063897  -2.7189591577232575 
H  0.08453614396663611  0.9992705847143055  -2.472930817472113 
C  1.281848520221051  1.296483317438669  -4.0439005928575815 
H  0.6418635839187359  1.112719829226609  -4.6925545003508615 
C  2.544894193713122  1.6290488193383466  -4.4110670722906224 
H  2.766275123729907  1.6882984714014038  -5.312350447570191 
C  3.497619873585175  1.8774630691239071  -3.4553556517164807 
H  4.3686099706329244  2.0857293926574396  -3.7059180178582913 
C  3.159410954590488  1.8176627438834903  -2.1298912092610545 
H  3.8044593461414524  1.9875511488550877  -1.4834852977675572 
C  1.0309123085879863  3.215043725450651  0.42669471073462767 
C  0.9253458630329263  4.17153878456786  -0.5164740214958574 
H  1.1085383568284666  3.969817907758446  -1.4061860356196527 
C  0.5362905164110574  5.47837625613694  -0.1525245467498706 
H  0.4349351574059773  6.134260377344827  -0.8031685436854619 
C  0.3138240165089863  5.770150825469216  1.1396549752641219 
H  0.08674648686496855  6.639482543087983  1.3809379424888737 
C  0.41467197691339275  4.830033784207149  2.0891590392948975 
H  0.26107532553876567  5.053699393658956  2.979312438341977 
C  0.7468117369223094  3.5317227386528565  1.7535172677393487 
H  0.7795183875062465  2.873484417252736  2.410048757652774 
C  2.864244762235619  1.0327935114418456  0.9278878649785183 
C  3.6440473415719032  2.0017073253957705  1.5555738759663054 
H  3.404471821190274  2.8973974061288654  1.486563809235312 
C  4.7668093704787395  1.6521518252546645  2.2774243374652343 
H  5.255195401414545  2.310779918312616  2.7144274264950954 
C  5.165106244984036  0.36266994101609074  2.360172207406682 
H  5.922069745414657  0.13273683444013984  2.8503295755204348 
C  4.428235413927695  -0.6154842187892449  1.6999877195596034 
H  4.716324334957086  -1.5004432164066404  1.7188245482771598 
C  3.277948399112753  -0.29144918624872096  1.0212855011790898 
C  2.5912748894963  -1.4002507469145815  0.3346621107706319 
H  3.1213844775170223  -2.133019377941835  0.11705581933825543 
C  0.9795959728026882  -2.7059316029809333  -0.7981528527548489 
H  0.008058684068374955  -2.663702216961328  -0.927717954364853 
C  1.2682562497112737  -3.9815054978948243  -0.022131736322258858 
H  1.0320063381279112  -4.741125732876522  -0.5577229310624412 
H  0.7525370366767641  -3.9871619547102473  0.7879865438093401 
H  2.202710322097714  -4.018258915252533  0.1943424490046437 
C  1.5922932071099551  -2.688867887832833  -2.140404102967573 
H  1.2603490543779539  -1.936440536722165  -2.6345537774262535 
H  1.3694325318973686  -3.4997914070453273  -2.603742883377564 
H  2.5462485509872677  -2.620556898881741  -2.056120352004193 
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

