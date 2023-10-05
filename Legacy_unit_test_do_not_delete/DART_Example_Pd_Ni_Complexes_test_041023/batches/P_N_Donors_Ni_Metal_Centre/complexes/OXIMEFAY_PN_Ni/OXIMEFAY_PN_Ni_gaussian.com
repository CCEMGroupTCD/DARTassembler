%chk=OXIMEFAY_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.4896503868357842  1.3642000311537896  9.264444108133585e-16 
N  1.489650386835783  -1.3642000311537896  -1.1102230246251565e-16 
N  3.5892025983189355  -1.921828599272656  -0.23522561780981366 
C  1.5468868637328308  -2.7023872055493317  -0.3435150098550179 
C  2.8203256470867246  -3.045592981321096  -0.47507602432012674 
C  5.0637774264938065  -1.8372156009117226  -0.2186688481523721 
C  2.74692301082223  -0.9301487875841815  0.059820205031143714 
C  3.04692513912703  0.47816992384048884  0.37858856930800755 
C  1.5094194224061885  1.6884025746956919  -1.7863473447215243 
C  0.3125468916495797  1.6363955937692989  -2.4838763470255674 
C  0.31857494785430696  1.868572593305149  -3.8629721788751725 
C  1.4798861280080073  2.1529800380637116  -4.498708545788352 
C  2.662654349220573  2.203224938690219  -3.8157165144653282 
C  2.691731596881801  1.9889558884595284  -2.4399632852189908 
C  1.6262450935827777  2.962298250046316  0.8342183939840037 
C  1.7091402574521348  2.995998146652473  2.2204603457283114 
C  1.7253001878701655  4.21033377351719  2.876524010303596 
C  1.6907817655540278  5.394769657625336  2.161134770128884 
C  1.6315270756226932  5.368971864635817  0.8042760576617191 
C  1.585345818533561  4.159451338086299  0.11833171782863004 
H  0.7814249842206753  -3.2793068155968093  -0.4517597074452321 
H  3.1615754419887776  -3.9085397032943696  -0.7033348084115497 
H  5.459735931957839  -2.6709521417605337  -0.47399924666408777 
H  5.3690108643134815  -1.1668839493088563  -0.8421567191620373 
H  5.384405193524985  -1.601010189211601  0.6473598341048606 
H  3.763483880890615  0.8056245087976426  -0.16588163593075067 
H  3.2854252557911066  0.5816985829628241  1.2933413561808145 
H  -0.4960202973550074  1.4348018181355275  -2.028865027945541 
H  -0.4925405609987681  1.835402395411422  -4.3659385289123325 
H  1.4729026046682034  2.3397157287120596  -5.441328953078504 
H  3.4751077053543207  2.3869498172096986  -4.2886989275844805 
H  3.5132970593282575  2.05638668152453  -1.9519147127779108 
H  1.735943799056542  2.17607847984673  2.7237109097047676 
H  1.7534233978510725  4.235248466541535  3.830285276623371 
H  1.7299898000314156  6.2240738096294645  2.627746228947825 
H  1.6254076069479855  6.188825763960288  0.31678952875491756 
H  1.5265715164777216  4.1524900239218585  -0.8475780768846105 
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

