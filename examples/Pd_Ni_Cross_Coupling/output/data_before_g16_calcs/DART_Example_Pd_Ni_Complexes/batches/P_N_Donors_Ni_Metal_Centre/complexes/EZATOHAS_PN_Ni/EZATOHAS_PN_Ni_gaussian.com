%chk=EZATOHAS_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.3792073275979981  1.4757666304331454  4.093512522797002e-15 
N  1.3792073275979948  -1.4757666304331458  1.1102230246251565e-15 
H  1.4090622586531918  -1.8630570313769912  0.8126737896772975 
H  1.2325009186285483  -2.126744434232676  -0.6034598151721806 
C  2.7086461622868265  -0.8711260411970932  -0.27759674470265083 
H  2.8068126490773544  -0.706831328198443  -1.2288491392604397 
H  3.4129543512210896  -1.476334823250835  0.002039867930898165 
C  2.809578181122938  0.44246554233873003  0.49283568532449207 
H  3.6407179957759146  0.8945939386856261  0.2787393945686395 
H  2.7888629612832263  0.27328355255033165  1.4486175139853235 
C  2.040086157048208  2.29824541533679  -1.51259097946801 
C  3.1555152281808754  3.1409715275507386  -1.4454212164426228 
H  3.539946915726974  3.329584252123798  -0.6186774666993307 
C  3.695125763292861  3.6985285985887897  -2.5902207110842186 
H  4.442817554060842  4.247369379574787  -2.5303279473741807 
C  3.133266832963173  3.446152204977908  -3.8015105023219853 
H  3.4948955980031755  3.8284723620845904  -4.569282261222754 
C  2.0315270556855216  2.629145260833819  -3.8938128662383162 
H  1.651039603581626  2.448324269562537  -4.724149929067134 
C  1.488870424002143  2.075292982875683  -2.7436527929047507 
H  0.7334903427564191  1.5396975820174794  -2.809985284466231 
C  1.3417860187416715  2.8653838488398895  1.18406000352643 
C  1.6289734575971813  2.7049311457383833  2.522927530592133 
H  1.9467330298355818  1.8846977531153781  2.8268222095140167 
C  1.4549268658997363  3.7341187209813684  3.4256304933013357 
H  1.687542231024472  3.61201986653911  4.317647966029207 
C  0.9453019049885789  4.936352785915645  3.018186555651104 
H  0.809014350110766  5.621424855679999  3.631923406112464 
C  0.6318229325587159  5.123698911175878  1.6929516771505908 
H  0.2712765846633618  5.933662143403976  1.413167486285667 
C  0.8532033756993693  4.100508320698297  0.7662020373054961 
H  0.6744776164058995  4.2463593892624525  -0.1344775776318509 
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
-H 0
6-31g(d,p)
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
-Br 0
6-31+g(d)
****

-Ni 0
lanl2dz

