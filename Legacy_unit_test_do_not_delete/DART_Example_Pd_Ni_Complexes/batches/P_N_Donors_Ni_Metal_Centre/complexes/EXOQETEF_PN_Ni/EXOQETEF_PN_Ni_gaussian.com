%chk=EXOQETEF_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.4490272564724238  1.4072739640880165  -1.4956705438325699e-15 
N  1.4490272564724243  -1.4072739640880165  3.3306690738754696e-16 
C  2.508693189932022  -1.0132205145456807  0.6270114100788206 
C  2.5598815451205463  0.40856210281173166  1.0764831013713865 
C  3.689614510509938  -1.878402576511156  0.9084780941778549 
C  4.0886994573589455  -2.0867878074494985  2.22691084207413 
C  5.113891366519061  -2.971349216078064  2.497111145104391 
C  5.78095338909141  -3.606704337267586  1.502638391319775 
C  5.43883025293381  -3.346483819198586  0.1688295213231884 
C  4.390419919073989  -2.5108778959680937  -0.1056829973678845 
C  1.3393103166450986  -2.7588993821228485  -0.47970399460198204 
C  1.3784042731215127  -3.851035774339845  0.37822552660723113 
C  1.1978248770085465  -5.116162796819633  -0.1342009351323111 
C  1.0017070282680391  -5.330879559661059  -1.473966754654637 
C  0.9632794540121807  -4.230969300891297  -2.3374400050740407 
C  1.1269895688645626  -2.956186972943818  -1.8359434301626814 
C  2.379872964458288  1.7504864034036287  -1.5222250127439834 
C  1.8132983468821195  1.449611384889097  -2.745614041611057 
C  2.4950177636965214  1.7103633915090737  -3.915175676490508 
C  3.76280963926602  2.258245758180769  -3.859182975803653 
C  4.34775040250037  2.5269871402080333  -2.6568237852566394 
C  3.6681568922001624  2.291082662137809  -1.4881626205111487 
C  1.240814392835737  2.9991794474327538  0.8396678959643392 
C  0.6552295401274371  3.008798896342814  2.1034087303450644 
C  0.3787369481074947  4.209690480445411  2.7472174440090407 
C  0.6946318774399789  5.401088774441802  2.130368792649146 
C  1.2386818046771757  5.407562183481173  0.8738335039888321 
C  1.5127061913731292  4.209693166489562  0.22573386378875368 
H  3.4678166102759973  0.7426851250511959  1.0146780573153147 
H  2.2719760417193786  0.4734518913890242  2.0003123964523195 
H  3.6688176930754657  -1.633104007213948  2.9210835608020056 
H  5.352313938588058  -3.1375025422724616  3.381865178337658 
H  6.4597761433751355  -4.21022898966205  1.7029756125713296 
H  5.9218056810388155  -3.7370505560430276  -0.5239058754319146 
H  4.143515003624473  -2.36410351139838  -0.9899611395613231 
H  1.5256541672858743  -3.7277368329624663  1.2877764870947899 
H  1.208571971408887  -5.844762050529604  0.44382663121313093 
H  0.8960335582268704  -6.193355908242244  -1.804247129287352 
H  0.8255107508249881  -4.357698199848493  -3.2486010614143885 
H  1.0954963712111554  -2.227174871405475  -2.410346254243671 
H  0.9667331206855172  1.0669653161479848  -2.7812372625727093 
H  2.1048475497798154  1.5187970036702234  -4.736942731818586 
H  4.218985819940089  2.446499564305636  -4.64700506654823 
H  5.212036106753351  2.8702258785961394  -2.62859968396574 
H  4.066096301203883  2.4931073192089315  -0.672171741458523 
H  0.44913086555042003  2.2037678512796197  2.5204593301792633 
H  -0.015075761923453168  4.2083197096863225  3.588967436824979 
H  0.537013975448504  6.205422521022226  2.5702332502763796 
H  1.425802891641778  6.21618390820515  0.45314584655447093 
H  1.879963062450101  4.221233748078112  -0.6281104402991569 
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
-C 0
6-31g(d)
****
-N 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****
-P 0
6-31+g(d)
****

-Ni 0
lanl2dz

