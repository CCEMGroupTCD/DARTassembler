%chk=QICESACU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.2320689865019734  1.600657993607629  -5.897770124808791e-16 
N  1.2320689865019734  -1.600657993607629  -2.220446049250313e-16 
C  3.278422288603483  -0.2576482353800642  -0.48594428928005007 
C  2.963605828975431  1.0448945491037778  -0.06603027268099354 
C  3.9758655042946938  1.8986451725708084  0.3440766922101392 
C  5.308584259126639  1.4948734595746167  0.3549253927485314 
C  5.6186500957037415  0.21089259157677853  -0.07200409511012773 
C  4.6431738053185505  -0.6341310020963882  -0.4645033473795547 
C  2.279501073681744  -1.2163827502661846  -0.9976036939091304 
C  1.7432729353655523  -0.8519159857669552  -2.3478680240509995 
C  1.8620673633252132  -1.976201306045771  1.293569776116902 
C  0.48957546258266293  -2.774677082947617  -0.5620630879209272 
C  0.9248175779247049  2.549032437781439  -1.5197093325913789 
C  -0.2803021212739156  2.420055460298136  -2.209091454409724 
C  -0.530798524050994  3.215106679910485  -3.3102208566094466 
C  0.38046305608194  4.141437630008164  -3.6920872252055776 
C  1.5717219852307607  4.2724513323964874  -3.052181065280741 
C  1.8434685763515821  3.482320389340062  -1.9500170742889176 
C  1.1920105840556217  2.864977819005685  1.3041078841487421 
C  1.5169314332478656  2.340339896130418  2.7037213742385493 
C  1.3155935275106971  3.380421754121687  3.787696035561792 
C  -0.13458845999714186  3.670038731546808  4.00825201370562 
H  3.7931031818952285  2.802193207676255  0.6708128422888311 
H  5.950304422752826  2.0792931942936064  0.671107959103112 
H  6.462030520508812  -0.09234416914187715  -0.0801272721037446 
H  4.769514353284993  -1.5013875869157727  -0.7207870027830197 
H  -0.9248516373331803  1.8177672551662156  -1.9214958888897655 
H  -1.3991431393853726  3.117885514851195  -3.7173472056515005 
H  0.25754504219391317  4.612292892662573  -4.3196768231853895 
H  2.265451384348247  4.959589730725452  -3.183123275213263 
H  2.568850576368808  3.5929384168681864  -1.4298059574941981 
H  1.8439319647508103  3.6331153755044623  1.038580808902977 
H  0.48108636403794003  3.1616564111199787  1.2806520453585473 
H  0.9061206217456734  1.6735009965886785  2.7976201208222697 
H  2.47106339802752  2.0189965691647913  2.7377291674688147 
H  1.702528048306612  3.0928995464808624  4.548933098985406 
H  1.9340050758596394  4.175142747038445  3.595253445977052 
H  2.6532665365654564  -1.8847350649220052  -1.1177850246120973 
H  -0.45746406447395027  4.323644295120095  3.2783555270037157 
H  -0.18357879481283979  4.021966121563695  4.703353063374356 
H  -0.6872695314102628  2.914185844901634  3.920731347278714 
H  1.2472554913361178  -1.5306302703540395  -2.7378141679169925 
H  2.281384864885177  -0.5929524085034856  -2.9027283460777067 
H  1.2955990821663135  -0.1292311151576997  -2.264903399335183 
H  1.0730170423873584  -2.1961669710448546  1.9468339931450576 
H  2.5608111553408  -1.172931831782814  1.888854451576459 
H  2.5700636154289622  -2.7085460009658724  1.0771018689587262 
H  -0.17262332648043488  -3.0487779381909936  0.08058934140478519 
H  1.09725230261991  -3.510360411118131  -0.617333231165788 
H  -0.024906137913503823  -2.4732865317326787  -1.3422609400003076 
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


