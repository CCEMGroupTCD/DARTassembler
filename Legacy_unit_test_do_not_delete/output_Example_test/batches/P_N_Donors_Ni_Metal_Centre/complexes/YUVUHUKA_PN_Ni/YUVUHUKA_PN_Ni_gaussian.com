%chk=YUVUHUKA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.4330550085394476  1.423535507987069  -2.2307064969184675e-16 
N  1.4330550085394476  -1.423535507987069  2.498001805406602e-16 
N  3.5519390338495045  -2.0484792343700993  0.18802863117197557 
C  2.6558362416343235  -1.3155018171257375  -0.518634387242004 
C  2.8646009802779964  -2.642881920676458  1.2511609355557427 
C  3.255092175776417  -3.4138883787815897  2.32486951696803 
H  4.158071014487294  -3.697272321603245  2.420226623631515 
C  2.313764427182613  -3.7650228444801743  3.2412325545933505 
H  2.5727379311789944  -4.273593256625155  4.000555483338803 
C  0.9631717066269103  -3.3892414171827094  3.091958403849755 
H  0.3224021902969343  -3.680911538981857  3.73005467949471 
C  0.5596091118452866  -2.604525395887548  2.0342721942318946 
H  -0.342188516722872  -2.3219078512300997  1.94109622865799 
C  1.5293589271192414  -2.2524464753492652  1.1127106827261337 
C  5.003372426355679  -2.218063771352503  -0.09689888932141494 
H  5.184194237551691  -1.7626274184393536  -0.910355194872589 
C  5.3255645770234405  -3.694749033648819  -0.3265289269143745 
H  4.758381120324934  -4.038102955807844  -1.0045102692723877 
H  5.188128174345888  -4.177882282620256  0.47963572063433757 
H  6.23227201936822  -3.779895950827861  -0.5968271076177321 
C  5.867838904381719  -1.5835675224481132  0.9738801628156405 
H  5.626793571416138  -0.6699309340365716  1.0748890325807026 
H  6.7806891175669515  -1.644267978822528  0.7197615037532419 
H  5.736982726377493  -2.0392364722990592  1.7971945309763528 
C  2.983762207238681  -0.33481048405596314  -1.5572464468312808 
C  2.647497168486309  1.0071056109638143  -1.3130340441164743 
C  3.140779519381206  1.976147529665663  -2.2108705838575506 
H  2.9213054951657322  2.888401644658765  -2.062844936171692 
C  3.9217072310642935  1.6547074549835532  -3.280224064944998 
H  4.310773679140279  2.3450125014586742  -3.8013415294312067 
C  4.159357230937072  0.3082281966541487  -3.6237661794654366 
C  4.866241464111821  -0.08613850993670757  -4.770308582354083 
H  5.232327916821687  0.5848913587532105  -5.336878520197863 
C  5.044130070624987  -1.3798895155608912  -5.092929721027947 
H  5.533407100926085  -1.616931326676053  -5.872112419374737 
C  4.496172169245105  -2.3857947191654145  -4.266536812251071 
H  4.596516117899749  -3.298495408898965  -4.506892388280193 
C  3.81376232671333  -2.050070532927314  -3.106608259455432 
H  3.447979961406782  -2.732728731158098  -2.556777780472073 
C  3.664375161351148  -0.6947496510982609  -2.7462587733961574 
C  2.290992036556016  1.1571817409945258  1.5624266852454425 
C  1.5322024175066142  0.7703896294450314  2.6653875351698324 
H  0.5892172778633132  0.6783305769909709  2.597207499432137 
C  2.205475284594727  0.5132154648766485  3.8992625310756415 
H  1.700676713826818  0.2353555254647317  4.654079288322521 
C  3.5336377347846684  0.6503212018989641  4.011515071274345 
H  3.9648991748103084  0.48962966753501735  4.843637268509102 
C  4.256074012188653  1.0266467146763119  2.9261685581627184 
H  5.1988426264117695  1.1064137976629298  3.001017222196698 
C  3.643520803008908  1.2924149535978504  1.7015106675237812 
H  4.169230072780259  1.5729021193620567  0.9607569948916171 
C  1.1462337325736152  3.212513675873509  -0.1003668295753396 
C  0.5442488192465214  3.733009069801351  -1.2160867671163433 
H  0.3105671964921861  3.158467974471376  -1.9364987318725866 
C  0.2697066607937204  5.1051064983803816  -1.3043456460968261 
H  -0.12054400484345273  5.471577715341931  -2.0892096311977237 
C  0.5802993335652504  5.923360597743937  -0.21799995503186725 
H  0.38549290000852054  6.853217592931742  -0.2519509777736896 
C  1.1604563534658159  5.392800038655402  0.8860973529078734 
H  1.362274431934173  5.9644100901325245  1.619438877883568 
C  1.4663424181185887  4.0601893595705825  0.9742777557213032 
H  1.8876042466449943  3.7148863184858762  1.7531964816768777 
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
