%chk=ETUKEVUB_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.4666599631475596  1.3888875233437727  -2.844703956731018e-15 
N  1.4666599631475594  -1.3888875233437732  2.220446049250313e-16 
N  2.8516102782861816  0.43537320505566046  0.3184595883578225 
C  2.6484985298422523  -0.9322763152070668  0.2963040615263236 
C  1.3147519045384457  -2.8470637694116983  -0.2509863684204104 
H  2.2330643744512013  -3.239712379060861  -0.30077659781665006 
C  0.6487516470411057  -3.069386970398198  -1.6080696051147052 
H  1.218325904157334  -2.686980907219432  -2.323257838011907 
H  -0.22363475731547888  -2.60198482365192  -1.6289119828769183 
C  0.4350131306280438  -4.561692700821975  -1.8670072492414336 
H  -0.04427921279167624  -4.68536013761087  -2.724908470022304 
H  1.3123778393760752  -5.015385299689127  -1.93591532987547 
C  -0.37680545615537087  -5.1842874973525  -0.7343834319542772 
H  -1.2688720455444242  -4.758938881952936  -0.6984738862393176 
H  -0.504397255672306  -6.150413718840538  -0.9130793939939177 
C  0.3220905551673787  -5.007594193628972  0.609093855451981 
H  1.1802796579372912  -5.500499197449465  0.5999215549932098 
H  -0.24326703879164735  -5.3902675066704955  1.3257237295422957 
C  0.5890677490506504  -3.525809495149393  0.9086093396894481 
H  -0.2711492490961671  -3.063073723391061  1.0731033871031215 
H  1.1387995000894455  -3.45107588407187  1.7291109541151348 
C  4.100871496780515  1.161305667587003  0.7122965569237864 
H  3.914780638596677  2.1176597010670988  0.4861847811429863 
C  4.358827404715923  1.1685868331457887  2.2191296910794573 
H  4.666264163194107  0.2740828476761936  2.512737200376189 
H  3.52271563965779  1.3841429994789922  2.7030913274589157 
C  5.427386742869576  2.216243058393515  2.533237800723719 
H  5.073978532213793  3.117101608110869  2.3249735694159543 
H  5.637278805701566  2.188910078319934  3.50055279086936 
C  6.7044468879307555  1.9829484974200553  1.7364373152089565 
H  7.136127221271487  1.1484503213477222  2.0497190956220637 
H  7.332144096846614  2.730442803150839  1.9011038473863058 
C  6.4328411857647225  1.8755175661947723  0.23698574961338015 
H  7.2737407271540775  1.6412311754461526  -0.230795778881681 
H  6.130100237248788  2.753938931852054  -0.1033454015445221 
C  5.372752309038299  0.8217510333394387  -0.07467380450766237 
H  5.700798190070305  -0.07647597139060314  0.1836567548763156 
H  5.178456588407737  0.815248408146104  -1.0464053460407519 
C  3.8096797787109535  -1.8229139796753242  0.6359773573655582 
H  3.4874300940242766  -2.606209611405702  1.1266283779814819 
H  4.448904564724526  -1.3301163099603586  1.1918955075184392 
H  4.25066047378023  -2.1134626441405864  -0.1895998710956509 
C  1.3231983426045963  2.5541080453798313  1.370531283282961 
C  0.6940481147551602  2.0983616984934477  2.535700845233835 
H  0.3549465892097401  1.2122498828749713  2.5766844467605328 
C  0.5659719477406088  2.9417054288626643  3.6312422035995544 
H  0.15420528945249057  2.628049094078114  4.427596171864743 
C  1.0422354529918159  4.245540307183318  3.558399471358026 
H  0.9378090471774847  4.828561823605688  4.300902809503559 
C  1.6673908677644884  4.6993318125477135  2.4114746364017154 
H  1.9918083461140093  5.591622400641113  2.3697619747595633 
C  1.8220928218088854  3.8540197118432094  1.3174606886764622 
H  2.2669215810806946  4.1619365105846935  0.5367724538621252 
C  1.8394820335283362  2.350178050877394  -1.4765742939355206 
C  2.794495750470433  1.9137825726166375  -2.393007448620592 
H  3.2533952641945953  1.093619873438838  -2.2503180195057277 
C  3.07190757157454  2.6838013120932747  -3.5187139078569913 
H  3.729440316067806  2.3917141213153443  -4.137678157957545 
C  2.3986813656394474  3.871904169427073  -3.743438018796433 
H  2.615281823631896  4.405866432414757  -4.499562243578645 
C  1.4065191294951878  4.282078870707758  -2.861503536196647 
H  0.9244908446208381  5.083620400377366  -3.0286930059488064 
C  1.1167484399856915  3.523230153813402  -1.7350629330714726 
H  0.42905334635938286  3.7991521910978783  -1.1400822295666415 
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


