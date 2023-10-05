%chk=IWARAQUV_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.4666599631475596  -1.3888875233437727  3.0147936227148825e-15 
N  1.4666599631475594  1.3888875233437732  -3.92134270908896e-16 
N  2.8516102782861816  -0.4353732050556605  -0.31845958835782245 
C  2.6484985298422523  0.9322763152070668  -0.2963040615263237 
C  1.3147519045384457  2.8470637694116983  0.25098636842041006 
H  2.2330643744512013  3.239712379060861  0.30077659781664967 
C  0.6487516470411057  3.069386970398198  1.6080696051147048 
H  1.218325904157334  2.6869809072194326  2.3232578380119064 
H  -0.22363475731547888  2.60198482365192  1.628911982876918 
C  0.4350131306280438  4.561692700821975  1.867007249241433 
H  -0.04427921279167624  4.68536013761087  2.7249084700223034 
H  1.3123778393760752  5.015385299689127  1.9359153298754694 
C  -0.37680545615537087  5.1842874973525  0.7343834319542766 
H  -1.2688720455444242  4.758938881952936  0.698473886239317 
H  -0.504397255672306  6.150413718840538  0.9130793939939169 
C  0.3220905551673787  5.007594193628972  -0.6090938554519817 
H  1.1802796579372912  5.500499197449465  -0.5999215549932104 
H  -0.24326703879164735  5.3902675066704955  -1.3257237295422963 
C  0.5890677490506504  3.525809495149393  -0.9086093396894486 
H  -0.2711492490961671  3.063073723391061  -1.073103387103122 
H  1.1387995000894455  3.45107588407187  -1.7291109541151353 
C  4.100871496780515  -1.161305667587003  -0.7122965569237862 
H  3.914780638596677  -2.1176597010670988  -0.486184781142986 
C  4.358827404715923  -1.168586833145789  -2.2191296910794573 
H  4.666264163194107  -0.27408284767619395  -2.512737200376189 
H  3.52271563965779  -1.3841429994789924  -2.7030913274589157 
C  5.427386742869576  -2.2162430583935153  -2.5332378007237186 
H  5.073978532213793  -3.1171016081108696  -2.324973569415954 
H  5.637278805701566  -2.1889100783199344  -3.5005527908693597 
C  6.7044468879307555  -1.9829484974200555  -1.7364373152089563 
H  7.136127221271487  -1.1484503213477224  -2.0497190956220637 
H  7.332144096846614  -2.7304428031508396  -1.9011038473863053 
C  6.4328411857647225  -1.8755175661947723  -0.23698574961337993 
H  7.2737407271540775  -1.6412311754461526  0.2307957788816812 
H  6.130100237248788  -2.753938931852054  0.10334540154452243 
C  5.372752309038299  -0.8217510333394387  0.07467380450766246 
H  5.700798190070305  0.07647597139060311  -0.1836567548763156 
H  5.178456588407737  -0.8152484081461039  1.0464053460407519 
C  3.8096797787109535  1.8229139796753242  -0.6359773573655584 
H  3.4874300940242766  2.606209611405702  -1.126628377981482 
H  4.448904564724526  1.3301163099603583  -1.1918955075184394 
H  4.25066047378023  2.1134626441405864  0.18959987109565066 
C  1.3231983426045963  -2.5541080453798313  -1.3705312832829608 
C  0.6940481147551602  -2.098361698493448  -2.5357008452338348 
H  0.3549465892097401  -1.2122498828749715  -2.5766844467605328 
C  0.5659719477406088  -2.9417054288626647  -3.631242203599554 
H  0.15420528945249057  -2.6280490940781145  -4.427596171864743 
C  1.0422354529918159  -4.245540307183318  -3.5583994713580256 
H  0.9378090471774847  -4.828561823605689  -4.300902809503558 
C  1.6673908677644884  -4.6993318125477135  -2.411474636401715 
H  1.9918083461140093  -5.591622400641113  -2.3697619747595624 
C  1.8220928218088854  -3.8540197118432094  -1.3174606886764617 
H  2.2669215810806946  -4.1619365105846935  -0.5367724538621247 
C  1.8394820335283362  -2.350178050877394  1.4765742939355209 
C  2.794495750470433  -1.9137825726166373  2.3930074486205926 
H  3.2533952641945953  -1.0936198734388378  2.2503180195057277 
C  3.07190757157454  -2.683801312093274  3.518713907856992 
H  3.729440316067806  -2.391714121315344  4.137678157957545 
C  2.3986813656394474  -3.8719041694270726  3.7434380187964336 
H  2.615281823631896  -4.405866432414756  4.499562243578646 
C  1.4065191294951878  -4.282078870707758  2.8615035361966474 
H  0.9244908446208381  -5.083620400377366  3.028693005948807 
C  1.1167484399856915  -3.523230153813402  1.735062933071473 
H  0.42905334635938286  -3.7991521910978783  1.140082229566642 
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


