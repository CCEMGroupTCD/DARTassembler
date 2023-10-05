%chk=LIGALIVI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.412242490686356  1.4441852884931343  1.0079459622788327e-16 
N  1.4122424906863558  -1.4441852884931343  -1.3877787807814457e-16 
C  2.6125608890490577  -1.001642852912135  0.061004797041473297 
H  3.2608585944546196  -1.5445323881033315  -0.09563446287052071 
C  2.9259389632375346  0.4156750849652866  0.4587051560918687 
C  4.243290163483007  0.9140504913170453  -0.1400336430147421 
H  4.169072643187428  0.9457281093931407  -1.1163340005862266 
H  4.967390412466268  0.3026232582357964  0.11167813277061853 
H  4.440538263449195  1.8106081528003508  0.2025780273867875 
C  3.052704259526122  0.3817703184577361  2.0002329555379377 
H  3.373355581645947  1.2508066914217046  2.320030519769496 
H  3.686783311766141  -0.31729952254127003  2.2605811318634634 
H  2.175370069400546  0.19135415289290725  2.3965562739641015 
C  1.1936273378853177  -2.845175124161037  -0.25953411205659427 
C  0.7322471232025201  -3.2253958197844796  -1.5260041393840953 
C  0.4860957102744372  -4.577306153282714  -1.7401051484760695 
H  0.18537873248391779  -4.859528402219531  -2.5959056751606693 
C  0.6634726500577898  -5.52634702145973  -0.7467853451081319 
C  1.1052685605053707  -5.1008842467758  0.5018695395057471 
H  1.2266969579365457  -5.744514705696944  1.1910915912759055 
C  1.3752716540518943  -3.7594335614406904  0.7750352003953859 
C  0.5388102261486583  -2.224646188476439  -2.625997517680104 
H  0.133453383989099  -1.4111023763622041  -2.260879747248307 
H  -0.049023852872678964  -2.603696383324481  -3.3136225138275783 
H  1.4059664984852511  -2.0059834985043685  -3.025755053395822 
C  0.3860638219062995  -6.979780307842343  -1.012960858817878 
H  -0.4863430416933605  -7.071702151983625  -1.4481363738849808 
H  0.38082112159217707  -7.472670978223509  -0.16463613847511802 
H  1.0836305878753714  -7.3430563771806545  -1.5980908567342171 
C  1.8144580014695377  -3.311278442317916  2.149744944830943 
H  2.7647002263319904  -3.075597305893145  2.128851021742796 
H  1.6736855303371192  -4.0383459961603965  2.7920133416253154 
H  1.2904347151393118  -2.5288712581407506  2.422642514500666 
C  1.588158786478906  2.8994235683352247  1.0529159015595133 
C  0.8985955631584622  2.90212732804744  2.2612755544474323 
H  0.2581316849663686  2.2258167377220732  2.444701659854817 
C  1.1504088536777728  3.900431179005725  3.2047770424588293 
H  0.696700165623575  3.8917894035358866  4.039252372292431 
C  2.0643348301880526  4.906222074795552  2.9237591899902453 
H  2.226617014561734  5.591063859212262  3.5631524672084414 
C  2.7406373271226663  4.915178463635029  1.7186681868958518 
H  3.3611625589379868  5.608614778891304  1.5305644831461211 
C  2.514266283709304  3.9173156225700527  0.7842715092873981 
H  2.987371884624995  3.9213921424230875  -0.03940158871211688 
C  1.6629094295141302  1.9778552380078325  -1.7144128472348112 
C  2.1141394569540553  1.043034470370289  -2.661740718000876 
H  2.34077637602672  0.16300095934755143  -2.3836162020188416 
C  2.231747137816673  1.3878187677643865  -3.991593447318133 
H  2.5321893295742743  0.7438672246972584  -4.621978400813508 
C  1.9120891755965765  2.6678962063934044  -4.411295756463903 
H  2.0000394981824865  2.9044751708256555  -5.327882388673412 
C  1.464734750748305  3.601417902495635  -3.496767155776587 
H  1.2425004242937618  4.479497149046573  -3.785394910716814 
C  1.3387097612741166  3.259334636631962  -2.1544150927006855 
H  1.0283312343478812  3.90619225070912  -1.5306990217729526 
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


