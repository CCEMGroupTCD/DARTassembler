%chk=ECEHUCIZ_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.3688394418265408  -1.4853883608336238  1.806183623491938e-15 
N  1.3688394418265428  1.4853883608336238  -1.5374156464789625e-17 
C  1.4164693073759458  2.446163079843188  -1.210864940949187 
C  0.7536478487828134  3.7580244161649574  -0.7908796456202329 
C  0.4793287010403928  3.566517711588144  0.6720210096821189 
C  -0.07161123310202755  4.452866251616506  1.632945936279801 
C  -0.13301884325146318  4.0220012523507105  2.981263242340762 
C  0.3476972204349775  2.772695273570487  3.3965827909991577 
C  0.8130025431154773  1.8553658694035569  2.4065896794393016 
C  0.8748914715442206  2.3121192565239737  1.12215554570881 
C  1.1027493412326148  -2.5698957007833787  1.4227750025780885 
C  1.692712825409931  -2.291999400472838  2.6681275685653634 
C  1.4350113285493347  -3.11740666925736  3.7458580168470945 
C  0.6024530229957854  -4.219553557685398  3.6135673362145257 
C  0.020482523567742383  -4.487256035693799  2.4154619731250997 
C  0.2590164743553087  -3.6671613102198855  1.3216728181360793 
C  1.7771804365723691  -2.554327078724821  -1.3988574549978827 
C  2.5142518551994373  -3.716938517951946  -1.2120744653512567 
C  2.8868799484005265  -4.496199905425312  -2.3465784204569373 
C  2.5042924799068773  -4.090412081921886  -3.6052957522070117 
C  1.8090335500840828  -2.946341982270455  -3.7862912979677406 
C  1.4223258162317547  -2.158243252761747  -2.6995008116615877 
C  2.8219348956515566  -0.49095486859121923  0.3451146499059249 
C  4.080457735807375  -1.1031232778614277  0.6178923240249293 
C  5.169631501519814  -0.31742605498496423  0.8291173747362566 
C  5.033269260206905  1.0648714566554824  0.8169379471328587 
C  3.838750700906847  1.716979902191159  0.5758017890213205 
C  2.711473425488885  0.8952794720659267  0.3045251636768214 
C  3.8547436272030247  3.2236257760985216  0.5873869927918558 
H  0.945774429414295  2.0666554205239547  -1.9437404818009205 
H  2.318296671328257  2.6053929528087068  -1.4625946237087475 
H  -0.05251923795744662  3.8990594818244366  -1.27184304184822 
H  1.3404296734380297  4.49256610827009  -0.9319266335839366 
H  -0.38498453942162114  5.315485179289331  1.381046876124931 
H  -0.5221350341023261  4.596102051498899  3.626652522105392 
H  0.36577284889204864  2.5473904483649226  4.321713093926815 
H  1.0637595591639801  0.967419151375024  2.6347393558467225 
H  2.2649782517782744  -1.539583621899773  2.769283576482641 
H  1.8308546356909106  -2.9301515665959257  4.587868184670789 
H  0.4364769396605188  -4.783303885263057  4.362355048397888 
H  -0.5518991288384729  -5.241429920139281  2.326677931714786 
H  -0.15841822865925126  -3.86212056305628  0.489958168479936 
H  2.7695399849297178  -3.988709558202621  -0.33628734938077764 
H  3.4034160912605267  -5.2864284180152765  -2.2353681645936456 
H  2.729346075287972  -4.624753100693948  -4.357131847951087 
H  1.579997192083289  -2.6734400189263368  -4.668582681308118 
H  0.9213145078725091  -1.3635962862690496  -2.83950580278733 
H  4.155996423764547  -2.0511631354039292  0.6580189377931263 
H  6.02150153298  -0.7107610207121139  0.9827804940474968 
H  5.805885507880614  1.594130760638986  0.9800027286692297 
H  4.454157354788709  3.5365126528606154  -0.07910723488901654 
H  4.1423546215178675  3.5248699333833375  1.4408654598025958 
H  2.9822665200683014  3.547745159752374  0.41157998876272867 
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

