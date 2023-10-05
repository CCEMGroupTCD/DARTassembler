%chk=IXEPIDOZ_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.3834487287210895  1.471791294647444  9.152119315094444e-17 
N  0.9723590338330448  -1.0286024251574617  -2.2802657123057224 
N  1.38344872872109  -1.4717912946474445  4.440892098500626e-16 
C  1.3603361853461395  1.4700145243882974  -1.8293304959537156 
C  1.4845156346105544  2.68068198211832  -2.525309447628199 
H  1.7563951167520138  3.4357346985200694  -2.0557260123956143 
C  1.2228190605679816  2.8055196175005914  -3.88158784011499 
H  1.3280662003469168  3.629841003564109  -4.301096376436563 
C  0.8048995305309115  1.7112300630934911  -4.610942813724657 
H  0.5642212378328202  1.8019466344481156  -5.504674262821977 
C  0.7513653895618762  0.4762956185088172  -3.9871045147007123 
C  0.4072464144242678  -0.8333688509102182  -4.491311772176854 
C  -0.05203653291390853  -1.2951677469030163  -5.730258390859138 
H  -0.18919466211865044  -0.7045336718399267  -6.4355434941808145 
C  -0.29640742892135674  -2.6405889130706948  -5.882985863216044 
H  -0.5879287152252226  -2.966136102908349  -6.704441167429545 
C  -0.11243663529514314  -3.5176117115437604  -4.82172218813286 
H  -0.28206586779474563  -4.422630096177615  -4.9501025096635285 
C  0.3160872149621028  -3.089342491381678  -3.5840439485154048 
H  0.42564258043715064  -3.683155274606052  -2.87757848554148 
C  0.5749266087662674  -1.7358017874483818  -3.438060041830998 
C  1.0657743042512835  0.34086136521154753  -2.6204371457342686 
C  1.7045751012067711  -1.700347946139652  -1.2754837082271615 
C  2.70406004197154  -2.6021928917016854  -1.6128005731311572 
H  2.9274020962061753  -2.74537361016794  -2.5039118224540156 
C  3.356094420238649  -3.2776870827705  -0.6200590924877187 
H  4.032202358922274  -3.879300092183022  -0.8329217527984961 
C  3.010446994135203  -3.065963449716505  0.6989061913747892 
H  3.4386393060537674  -3.5284676427196704  1.3820669231709815 
C  2.0254957922632704  -2.161872360280661  0.9730738261054009 
H  1.7864902415715997  -2.015582495146232  1.8597434444512388 
C  1.041354850082259  3.227498900336867  0.4695141760609155 
H  1.7326715735397142  3.777567230639044  0.04404120410307827 
C  -0.31781634540294257  3.7279572755689605  -0.05739324087629838 
H  -0.3620661803232186  3.600442392390529  -1.017431425886419 
H  -1.0335160573752404  3.215354484844667  0.34816558497323236 
C  -0.4982821539351161  5.2063493294649  0.27102514613663553 
H  0.17317417286405012  5.725307144868315  -0.198115096256953 
H  -1.3699009626906913  5.499383098293729  -0.03618563256336782 
C  -0.3792587399203835  5.460999215194644  1.7613360373986335 
H  -1.108877254093673  5.019523666070917  2.2231449516710673 
H  -0.4490589708369257  6.413314027150196  1.9323615229090099 
C  0.9471310573809515  4.945303880020227  2.3032125887406205 
H  0.9707604725869383  5.071786368180747  3.2646552136152174 
H  1.6718631117927334  5.460042398147878  1.9156791743569883 
C  1.1473006704995528  3.4652609806485497  1.9826239395288563 
H  0.47384107844869594  2.939304421579874  2.4422987439801496 
H  2.019068208267054  3.1810497466772283  2.296866130488782 
C  3.152455519428706  1.1169430957898048  0.4171335519925433 
H  3.3722425048073434  0.2915243530522307  -0.06604720113789442 
C  3.406018374756693  0.8012951967367787  1.8918561334367126 
H  3.283314028059694  1.6054025245243724  2.4224600334411672 
H  2.7676857547669247  0.13899348755762575  2.1975270876920683 
C  4.821267752098788  0.2756199479337005  2.0819473292578206 
H  4.983968259458949  0.11760791467479481  3.0250577799740404 
H  4.914899988397733  -0.5715260533603665  1.6169268492602944 
C  5.848685867253744  1.25260597185528  1.5531872459370863 
H  6.731407735022337  0.8528093993793224  1.6145977890026106 
H  5.843991410061992  2.050412562552143  2.1048448081188424 
C  5.58037778118192  1.6431059952554536  0.11287276893112463 
H  5.7267050944701445  0.8734114916122862  -0.4576024614228025 
H  6.205947263496297  2.3349065166356935  -0.15284405668397882 
C  4.146769307417996  2.1544443030284146  -0.08622818273982413 
H  4.02662450818668  2.985641402405907  0.3995936562136221 
H  3.990498174496831  2.3288431739372717  -1.0260959962334306 
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

