%chk=IKOZIKUB_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.374021636474478  1.4805960092138564  8.476281519861388e-17 
N  1.374021636474477  -1.4805960092138564  -4.440892098500626e-16 
C  2.20609067053302  1.5562740600141518  -1.600641373071445 
C  3.583420718189207  1.5800114602106168  -1.721336052553754 
H  4.111979426067977  1.5376990572255085  -0.957366718824213 
C  4.18245889719287  1.6661217017664447  -2.9629844312806877 
H  5.10932807386094  1.6780695017671976  -3.034278816878614 
C  3.4077409705695927  1.730484170938977  -4.091121327552829 
H  3.8062889711443915  1.7991852375141926  -4.928056883428685 
C  2.037128984688123  1.6945065648475734  -3.982547085589055 
H  1.514602621841179  1.7438330384740814  -4.751001241104845 
C  1.4270026255179264  1.587735191011481  -2.750826943034378 
H  0.5002160777066923  1.5370938482458663  -2.690112865597168 
C  0.9169346760463541  3.169655857030036  0.3907978775553553 
C  0.30831111803377187  3.458237406177169  1.612107070235072 
H  0.12486478988134397  2.7631001279492766  2.2015923619135953 
C  -0.026108333061692024  4.7387425326273345  1.9640162372993364 
H  -0.4097098931447183  4.913228059666658  2.7924609391085107 
C  0.2140528204206542  5.772384323146944  1.0683042750874507 
H  -0.013058435786079903  6.646250662748239  1.294102941638056 
C  0.7783674987385426  5.512060083718085  -0.1326290095711103 
H  0.9300546079798873  6.211866187849183  -0.7256205925728185 
C  1.1341032829489262  4.226935064361663  -0.4962569876701358 
H  1.5162193399094184  4.064083821686746  -1.3283645719236468 
C  2.646240499664511  1.0585839234648466  1.2195112355469688 
C  3.063562606485437  -0.2908268358285202  1.3264538456789434 
C  4.103065657407074  -0.6033680166234776  2.21186364200603 
H  4.408483160584693  -1.4797151072042047  2.271275893128488 
C  4.677245873989557  0.36683066084203464  2.9967566443505826 
H  5.35659528966692  0.13929903985727132  3.5891853597846377 
C  4.258908710220549  1.670895252387963  2.911102639856678 
H  4.637303031368539  2.3211088616433058  3.456076198905229 
C  3.2599751518176783  2.0112975308377097  1.9989646652482147 
H  3.0055609518618165  2.9022246215438336  1.9155145969034795 
C  2.514429252773003  -1.3983384256177303  0.5485486018843644 
H  3.0699301107870696  -2.1368067236516697  0.44461449033047523 
C  1.088855406433456  -2.7239549494945727  -0.7710938573664181 
H  0.25663779556163924  -3.108863453538558  -0.45372989963106386 
H  1.7952311876296325  -3.368445073108145  -0.6070331176872823 
C  0.9816226421012548  -2.477779468315525  -2.2647235125085463 
H  0.31634263132596496  -1.7699924227753836  -2.3937008183743234 
C  2.2517387095868604  -2.006690234210171  -2.897292345620779 
H  2.646541564334963  -1.2165018285145481  -2.473296827499465 
C  3.28450650171506  -3.1152513270308164  -3.2823565104387575 
C  2.3909151178227805  -3.3274524798400957  -4.546220648461476 
H  2.881055224994434  -3.5573575976991796  -5.362947972679317 
C  1.1908737553436892  -4.178040639718144  -4.274606509259197 
H  1.4683964924992154  -5.103771525436092  -4.183385253021568 
H  0.5837276686417288  -4.123975874919474  -5.028558563842083 
C  0.44627111079674986  -3.7444193270814274  -2.9921429624240834 
H  0.46866441329014463  -4.483099141287266  -2.3658947272263444 
H  -0.4831199235444328  -3.588988099678491  -3.222362094631364 
C  1.997248184330293  -1.8582360315565336  -4.424596376345875 
H  1.0749738662934814  -1.6730798489827619  -4.657642059301901 
H  2.6044338654415284  -1.2466466416318671  -4.870169329727327 
C  4.62229428005284  -2.5021632602313018  -3.578756143498611 
H  5.248550017837964  -3.19325636978473  -3.8049000685088052 
H  4.935811022833114  -2.0285661475575254  -2.80478006708314 
H  4.537921347682979  -1.8911231322220656  -4.314733663702729 
C  3.479079145408422  -4.3400084982748695  -2.3891914765214763 
H  4.1478112113686  -4.911633678089707  -2.7736765058247017 
H  2.6513824362642016  -4.819600142843078  -2.318629431963568 
H  3.762336635247728  -4.058457756272439  -1.5163020080732568 
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

