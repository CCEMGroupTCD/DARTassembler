%chk=ZIKADOHE_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3692743497926192  1.4849874595430097  -3.024340535318518e-16 
N  1.3692743497926188  -1.4849874595430097  4.996003610813204e-16 
C  2.5027630766117612  -1.439462172736342  -0.5851587002131331 
C  1.0661151992463016  -2.706236493227173  0.7085575891665752 
C  0.7225858431666123  -2.606697329411278  2.0522484823329905 
C  0.3900058906538063  -3.7538292927507055  2.7371042509964485 
C  0.3569770557119909  -4.9677555070197155  2.0971091163486073 
C  0.6905290913705113  -5.04965102015585  0.7670018408116817 
C  1.0625163970034062  -3.9097433838173408  0.07261718871665118 
C  2.175532387198646  1.3536174882032959  1.629708761634566 
C  1.3761252210683406  1.0825169159978647  2.726802537247855 
C  1.940331483119353  0.8650179041884352  3.9621717669619634 
C  3.289697358019708  0.9409620909728176  4.108492360687123 
C  4.091813442229106  1.237793557236949  3.035589551441061 
C  3.5326463066857796  1.4541123839958332  1.7934813034538029 
C  1.064609645925644  3.2461006381046738  -0.2958507233624465 
C  0.3429113016944474  3.6192553616944125  -1.4288704858541066 
C  0.015060032638044785  4.936261537357243  -1.6414467042889873 
C  0.40847756053874695  5.895777350668212  -0.7702494116240526 
C  1.1549917870548438  5.569822618263297  0.30257547717744276 
C  1.476161643432843  4.228942220786559  0.5593907084562129 
C  2.6768816130453166  1.029746581476576  -1.1919381687093293 
C  3.0667915878720358  -0.3074613338085077  -1.3346802537072318 
C  4.092288698608236  -0.629640451837716  -2.213173100449691 
C  4.692371328106385  0.345163891254499  -2.9835400912165544 
C  4.322511766146632  1.6473626994485384  -2.8439210439308393 
C  3.3257035702945914  1.9937059137934434  -1.9476438610462818 
H  3.0285028497577984  -2.2056478214354502  -0.5316835409343126 
H  0.7165494952002875  -1.7810965863738302  2.481127038102848 
H  0.18577360161357293  -3.705415654953302  3.6438432386817823 
H  0.10981262884480736  -5.732981824385185  2.563437551362723 
H  0.6667546067218166  -5.871064727834824  0.33184610464728886 
H  1.3077583327835873  -3.9671565386908765  -0.8223494769096611 
H  0.45181881122990486  1.0467383983689342  2.627786384990379 
H  1.4016557687892477  0.6684470623782084  4.694663459024738 
H  3.6723439505479174  0.790182645960599  4.942951032035434 
H  5.013527828646352  1.291853343724779  3.1458184920891603 
H  4.0748829258387635  1.6676532782048565  1.0689694158011784 
H  0.08237121109025569  2.9718476049340565  -2.0440244077911163 
H  -0.48265092525884734  5.170817602300167  -2.3914948326033403 
H  0.16226048214694888  6.7807618881163805  -0.9132263985340561 
H  1.4593947083197425  6.236613663289843  0.874280115778364 
H  1.9707451546649786  4.008301307023668  1.3159409629609464 
H  4.3779623522076045  -1.511105579093957  -2.283213325109069 
H  5.350539380106075  0.11243491244793358  -3.5992595466104738 
H  4.739470230754648  2.3053926401867  -3.35064925104374 
H  3.0893546690380753  2.887563186102319  -1.8500523469921752 
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
