%chk=VEPODOSE_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
C  3.245688615857812  -0.9539552521268833  -2.0031474856873754 
C  2.450188834059115  -1.6478390128220908  -2.8998392966245063 
H  1.6713382124958254  -2.0608156774732254  -2.603248865652588 
C  2.806633619535628  -1.729546365361125  -4.238261126036936 
H  2.266762390010161  -2.1959762074170928  -4.8353497191036 
C  3.9402355326398104  -1.1302928483549632  -4.676901496913005 
H  4.18881962149451  -1.1982417065752677  -5.57004663261876 
C  4.723643867253989  -0.418894550384443  -3.7947261555540854 
H  5.487440153320813  0.010203755821133583  -4.104359579140855 
C  4.38892807001274  -0.33356016545870176  -2.45985196004856 
H  4.931690134912998  0.13767089157493984  -1.8714144597206306 
C  2.8438086066408905  0.6327086287124574  0.48074051449134286 
H  3.624958013395891  1.1154367669902383  0.17112481908624996 
H  2.8951944618795955  0.5569438508043099  1.4457220854604185 
C  4.062859926012036  -1.9992725518668624  0.5300440216261436 
C  4.073901788287287  -3.3703839900680275  0.3097966350095814 
H  3.4315624037143877  -3.757936630223362  -0.23871568585035055 
C  5.041873067531089  -4.158446985388203  0.9067935595964117 
H  5.054110141036353  -5.074231357358906  0.7494547910546125 
C  5.98383083559457  -3.5945037278374676  1.7309596929991424 
H  6.622273162966041  -4.134902048008168  2.137905180966753 
C  5.992922555940431  -2.2349233278911043  1.9600718417447534 
H  6.631422393116858  -1.8600748657558748  2.523310049220575 
C  5.038315609118087  -1.4213830262013916  1.3418355320660358 
H  5.054523140235831  -0.5005655865197489  1.4712395069401274 
C  1.7887186524977956  2.4283840954407028  -1.5211808543083263 
C  3.0912591741367104  2.885048833255502  -1.7161289672132667 
H  3.724385197799077  2.7764787152765695  -1.0446051719211495 
C  3.4464655745069903  3.4992964365676102  -2.906626201549132 
H  4.326146768644963  3.7650807746541974  -3.041628654519068 
C  2.5103805770338368  3.7190941051633835  -3.889023688509284 
H  2.756051315357289  4.118980857099202  -4.693016234468403 
C  1.2023961462655144  3.33974403186356  -3.6717192229427655 
H  0.5518709074537547  3.5294785551999053  -4.308989600665803 
C  0.8547618826945014  2.6786256970868743  -2.510031835469889 
H  -0.02283516183991363  2.3973419975898067  -2.390166294876576 
C  1.2195620221833237  2.8550957908331425  1.258976195065787 
C  1.8090330262239434  2.7514230168050644  2.4938454396907885 
H  2.3002485469610456  1.991043353850974  2.709401674471447 
C  1.6783750383484797  3.7708605316107358  3.4229370786646713 
H  2.0925958026369  3.694801070723107  4.252869947541331 
C  0.9504455017122185  4.878475754172209  3.1330435664181633 
H  0.8648791564096052  5.556102866495533  3.763952907921292 
C  0.34187462186934436  4.999456475194032  1.9100288447780913 
H  -0.16263108444822705  5.755972097558777  1.7116228332484662 
C  0.48052025085049377  3.984959933557903  0.9592733101779571 
H  0.0768946859360693  4.0692642022580685  0.1261008588652353 
N  1.313790443716196  -1.5342928892489853  -1.0721608891271415e-15 
P  2.794287110126361  -0.9982861425824159  -0.2616818837374239 
P  1.3137904437161962  1.5342928892489858  -1.878966875773315e-16 
H  1.2084523041018074  -2.2411767558169133  -0.0818748162047195 
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

