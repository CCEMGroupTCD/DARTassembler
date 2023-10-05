%chk=ISULAYUB_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
C  3.245688615857812  0.953955252126883  2.0031474856873754 
C  2.450188834059115  1.6478390128220903  2.8998392966245063 
H  1.6713382124958254  2.060815677473225  2.6032488656525885 
C  2.806633619535628  1.7295463653611245  4.238261126036936 
H  2.266762390010161  2.1959762074170923  4.8353497191036 
C  3.9402355326398104  1.1302928483549626  4.676901496913005 
H  4.18881962149451  1.198241706575267  5.57004663261876 
C  4.723643867253989  0.41889455038444257  3.7947261555540854 
H  5.487440153320813  -0.010203755821134086  4.104359579140855 
C  4.38892807001274  0.3335601654587015  2.45985196004856 
H  4.931690134912998  -0.13767089157494006  1.8714144597206306 
C  2.8438086066408905  -0.6327086287124573  -0.4807405144913429 
H  3.624958013395891  -1.1154367669902383  -0.1711248190862501 
H  2.8951944618795955  -0.5569438508043096  -1.4457220854604185 
C  4.062859926012036  1.9992725518668624  -0.5300440216261434 
C  4.073901788287287  3.3703839900680275  -0.309796635009581 
H  3.4315624037143877  3.757936630223362  0.23871568585035102 
C  5.041873067531089  4.158446985388203  -0.9067935595964112 
H  5.054110141036353  5.074231357358906  -0.7494547910546119 
C  5.98383083559457  3.5945037278374676  -1.730959692999142 
H  6.622273162966041  4.134902048008168  -2.1379051809667526 
C  5.992922555940431  2.2349233278911047  -1.9600718417447531 
H  6.631422393116858  1.860074865755875  -2.5233100492205747 
C  5.038315609118087  1.4213830262013918  -1.3418355320660356 
H  5.054523140235831  0.5005655865197491  -1.4712395069401274 
C  1.7887186524977956  -2.4283840954407028  1.521180854308326 
C  3.0912591741367104  -2.885048833255502  1.7161289672132662 
H  3.724385197799077  -2.7764787152765695  1.044605171921149 
C  3.4464655745069903  -3.4992964365676107  2.9066262015491318 
H  4.326146768644963  -3.765080774654198  3.0416286545190676 
C  2.5103805770338368  -3.719094105163384  3.8890236885092837 
H  2.756051315357289  -4.118980857099203  4.693016234468402 
C  1.2023961462655144  -3.3397440318635603  3.671719222942765 
H  0.5518709074537547  -3.5294785551999057  4.308989600665803 
C  0.8547618826945014  -2.6786256970868747  2.5100318354698885 
H  -0.02283516183991363  -2.397341997589807  2.3901662948765754 
C  1.2195620221833237  -2.8550957908331425  -1.2589761950657874 
C  1.8090330262239434  -2.751423016805064  -2.493845439690789 
H  2.3002485469610456  -1.9910433538509738  -2.7094016744714473 
C  1.6783750383484797  -3.7708605316107353  -3.422937078664672 
H  2.0925958026369  -3.6948010707231065  -4.252869947541332 
C  0.9504455017122185  -4.878475754172209  -3.1330435664181637 
H  0.8648791564096052  -5.556102866495532  -3.763952907921293 
C  0.34187462186934436  -4.999456475194032  -1.910028844778092 
H  -0.16263108444822705  -5.755972097558777  -1.7116228332484669 
C  0.48052025085049377  -3.984959933557903  -0.9592733101779576 
H  0.0768946859360693  -4.0692642022580685  -0.1261008588652358 
N  1.313790443716196  1.5342928892489853  1.2600575767044729e-15 
P  2.794287110126361  0.9982861425824159  0.261681883737424 
P  1.3137904437161962  -1.5342928892489858  0.0 
H  1.2084523041018074  2.2411767558169133  0.08187481620471977 
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
