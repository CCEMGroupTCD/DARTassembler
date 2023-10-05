%chk=SAWAQUNI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.42457821652586  1.4320184722970575  4.917138520426784e-17 
N  1.42457821652586  -1.432018472297058  -5.551115123125783e-17 
C  1.6832901197592478  2.420523161576204  1.5081094789167455 
C  1.619120654494044  3.7184131371284534  1.1118627867052415 
H  1.6836626991930042  4.430193557366329  1.7376819761138795 
C  1.453495023351802  3.9661262161580755  -0.2827456998552102 
H  1.429275963611751  4.847757591328008  -0.6361938791704917 
C  1.3308793340397966  2.866211260121627  -1.0790150401050358 
C  2.893678571685479  0.46601488940958746  -0.495371985803105 
H  3.0816629940514533  0.6183240908754997  -1.4560294034586374 
H  3.6795418143064476  0.7714678451234598  0.022428589460145892 
C  2.6786845131886414  -0.9924348445871087  -0.2571140620817159 
C  3.741825364007264  -1.8724210838558115  -0.32946280096808706 
H  4.618148887516898  -1.5462553192125819  -0.5010761341283705 
C  3.525510847859599  -3.232909476770799  -0.15208658258635696 
H  4.2441738481887565  -3.852631950445639  -0.20305605953343978 
C  2.239201634198604  -3.660319501041095  0.10166957596776749 
H  2.0652371917638463  -4.586541250731526  0.22506687402879905 
C  1.2096171367782391  -2.756991628513723  0.17745478156135697 
H  0.32992913250208944  -3.065757345265091  0.3571735821685286 
C  1.7230888132160904  2.681698559963394  3.9903016079181706 
H  1.6147211612100734  3.618795567411513  3.8778659346179047 
C  1.7949146882303317  2.1664887080845068  5.252573394105827 
H  1.7393449263542238  2.746205570024391  6.003161251420087 
C  1.94736952829129  0.8040499073682004  5.447939271486473 
H  1.9840233911653276  0.44666093166661236  6.3269183749828874 
C  2.04360210411831  -0.03118421016480699  4.35038582847853 
H  2.140725651469597  -0.9682858334751703  4.471885722054532 
C  1.9963094436998605  0.5078964266664339  3.060865611008409 
H  2.098679266049075  -0.06825421730993098  2.312620066472519 
C  1.805278691087294  1.8580189372362115  2.8519565446445014 
C  1.2047794828190437  2.8053032713719492  -2.552609090184185 
C  0.8302167612325828  1.6076625054334006  -3.1974136523267758 
H  0.6629392334589654  0.8222119880278711  -2.6892807628093234 
C  0.709101348155223  1.5877314147721056  -4.575420364335137 
H  0.464607946816058  0.7783436817069103  -5.009043477493552 
C  0.934368232435172  2.715564418496772  -5.332651482657827 
H  0.8453771191843285  2.687010485156505  -6.279233420290897 
C  1.2916154787208303  3.8883791673917187  -4.696919031232081 
H  1.4367764637691884  4.673538791203307  -5.21110451678745 
C  1.4392641005610638  3.937135828950704  -3.338488715602188 
H  1.7039919266489851  4.7505073549001136  -2.9234666091192913 
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


