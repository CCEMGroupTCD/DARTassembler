%chk=NOSIGULO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
C  3.2803541224828754  0.9539552521268835  1.9691261245271825 
C  2.5006249290874942  1.6478390128220908  2.8795647504710047 
H  1.716716713499585  2.0608156774732254  2.596612309259745 
C  2.8803741080481386  1.729546365361125  4.211561912487436 
H  2.351005736331118  2.1959762074170928  4.817981618190788 
C  4.0214586981225615  1.130292848354963  4.630351394984842 
H  4.285592458332304  1.1982417065752675  5.519022109807851 
C  4.789351633246509  0.4188945503844428  3.734638052541283 
H  5.558435437794439  -0.010203755821133864  4.03089423423817 
C  4.431390047839809  0.3335601654587017  2.405808760975774 
H  4.963799797029993  -0.13767089157493984  1.807988378373927 
C  2.8351854986140963  -0.6327086287124571  -0.5073697942182547 
H  3.621619471486438  -1.1154367669902383  -0.21143419171816047 
H  2.869722276957115  -0.5569438508043095  -1.4731012006508322 
C  4.053190685719395  1.9992725518668628  -0.5779411712836837 
C  4.06807471317352  3.3703839900680284  -0.3579200364791473 
H  3.435406019959637  3.7579366302233628  0.20171911130777959 
C  5.025479532733322  4.158446985388203  -0.9717194637277236 
H  5.040460682609952  5.074231357358906  -0.8146182250228773 
C  5.952910153734262  3.5945037278374685  -1.8121995024513995 
H  6.584153065193081  4.134902048008168  -2.23022536563416 
C  5.95800193102745  2.234923327891105  -2.041435428705717 
H  6.586574659468157  1.8600748657558754  -2.615731210942977 
C  5.01433008681765  1.4213830262013922  -1.4066330909876603 
H  5.028276738685799  0.5005655865197496  -1.5363002174300089 
C  1.815194585165235  -2.4283840954407028  1.5126605309680228 
C  3.120939037274306  -2.885048833255502  1.6848464857464416 
H  3.742248926568909  -2.7764787152765695  1.002375394089969 
C  3.4968683795854205  -3.4992964365676107  2.8689631951188286 
H  4.378771711695957  -3.7650807746541983  2.9885925328259475 
C  2.578071152616806  -3.7190941051633843  3.8675479940191124 
H  2.8377360787815036  -4.118980857099203  4.667130542444288 
C  1.2664934483589507  -3.3397440318635607  3.6731041008746455 
H  0.627189189334415  -3.529478555199906  4.321630650097921 
C  0.8986378907308408  -2.678625697086875  2.5176606984770067 
H  0.01908256634478578  -2.397341997589807  2.4131295942954814 
C  1.1976042093771757  -2.855095790833143  -1.2571399341827556 
C  1.7654339941746164  -2.751423016805064  -2.5021087897860887 
H  2.2528127253891723  -1.9910433538509738  -2.726205087224325 
C  1.618581021242708  -3.7708605316107358  -3.4287786272945193 
H  2.0182543719591077  -3.694801070723107  -4.265814242508504 
C  0.89582169119994  -4.878475754172209  -3.126225145102977 
H  0.7992574918131747  -5.556102866495534  -3.7555450574219438 
C  0.3086880496505888  -4.999456475194032  -1.8927756681879024 
H  -0.192278155651467  -5.755972097558779  -1.6855950362182515 
C  0.4639055342432379  -3.9849599335579033  -0.9445846381285072 
H  0.07488230771001136  -4.0692642022580685  -0.10449484560831679 
N  1.313790443716196  1.5342928892489858  1.2633392800161238e-15 
P  2.7986286019053805  0.9982861425824163  0.23580379876781263 
P  1.3137904437161962  -1.5342928892489858  0.0 
H  1.20989726018007  2.241176755816914  0.08370075029953987 
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


