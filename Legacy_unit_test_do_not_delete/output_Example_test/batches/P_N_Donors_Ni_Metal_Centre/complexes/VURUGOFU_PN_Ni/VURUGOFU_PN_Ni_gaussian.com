%chk=VURUGOFU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3816794038415714  -1.4734524169446388  1.8044587861071833e-16 
N  1.3816794038415714  1.4734524169446392  -1.8044587861071838e-16 
C  1.6701195537471754  -2.741669112312962  -1.247070336232692 
C  1.239961766437022  -2.51678898982683  -2.530585111576974 
H  0.7913514257915537  -1.7362693057336243  -2.759205943980775 
C  1.505918944260528  -3.521051246947004  -3.4912633555734067 
H  1.1912349165732745  -3.416504749056405  -4.360155133930935 
C  2.194335386127852  -4.604332978100436  -3.1760912858781825 
H  2.3928722215863782  -5.227417927734405  -3.838863894323749 
C  2.606086454958323  -4.823948831506144  -1.9375967302785737 
H  3.0638159822275157  -5.6067344116416  -1.7323512299675405 
C  2.3564845996712176  -3.8807520549982164  -0.938660408673277 
H  2.655851952982347  -4.0264864455564595  -0.07141087988693456 
C  1.2080353061038505  -2.316494764776813  1.5912793554503648 
C  0.36569769191569557  -3.3602067349189606  1.7166790249308388 
H  -0.07242795927351131  -3.680866940546593  0.9604700086549818 
C  0.13267191827466873  -3.971383173088144  2.9362679762771045 
H  -0.45627305521179595  -4.6880907354468935  2.9969971971086165 
C  0.7560920650126128  -3.527582275631291  4.0183586120640165 
H  0.6038061520513607  -3.9312755867688356  4.842472397445845 
C  1.6263529492449817  -2.4721723980428068  3.920638175799375 
H  2.0713475006272417  -2.161637466790126  4.677632703990927 
C  1.840682899796345  -1.8738131705163406  2.6991930697341635 
H  2.4279704303945424  -1.1540914490650225  2.6363811275970512 
C  2.8774311523767993  -0.48030306364931574  0.09909489796622034 
C  4.149972853666858  -1.0579042065248885  0.19135067408179976 
H  4.245145827064391  -1.9774930642429247  0.293140152396168 
C  5.26272730508687  -0.24032525986436568  0.1272547551601191 
H  6.113590942848036  -0.6093331244871347  0.19579186760564501 
C  5.118795469635929  1.079348665563154  -0.03432905371455899 
H  5.879479631699077  1.609810247499238  -0.10561317309717674 
C  3.873326566690822  1.6764098396350229  -0.09623798855110664 
H  3.8005060836286035  2.598711428320162  -0.1856404871484861 
C  2.7289013982963026  0.883529905702352  -0.02238120611272387 
C  1.0280200032920361  1.8920390797422966  1.4145373392494451 
H  1.2403761799686532  1.1697375464448114  2.026069005823961 
H  0.07489565853996316  2.062267410237087  1.471397773102028 
C  1.7916698084379548  3.1485001664456873  1.8302903581209249 
H  1.4173135313108105  3.4781950038312996  2.6622148390546214 
H  2.7165885346136767  2.9078490087141913  2.000480752244639 
C  1.7644723103012543  4.1352416881127745  0.9247223565256629 
H  2.6812163525451433  4.333401589065293  0.676654876462605 
H  1.4096420192765797  4.925290197051991  1.3612679647893022 
C  1.068949563270905  3.947272835630714  -0.2075072839890071 
H  1.3053519391070747  4.657664612088961  -0.8229182853852438 
H  0.12733513125483809  4.050884043732411  0.0015378761434181078 
C  1.2341275508469796  2.6583357697577856  -0.9190592586466623 
H  2.019232033473979  2.712174960007647  -1.4863679779100352 
H  0.46528832896768657  2.5159213133512948  -1.4922721628162394 
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
