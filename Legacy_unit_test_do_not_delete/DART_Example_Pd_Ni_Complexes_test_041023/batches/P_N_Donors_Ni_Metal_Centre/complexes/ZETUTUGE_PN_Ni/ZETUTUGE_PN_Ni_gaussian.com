%chk=ZETUTUGE_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.3816794038415714  1.4734524169446388  0.0 
N  1.3816794038415714  -1.4734524169446392  0.0 
C  1.6701195537471754  2.741669112312962  1.2470703362326925 
C  1.239961766437022  2.5167889898268294  2.5305851115769746 
H  0.7913514257915537  1.7362693057336238  2.759205943980775 
C  1.505918944260528  3.5210512469470037  3.491263355573407 
H  1.1912349165732745  3.4165047490564047  4.360155133930935 
C  2.194335386127852  4.604332978100436  3.176091285878183 
H  2.3928722215863782  5.227417927734404  3.8388638943237496 
C  2.606086454958323  4.823948831506144  1.9375967302785744 
H  3.0638159822275157  5.6067344116416  1.7323512299675412 
C  2.3564845996712176  3.8807520549982164  0.9386604086732775 
H  2.655851952982347  4.0264864455564595  0.07141087988693506 
C  1.2080353061038505  2.316494764776813  -1.5912793554503646 
C  0.36569769191569557  3.3602067349189606  -1.7166790249308383 
H  -0.07242795927351131  3.680866940546593  -0.9604700086549813 
C  0.13267191827466873  3.9713831730881446  -2.936267976277104 
H  -0.45627305521179595  4.6880907354468935  -2.996997197108616 
C  0.7560920650126128  3.5275822756312913  -4.0183586120640165 
H  0.6038061520513607  3.931275586768836  -4.842472397445844 
C  1.6263529492449817  2.472172398042807  -3.9206381757993745 
H  2.0713475006272417  2.1616374667901264  -4.677632703990927 
C  1.840682899796345  1.8738131705163408  -2.699193069734163 
H  2.4279704303945424  1.1540914490650227  -2.6363811275970512 
C  2.8774311523767993  0.48030306364931574  -0.09909489796622029 
C  4.149972853666858  1.0579042065248885  -0.19135067408179962 
H  4.245145827064391  1.9774930642429247  -0.29314015239616775 
C  5.26272730508687  0.2403252598643657  -0.12725475516011908 
H  6.113590942848036  0.6093331244871347  -0.19579186760564493 
C  5.118795469635929  -1.079348665563154  0.034329053714558855 
H  5.879479631699077  -1.609810247499238  0.10561317309717655 
C  3.873326566690822  -1.6764098396350229  0.09623798855110643 
H  3.8005060836286035  -2.598711428320162  0.1856404871484858 
C  2.7289013982963026  -0.883529905702352  0.02238120611272376 
C  1.0280200032920361  -1.8920390797422963  -1.4145373392494454 
H  1.2403761799686532  -1.1697375464448112  -2.026069005823961 
H  0.07489565853996316  -2.062267410237087  -1.4713977731020282 
C  1.7916698084379548  -3.148500166445687  -1.8302903581209253 
H  1.4173135313108105  -3.478195003831299  -2.662214839054622 
H  2.7165885346136767  -2.907849008714191  -2.0004807522446395 
C  1.7644723103012543  -4.1352416881127745  -0.9247223565256635 
H  2.6812163525451433  -4.333401589065293  -0.6766548764626056 
H  1.4096420192765797  -4.925290197051991  -1.361267964789303 
C  1.068949563270905  -3.947272835630714  0.20750728398900664 
H  1.3053519391070747  -4.657664612088961  0.8229182853852433 
H  0.12733513125483809  -4.050884043732411  -0.001537876143418604 
C  1.2341275508469796  -2.6583357697577856  0.9190592586466619 
H  2.019232033473979  -2.712174960007647  1.486367977910035 
H  0.46528832896768657  -2.5159213133512948  1.4922721628162392 
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


