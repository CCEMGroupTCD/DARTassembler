%chk=RENESISI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
C  2.6402806186002916  1.1357138486741667  -1.1676027728652851 
H  2.206184906981231  1.123948391871048  -2.0464734670974627 
C  3.1913907076912116  -0.28364356847199257  -1.0030844392535492 
H  3.9042319693109535  -0.4007348066240288  -1.6498051190530987 
H  3.592361056132563  -0.34425247372432977  -0.12167169242324023 
C  2.230810375520799  -1.451029907942258  -1.1566744775555644 
H  2.772426001647939  -2.2667223019053635  -1.1537879630832997 
C  3.754305256655849  2.149897960394895  -1.3159636617087171 
H  4.1778131106933545  2.2861075350258897  -0.46675521300639883 
H  3.3901501984246636  2.9801592545998767  -1.6297876361111026 
H  4.401839901296614  1.824322450672276  -1.9470912569371737 
C  1.907049374706367  -1.8315727687941585  1.3479206205404104 
H  2.4326591444442416  -1.0447646646535689  1.6015696245100641 
C  2.8539392984070906  -3.0311342596529505  1.295874083261855 
H  3.158239259272834  -3.2361670282515953  2.1833029635772685 
H  3.606754312969655  -2.8197324198401663  0.7400020512500556 
H  2.3903823609033674  -3.788975379217674  0.9318377620120912 
C  0.8306081695632584  -2.0225700671739557  2.3938338788375506 
H  0.16973031765309443  -2.637679106828474  2.065436582579228 
H  0.41479340727788705  -1.178370383252918  2.5850173957352043 
H  1.222423664181385  -2.374527144917781  3.196107437920222 
C  0.7139546006312375  3.1264288047924613  -0.6832398449047287 
C  0.878907934015439  4.352164451703374  -0.03159535463619729 
H  1.2926404336992117  4.377610415578973  0.8004781887647916 
C  0.44135524476983434  5.514090002432339  -0.6015014422481364 
H  0.5614458957305009  6.319908464740847  -0.15305141007665563 
C  -0.17468834691824808  5.507182386619113  -1.8395253911214289 
H  -0.4814383971735279  6.299653149340063  -2.2165045774182577 
C  -0.32849942519863484  4.328041562091687  -2.493576679730381 
H  -0.7174756220942207  4.3185042794797095  -3.3390401197968775 
C  0.0900737775784699  3.1321803689123384  -1.915414608539004 
H  -0.0500699271938958  2.329700731384097  -2.3641240282810156 
C  1.474945000617371  -1.4326527284952175  -2.4647663796395385 
H  0.9555619507103126  -0.627110954782118  -2.523600793653805 
H  0.8915618717213191  -2.1934269304089224  -2.5095197998683783 
H  2.098441473345229  -1.4642907015291124  -3.1951096531474485 
C  1.9052814784298944  1.9178007876969705  1.66155577231364 
C  3.1606751341814716  2.434030491653945  1.93968781263117 
H  3.7755351799384345  2.548178231852972  1.2512882865279065 
C  3.5059374262690657  2.7794073603857727  3.226486269347467 
H  4.358206425250594  3.1052109911274237  3.4048106859892227 
C  2.5956994442051102  2.6436870073245973  4.249683939288681 
H  2.819329790033936  2.911261614780052  5.111765999112218 
C  1.3743788214524155  2.1224122418930444  3.9955201232894875 
H  0.767501787259431  2.0034278841048283  4.6905414812137405 
C  1.0260667828907626  1.765746436677311  2.70800825259356 
H  0.18044187313113302  1.4150181591211957  2.545593206133778 
N  1.2887821664268937  -1.555358649154592  -1.143979985761824e-14 
P  1.2887821664268937  1.555358649154592  1.6018779660405954e-14 
H  0.817691149145559  -2.271792701587245  -0.016670991325471063 
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

