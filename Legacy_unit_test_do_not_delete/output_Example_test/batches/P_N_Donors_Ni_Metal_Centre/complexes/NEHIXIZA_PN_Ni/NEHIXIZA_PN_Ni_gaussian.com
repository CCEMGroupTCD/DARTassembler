%chk=NEHIXIZA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.2090867845196254  -1.618088114875082  5.124654597549572e-15 
N  1.2090867845196231  1.6180881148750825  -1.1973593652246553e-15 
N  2.5691307196975477  1.4956651527514837  -0.16136100126940153 
C  1.004815146608064  2.9338796749437153  0.22112918938277085 
C  2.2122679490674835  3.6403323795223215  0.1598488823788402 
H  2.3329976716256886  4.576400487356443  0.2656094212715003 
C  3.1897548734739436  2.6999662919081953  -0.08586659104463132 
C  -0.3217566583675706  3.531761317548007  0.5373046591468201 
H  -0.9649724460281144  3.285322497971867  -0.15989788582139763 
H  -0.23894267289739868  4.508333350599489  0.5755628356877992 
H  -0.6338174378221215  3.196134833422991  1.4030627218179825 
C  4.663494296221817  2.86440509629342  -0.2545783080389026 
H  5.131013101380383  2.301639174197578  0.3972605246112631 
H  4.905421444702023  3.803034022640395  -0.10907451615372063 
H  4.920222947957145  2.597270846294082  -1.162446856814259 
C  3.225341535880428  0.22703299599253407  -0.4685954159461212 
H  4.2059830239495115  0.3522582467914645  -0.41819096733264133 
H  2.9999238606968266  -0.035658550958735886  -1.3963232655201145 
C  2.819746727209086  -0.8991740895472312  0.4770548744143911 
H  2.7619181413228766  -0.5481382159268071  1.4003189211994735 
H  3.512444686124465  -1.6068063745099046  0.46211141454089694 
C  1.0387685856306756  -2.9577491290512756  1.2494285766796183 
C  0.8168803228761492  -4.291109512930346  0.8821868775853544 
H  0.7403287855539331  -4.518758549230469  -0.03637388709232385 
C  0.7067042018737795  -5.286029953418396  1.8505289243126242 
H  0.5433114116540632  -6.18445887186591  1.5895897655007538 
C  0.8363269911438709  -4.967698324711165  3.197703804076381 
H  0.7564431813313557  -5.645021262148727  3.8585846143219205 
C  1.0835098122642381  -3.6540888367852915  3.5708129920818434 
H  1.1935255741037714  -3.437283171645678  4.4899278139132734 
C  1.1714789059975508  -2.649744901457975  2.6076570503367003 
H  1.3222947977662112  -1.7511153140870488  2.8746823702814064 
C  1.571204674289576  -2.5827143120834433  -1.5101406899914764 
C  2.8647843170647604  -2.8222638035089593  -1.9908016504984114 
H  3.615593664004203  -2.480782696936059  -1.520319814077354 
C  3.0525488523718014  -3.5614834195509077  -3.1591959711538027 
H  3.9311298233176557  -3.6957311765197973  -3.4973544506284835 
C  1.9635505393268105  -4.102426333569885  -3.8320520909579043 
H  2.0993551042255874  -4.607969781717114  -4.624930033763761 
C  0.6721522897668105  -3.9035839206986678  -3.3455263815889396 
H  -0.07380307939250796  -4.287137488823055  -3.7914113304282577 
C  0.4881389233421558  -3.1362292206761326  -2.1973299762076923 
H  -0.39190985121122157  -2.985361525180549  -1.8755505442003528 
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