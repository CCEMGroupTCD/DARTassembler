%chk=OMEMIBID_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.2724242000999508  -1.5687691528711287  -1.6548682850569259e-15 
N  1.272424200099951  1.5687691528711292  -3.586522658602468e-16 
C  2.563693581820453  1.4779289424174835  -0.05712127445988356 
H  3.037171066697476  2.299706797776616  -0.11159124559882905 
C  3.376751328147245  0.27324353568240656  -0.04822285114066625 
C  2.992724650834126  -1.0215183203131502  -0.05630724170517819 
C  4.211991984255578  -1.9212830291083605  0.02233090716748072 
H  4.304394107703291  -2.312798453080019  0.9268810386553251 
H  4.153207941777085  -2.6543211137082205  -0.6404844916247202 
C  5.364897474128425  -1.010770081817295  -0.287338455599311 
H  5.628890997369309  -1.0995843664078626  -1.237378850813924 
H  6.1445006101467055  -1.2364872654780883  0.2794646524644353 
C  4.890607373269458  0.4148467407132028  0.002617654088457604 
H  5.193128239416993  0.7196956976962869  0.8944988115145508 
H  5.2154530465181566  1.0483817564808438  -0.6852483447571908 
C  0.985417099619455  -2.595318117741227  -1.509518090980687 
H  0.026455972630045732  -2.877442348186469  -1.4798117697865194 
C  1.132456372739754  -1.695890918769583  -2.742645422569092 
H  0.5397446905550543  -0.9081883541580827  -2.6478512164895016 
H  2.0654188464973333  -1.3714885080843664  -2.8043997625121015 
C  0.7702465004344791  -2.4546992927240066  -4.023324226961988 
H  -0.1943074546127168  -2.6772610822396885  -4.0105865799891545 
H  0.9392183059211048  -1.8737499648181004  -4.806851570043442 
C  1.579104950102692  -3.7327448619720807  -4.161531250393201 
H  2.5350106798328387  -3.5050226764955354  -4.280929332919222 
H  1.2790561175100452  -4.225949453982284  -4.965646390457015 
C  1.4206037765705737  -4.622001048866645  -2.9341288949861934 
H  1.9905794657689344  -5.425160451060412  -3.035388454488515 
H  0.47950600651274267  -4.921927382284709  -2.8657173337219826 
C  1.8103287045378438  -3.8831604311258574  -1.6559037457666241 
H  2.774023857454031  -3.658440822536664  -1.6833609809044363 
H  1.6529755882829842  -4.468597353645103  -0.8732175656432638 
C  1.3090329148405613  -2.7217379256760235  1.4379160688099497 
H  2.038945186605299  -3.3838218614157527  1.2680419166777304 
C  1.6567827409848348  -1.9860389365113562  2.739540598295626 
H  2.516988814209822  -1.5083264900188784  2.6301292640881972 
H  0.955936468784775  -1.3157317003731008  2.938429315733252 
C  1.764838446115555  -2.9757356086816698  3.90414622325357 
H  2.5286374022000047  -3.585584140729261  3.746478431223219 
H  1.9390519409856908  -2.478200940453707  4.742174756250174 
C  0.4839401288936168  -3.798058850210867  4.060232344926732 
H  0.609215832427237  -4.466591468971605  4.7795731582809315 
H  -0.2592517992098693  -3.200023050435634  4.3249393342357365 
C  0.11952561685946961  -4.513098811920356  2.7649191408801457 
H  -0.7421547766161634  -4.986487339888821  2.880258384364556 
H  0.8136907952589937  -5.186420953082301  2.552886857209321 
C  0.003854822131156288  -3.5165341694178367  1.6066642161935307 
H  -0.7432899366759214  -2.891167028430591  1.7824720327449641 
H  -0.19577897043042913  -4.005363324807913  0.7691856936581078 
C  0.7778044194310565  2.943643847501984  -0.08216457043502522 
C  0.6851629981030963  3.522790427696525  -1.3471595963949465 
C  0.15438653349167208  4.816037787441532  -1.4357710579965204 
H  0.08112686484818954  5.2375264900909  -2.2838965584321818 
C  -0.26245686036578175  5.481703987015118  -0.3049722183514317 
H  -0.6193869588647993  6.358605739884432  -0.3815322640685867 
C  -0.16747960817512175  4.891941115717762  0.933655626162592 
H  -0.4659812293293051  5.362774181658773  1.702798503776978 
C  0.3676636980937755  3.5973433947700735  1.073287938585773 
C  1.1097716340475523  2.7749596950148656  -2.584830651860926 
H  0.8877105981428947  3.304704248876749  -3.3789894609077296 
H  0.6407441104426234  1.9153542058758426  -2.624162002499945 
H  2.0767485208517975  2.6184845759377606  -2.55658304140964 
C  0.4972990611910112  2.9564207496760972  2.423934678351754 
H  -0.03922444237118117  3.4544131525546216  3.074928777250229 
H  1.4377553399113152  2.9652106853876474  2.700579730441565 
H  0.17947300043677106  2.030426171479969  2.378053023446754 
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

