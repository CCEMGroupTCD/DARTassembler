%chk=AFOBOYAM_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.3896009157308442  -1.4659840705137286  2.295009237842896e-16 
N  1.3896009157308442  1.4659840705137286  -1.7953126995556455e-16 
C  2.5377346396227365  0.9842149816091663  0.30936814160194104 
C  3.9185621170005684  1.5331193054590242  0.10184577900614669 
C  4.642665587794815  1.5557216895274877  1.4756367368371437 
C  4.91801937067229  0.05459638109979758  1.7677068357054417 
C  4.226335809209166  -0.6588551331676668  0.6047345885769033 
C  2.7024334683202373  -0.43717522781499973  0.7920971114999591 
C  4.547463434284758  0.2695311926967489  -0.5761964205685585 
C  1.2157205397183437  2.8417070306501584  -0.40992449523200647 
C  1.1157144571817796  3.8382433535094185  0.5797252815465529 
C  0.7905637201971512  5.128478607796438  0.15408232455881388 
C  0.5979106987138214  5.420730736327317  -1.1716238558154626 
C  0.7476569109466712  4.424726816933893  -2.1326742262304927 
C  1.0478734326249437  3.133540941812698  -1.7770181339180973 
C  1.3664863477594626  3.578402307033033  2.058541777756347 
C  0.13332888240654128  3.8561333841908643  2.9217179906034874 
C  2.542700175800216  4.4174243944756135  2.5686423438305805 
C  1.238245452468167  2.0598840640483753  -2.8339778459649807 
C  0.21880949869488453  2.1205673564451515  -3.9662776226046934 
C  2.6634379211521466  2.1246469914973285  -3.396891916169957 
C  2.038057179207094  -2.0876473773193984  -1.5649644044929372 
C  3.0758940222366853  -3.045934604483948  -1.591924884898233 
C  3.5616899478675488  -3.520059800644115  -2.814332651155692 
C  3.0131140573001276  -3.051073667828177  -3.985119326649816 
C  1.9869039710811227  -2.1081801959259825  -3.9865227462759045 
C  1.5159245605174105  -1.6472333694279302  -2.7677318093361487 
C  1.156839793537236  -2.927345273605969  1.0523387021734052 
C  -0.08382310893992062  -3.5794281674296426  0.9615820936011462 
C  -0.3559710118537067  -4.67719744908087  1.7195404191063 
C  0.5988870835045139  -5.166583070268819  2.5995581640997694 
C  1.8268131306824862  -4.538617576975708  2.6876141229842534 
C  2.1103223796304165  -3.4148829620267818  1.918401407978459 
H  3.9812816503009136  2.3332113962202623  -0.40409582478218387 
H  4.065612521074554  1.912966744595837  2.1295235915895345 
H  5.475164995285185  2.0131949141614927  1.378415878626667 
H  4.4867348204551565  -0.18431740182914194  2.581817719474125 
H  5.8491853469120105  -0.11137075312338497  1.7362734927182757 
H  4.486757800004579  -1.5627904600608324  0.477529315475883 
H  2.512607081906488  -0.46982463442364564  1.7226028499396622 
H  5.484597085975102  0.3859192650038427  -0.647958217468106 
H  4.016899809540492  0.022153951501409145  -1.3270571803881266 
H  0.7030475852619154  5.816586695250963  0.8008126285539832 
H  0.3615838299485352  6.302881409729119  -1.4381745273890818 
H  0.6360756833343942  4.645154403085466  -3.0472380232987817 
H  1.5977342532495338  2.6681079888059887  2.160925128014178 
H  0.3695827803092355  3.8094059664938933  3.8343103005432564 
H  -0.5278603573891325  3.188271267069154  2.732556044926922 
H  -0.21576055342024558  4.706073363118031  2.712249405809915 
H  2.753684867460524  4.152242485054342  3.453224615188467 
H  2.288661367668875  5.339488577974081  2.5674997710702003 
H  3.2849891003024174  4.293651143012919  2.0024239453448476 
H  1.1393460456341993  1.2126071734033568  -2.4043132466473085 
H  -0.660399579755135  2.0622746664229106  -3.6017627550601 
H  0.3619054194672786  1.4027958259323676  -4.56310384307732 
H  0.3103776767542463  2.949532022189532  -4.424220579439052 
H  3.2912969146084787  2.086058134266311  -2.6796890198891394 
H  2.7814442065773397  2.947495990516314  -3.8635455044060585 
H  2.8114502875802763  1.3999108286905366  -3.9856538661161776 
H  3.441599341414688  -3.3686877261015606  -0.7776983734095557 
H  4.260896835650751  -4.1586960954245304  -2.8347498269075557 
H  3.34051232061963  -3.375773809237758  -4.81237987529458 
H  1.6234950367422585  -1.7913269883483425  -4.808276384108121 
H  0.8090553221341756  -1.0100464662491946  -2.7578677206808098 
H  -0.7426724904129975  -3.2429899141632452  0.3608018889741725 
H  -1.209875856220219  -5.111880119116085  1.654614745151516 
H  0.4100621861493592  -5.925842442446  3.13642445462277 
H  2.488558098378685  -4.8762573814972825  3.2838093679597513 
H  2.9578124430315276  -2.994703972780234  1.9851989336224487 
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

