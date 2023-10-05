%chk=ACUJUVAQ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.6009657858617716  1.388887523343773  0.0 
N  1.6009657858617716  -1.388887523343773  0.0 
N  2.98591610100039  0.4353732050556609  0.3184595883578246 
C  2.7828043525564663  -0.9322763152070668  0.2963040615263258 
C  1.449057727252656  -2.8470637694116983  -0.250986368420413 
H  2.367370197165411  -3.239712379060861  -0.30077659781665067 
C  0.7830574697553234  -3.069386970398198  -1.6080696051147076 
H  1.352631726871547  -2.6869809072194313  -2.3232578380119184 
H  -0.08932893460126601  -2.60198482365192  -1.628911982876923 
C  0.569318953342256  -4.5616927008219745  -1.8670072492414393 
H  0.09002660992253975  -4.6853601376108704  -2.7249084700223145 
H  1.4466836620902936  -5.015385299689128  -1.935915329875468 
C  -0.24249963344115444  -5.184287497352499  -0.734383431954281 
H  -1.1345662228302107  -4.758938881952939  -0.6984738862393244 
H  -0.37009143295809  -6.150413718840539  -0.9130793939939221 
C  0.45639637788159115  -5.0075941936289725  0.6090938554519761 
H  1.314585480651501  -5.500499197449467  0.5999215549932106 
H  -0.1089612160774418  -5.390267506670496  1.3257237295422908 
C  0.7233735717648616  -3.525809495149394  0.9086093396894477 
H  -0.13684342638195268  -3.063073723391061  1.0731033871031161 
H  1.2731053228036548  -3.4510758840718694  1.7291109541151375 
C  4.235177319494742  1.1613056675870037  0.7122965569237985 
H  4.0490864613108934  2.117659701067099  0.48618478114299246 
C  4.493133227430142  1.1685868331457894  2.219129691079467 
H  4.80056998590833  0.27408284767619384  2.512737200376203 
H  3.657021462371993  1.384142999478992  2.7030913274589135 
C  5.561692565583806  2.2162430583935153  2.5332378007237466 
H  5.208284354928021  3.11710160811087  2.3249735694159734 
H  5.771584628415792  2.1889100783199344  3.5005527908693983 
C  6.8387527106449895  1.982948497420056  1.73643731520898 
H  7.2704330439857126  1.1484503213477228  2.049719095622076 
H  7.466449919560837  2.73044280315084  1.9011038473863149 
C  6.567147008478946  1.8755175661947734  0.23698574961339003 
H  7.408046549868308  1.6412311754461537  -0.2307957788816708 
H  6.264406059963037  2.7539389318520557  -0.10334540154452365 
C  5.507058131752516  0.8217510333394398  -0.07467380450765482 
H  5.835104012784517  -0.07647597139060225  0.18365675487632072 
H  5.312762411121964  0.8152484081461051  -1.0464053460407485 
C  3.9439856014251697  -1.8229139796753244  0.6359773573655619 
H  3.6217359167384866  -2.606209611405701  1.1266283779814834 
H  4.583210387438742  -1.3301163099603588  1.1918955075184519 
H  4.384966296494435  -2.1134626441405855  -0.1895998710956478 
C  1.4575041653188068  2.5541080453798326  1.3705312832829586 
C  0.8283539374693683  2.0983616984934477  2.535700845233828 
H  0.4892524119239534  1.2122498828749715  2.5766844467605243 
C  0.7002777704548042  2.941705428862665  3.6312422035995575 
H  0.28851111216669856  2.628049094078115  4.4275961718647485 
C  1.1765412757060234  4.245540307183318  3.558399471358033 
H  1.0721148698916765  4.828561823605688  4.30090280950357 
C  1.8016966904786957  4.699331812547713  2.4114746364017186 
H  2.1261141688282166  5.591622400641113  2.3697619747595713 
C  1.956398644523094  3.8540197118432102  1.3174606886764608 
H  2.401227403794906  4.161936510584694  0.5367724538621281 
C  1.9737878562425497  2.3501780508773944  -1.4765742939355184 
C  2.9288015731846473  1.9137825726166386  -2.3930074486205872 
H  3.3877010869088067  1.0936198734388396  -2.2503180195057215 
C  3.2062133942887554  2.6838013120932755  -3.5187139078569847 
H  3.863746138782048  2.3917141213153457  -4.137678157957591 
C  2.532987188353666  3.871904169427075  -3.743438018796434 
H  2.7495876463461215  4.405866432414759  -4.4995622435786835 
C  1.5408249522093986  4.282078870707759  -2.8615035361966514 
H  1.0587966673350566  5.083620400377366  -3.0286930059488117 
C  1.2510542626999062  3.5232301538134037  -1.7350629330714726 
H  0.5633591690736031  3.7991521910978783  -1.1400822295666424 
Br  -1.654629867976521  1.654629867976521  0.0 
C  -1.352040782177733  -1.3199206589279684  -0.4819701001757076 
C  -1.6380700204479186  -1.599217458033982  -1.8139382277462257 
H  -1.1579680913617836  -1.159339363685303  -2.503515938957912 
C  -2.63068216563725  -2.52632649780156  -2.142289956238523 
H  -2.8368222554424674  -2.694399528648736  -3.0538898079700485 
C  -3.3101544432602585  -3.1961573187854975  -1.1619874135498987 
H  -3.9960123560788676  -3.811990145619169  -1.3868663427944008 
C  -2.9974710358883443  -2.9701089012107795  0.12423560353005347 
H  -3.438659168087169  -3.468581175155742  0.8017084178296768 
C  -2.0389854250767776  -2.015940632791841  0.4912562833419403 
H  -1.8619130127668426  -1.8492165108145995  1.4101568590026907 

-Pd 0
lanl2dz
F 1 1.0
1.472 1.0
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

-Pd 0
lanl2dz
