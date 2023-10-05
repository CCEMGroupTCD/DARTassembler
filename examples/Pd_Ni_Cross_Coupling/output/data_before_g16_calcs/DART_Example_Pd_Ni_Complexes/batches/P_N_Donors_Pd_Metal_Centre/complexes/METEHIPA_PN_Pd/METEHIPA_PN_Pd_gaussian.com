%chk=METEHIPA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.4745052721167187  1.522476338896602  2.733416555014215e-15 
P  2.9207766273154165  -1.571043438503278  0.6831828042095158 
N  1.4745052721167182  -1.522476338896602  -5.551115123125783e-17 
C  0.892799112568772  3.200264699043642  0.40208419827083713 
C  -0.03552958909350168  3.3351942541058137  1.4261793778517204 
H  -0.3629839265211021  2.563344032215561  1.873025595098578 
C  -0.4843311514015154  4.585967505124019  1.7996390305683545 
H  -1.1221805027960854  4.671216790570428  2.4984524278758276 
C  -0.008318162745101443  5.711417337232293  1.1597470083975394 
H  -0.3043821425618205  6.572405184143603  1.430856870494287 
C  0.8909887176835657  5.586392458265731  0.13346979221611724 
H  1.2058237451730187  6.363014580260039  -0.31407784222605406 
C  1.3465710005061633  4.333657475610918  -0.2585226271123688 
H  1.9660503739356678  4.25304127037718  -0.9745083415991309 
C  2.075383783126016  1.5760738680766875  -1.7159118974682637 
C  1.2198587862610362  1.1962902187497928  -2.7384339360782683 
H  0.33573453647125007  0.914882163954994  -2.532868964280364 
C  1.64053898896045  1.223009353202472  -4.060575492384141 
H  1.0443985040992438  0.974876344269282  -4.757513618896972 
C  2.9354226169157576  1.614523345749638  -4.359288119658434 
H  3.2269406422874307  1.6314293457592361  -5.263033530737851 
C  3.792726282457936  1.9750278842949003  -3.3647003292968956 
H  4.678542219719409  2.240978211177122  -3.5789771804143613 
C  3.3768771365618546  1.955702427965429  -2.034536029446982 
H  3.981581934653556  2.20192664204529  -1.3440909023256096 
C  2.9878863604266632  1.3002839104446982  0.9930679263195135 
C  3.6833290428401955  2.4238656051459495  1.4471526770101215 
H  3.34702812890892  3.2891566289491565  1.2454293530869978 
C  4.841304792402244  2.3100510744547327  2.1797887572230557 
H  5.302325342302363  3.091542125975636  2.461484193754342 
C  5.332169700779286  1.0682802983189488  2.503843384190691 
H  6.123592730114518  0.9882954395143417  3.0238017584860675 
C  4.668362679974029  -0.06327620668956824  2.0688268613980925 
H  5.015090627886745  -0.9195032899588947  2.286789087511583 
C  3.495785561129215  0.02798835122185439  1.314196719652576 
C  4.217348102972826  -2.00495206691505  -0.512191279791588 
C  4.064431434085195  -1.5464331644105096  -1.8096940614542627 
H  3.2541852681705796  -1.125425146995579  -2.074007946623114 
C  5.098238478749678  -1.704294442027594  -2.7273358932335303 
H  4.994392585058577  -1.3898148387866036  -3.618262483648756 
C  6.272384795823277  -2.3163629229573903  -2.3465011117930783 
H  6.9793903847109515  -2.4137043634867092  -2.9742106793777445 
C  6.4271016023477845  -2.78417378595982  -1.075825163238205 
H  7.232140625114123  -3.2240792005066172  -0.8262185566032538 
C  5.403948384061685  -2.617003946601624  -0.14123822414705503 
H  5.521240802039076  -2.9232287089108286  0.7508794102676482 
C  3.0562240024942873  -2.7362727065559502  2.0661256994263226 
C  3.264063243093564  -4.089735323579876  1.8187319962205517 
H  3.4152347694132796  -4.389003495018397  0.9291926934384638 
C  3.2532695270512804  -5.002536182611409  2.857521376915843 
H  3.3923133634587543  -5.9269093035887455  2.6863914707868743 
C  3.0363745461265212  -4.554153979512775  4.145399797883476 
H  3.0470858753386634  -5.175217979551871  4.864109082143795 
C  2.8076959593405624  -3.231746771227921  4.398674252619523 
H  2.6452460554209143  -2.9405842988303093  5.28715249848932 
C  2.8108552655975014  -2.322028810682288  3.3652877900390066 
H  2.643826221053225  -1.4036147642522607  3.5440706404432127 
C  -0.0708896797734575  -3.4409971200582277  0.34176287821147117 
H  -0.875085159072464  -2.881534070788568  0.2947919555364238 
H  0.25048339427151256  -3.473544476449928  1.2661639961634312 
H  -0.28479901421740017  -4.347056035233056  0.038641650202015716 
C  1.01110435908859  -2.8489235475695454  -0.549977714465292 
H  1.7902124624600864  -3.475416925910828  -0.5587035471142427 
C  0.5141468772499065  -2.676360323900398  -1.9830407364854443 
H  1.2456691393457606  -2.3445247629181187  -2.5443899687159384 
H  -0.22513457265448977  -2.0337819435302826  -1.996924644731616 
H  0.20352824304311756  -3.539010758226258  -2.3281974867322193 
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

-Pd 0
lanl2dz


