%chk=URIJIFAG_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5828346202304269  1.4095157200258526  2.436422548612611e-16 
C  3.0612587535532825  0.431289237404779  0.4869661597148438 
H  3.8510626617412287  0.7736307119582868  0.04051108760615346 
H  3.1999084434300853  0.4900118957079309  1.4463245926787083 
C  2.8109176108025182  -1.0086132171508688  0.08089846386738245 
N  1.5828346202304269  -1.409515720025853  -2.220446049250313e-16 
C  4.047071730529008  -1.8253386879238267  -0.1767041131265356 
H  4.802035938919045  -1.2264866864828514  -0.28887681221219397 
H  3.934440702612604  -2.317247628809766  -1.0033473102962527 
P  4.432977002260565  -3.0378829048614358  1.2074293534224683 
N  1.7260087942796214  2.825693111054444  0.8729115264588938 
N  1.6695918418632005  1.7449740858192668  -1.6108295562578319 
N  5.2373996292310085  -2.072471749337014  2.3361485115867433 
N  5.516105515552666  -4.032871304121203  0.3621453898739941 
C  1.295809305124454  -2.812546106209645  -0.4724785539400511 
H  2.1325983377929507  -3.185992615192421  -0.8179006576542187 
C  0.852489246314104  -3.677809066865117  0.6742061387170627 
H  1.52817001178067  -3.6729079260677153  1.3576516903138331 
H  0.7181103920547871  -4.575311302581079  0.3634035479736762 
H  0.03146485927485476  -3.3353805730621917  1.0334426798966831 
C  0.31222793423984707  -2.730436670780377  -1.6318074492334855 
H  0.6713055041868409  -2.161823353019509  -2.3171928762932925 
H  -0.5210219585137572  -2.3690617972949948  -1.3215334621840886 
H  0.1655676587846331  -3.6089446060996693  -1.9914518585942025 
C  1.8165075758319467  2.6259348635942183  2.369097996491676 
H  1.7133474942561984  1.6567155171260985  2.4749279827960953 
C  0.7488511526180728  3.1666104156100703  3.1884254337366986 
H  -0.09843449697549267  2.9555211806531574  2.7896128801753655 
H  0.8442997437858784  4.119468738449749  3.2505894348887443 
H  0.7930391383236057  2.7827761367060537  4.066331219261318 
C  3.2263084687372654  2.8889339701194006  2.8871443383669595 
H  3.86676320706124  2.497400355304608  2.28735316918518 
H  3.325008935294812  2.497165455476747  3.7579130842668453 
H  3.376247320453013  3.835730413067565  2.942908720967897 
C  1.5849784182454651  4.13866283087144  0.2559909011922869 
H  1.5527572770234066  3.9659677508935216  -0.7090606090442292 
C  0.27898406080405147  4.883493778894315  0.5743054110698054 
H  -0.46708293418703417  4.295732980465758  0.43548314046588515 
H  0.19685610532720133  5.647479284483602  5.0549376252373435e-05 
H  0.2907254910879673  5.171069230690314  1.490901564421562 
C  2.8050708469583157  5.022445386327544  0.46140324276182587 
H  3.5987980806737765  4.525441099644492  0.25192420448572533 
H  2.8398482869053447  5.312079970887204  1.374929830824645 
H  2.744681426462705  5.788133646959504  -0.11505541613975892 
C  0.4624737371681982  1.525127396209934  -2.4636141037171315 
H  -0.10202796094828726  0.9094807434829364  -1.951686489879337 
C  -0.3633396176325896  2.7793531312500326  -2.614432566100925 
H  -0.4990409929527504  3.178347670303265  -1.751703063567293 
H  -1.213633650359369  2.557761629599898  -3.001193479722926 
H  0.09765241406881953  3.3975301530600612  -3.185821143527784 
C  0.7253288042874504  0.8382918354243878  -3.7271916383900914 
H  1.2525118734373213  0.05385797683661764  -3.5616369854720293 
H  1.2039140376321187  1.4260187498746966  -4.317274605823627 
H  -0.10721836095223547  0.5845405328185262  -4.133305223581749 
C  2.985152487488037  2.224983492871208  -2.146895416786939 
H  3.520164742225221  2.5024292314634864  -1.3757878306369657 
C  2.858448799178732  3.463237029318454  -3.0717491646004573 
H  2.3594375416659843  4.149339593667479  -2.6205230172317115 
H  2.4034384814324126  3.2160836098492425  -3.8790373461001506 
H  3.733908264311962  3.794739217170201  -3.28489344521378 
C  3.785711622891326  1.1245867197358583  -2.8695042752959354 
H  3.8649459599683684  0.3585775722960576  -2.29727718238043 
H  4.663182898147786  1.4554279784565818  -3.082862548226525 
H  3.3327856889444223  0.875040403932053  -3.677610620170854 
C  6.466817989023397  -1.3493714777193944  2.0874398993675998 
H  6.640513569595977  -1.4001566651302064  1.1235638803541514 
C  6.37756257191511  0.13754968975900583  2.4504227119379545 
H  5.625011663813753  0.5304112074651317  2.001793724219284 
H  7.181540201911557  0.5837210066624192  2.1743974478080856 
H  6.26914421308186  0.23071375180833886  3.3994661551008125 
C  7.646448078475326  -1.9887402035846642  2.7820871367802145 
H  7.68967313450924  -2.9181461842259955  2.5469015554215932 
H  7.54248568140357  -1.902747838776008  3.731899451498064 
H  8.454828847938455  -1.5496729330269097  2.5067794298645314 
C  4.751156239813017  -2.0944949623874023  3.7360418054967885 
H  5.370472207015661  -1.5280511557263776  4.240905348113259 
C  3.41072260026801  -1.459053117167278  3.8869351099187224 
H  3.41492446296773  -0.600126556222082  3.460255139387286 
H  3.2123088423896453  -1.3521070981493182  4.819143595023224 
H  2.745125212040076  -2.0168577864089787  3.478713873630954 
C  4.820572229823674  -3.4606176190149416  4.389345875115037 
H  5.698649509137028  -3.8295592403190746  4.267327591144395 
H  4.170460419713509  -4.04113667486839  3.9862960517105925 
H  4.638470079662817  -3.375576314901073  5.328112173775881 
C  6.510471470945083  -3.6329396985413354  -0.652755313453626 
H  6.529681740652834  -2.6527791473878803  -0.6579388933276106 
C  6.151531362678821  -4.071145375302905  -2.062789328480159 
H  5.265878214295928  -3.770247773400161  -2.275824827701907 
H  6.186759983693072  -5.0284817787417015  -2.1203857251755824 
H  6.775710314788982  -3.6895294434209003  -2.6851086406314817 
C  7.919912387024774  -4.098569113435667  -0.29028647056421986 
H  8.13023996743886  -3.815413796804876  0.6027528709996867 
H  8.551510255625494  -3.717368419916517  -0.9046153349679309 
H  7.962478832362015  -5.056343028844882  -0.3398383090069823 
C  5.344890442433681  -5.503067290982554  0.48207090406357755 
H  6.002209325724943  -5.916936641080544  -0.11346873082079856 
C  5.632797959616018  -6.005853052738373  1.8625457641481735 
H  6.499863619047999  -5.701592625623639  2.1401112964567126 
H  5.613925175849879  -6.96610452121804  1.8647874857102242 
H  4.967524485966031  -5.673202991913652  2.468913212393071 
C  3.968547213532153  -5.976138942845911  0.022314187462013302 
H  3.7964175860219793  -5.6474399092576215  -0.8632264534186462 
H  3.2973920952459763  -5.643301008983496  0.6229014142496465 
H  3.9452241742741343  -6.935891537854348  0.017907123321683316 
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

-Pd 0
lanl2dz


