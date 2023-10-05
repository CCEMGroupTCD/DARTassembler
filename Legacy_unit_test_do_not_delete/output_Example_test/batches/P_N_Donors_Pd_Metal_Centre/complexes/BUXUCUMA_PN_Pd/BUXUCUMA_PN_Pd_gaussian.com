%chk=BUXUCUMA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.4356394176463685  -1.5591790989171193  -2.2136924238858786e-15 
C  2.133114183250996  -1.9274494890209832  1.6453574630205978 
N  1.4356394176463674  1.5591790989171193  -3.019666717411461e-16 
C  1.1562791818156921  2.8836958727978708  0.5424603067000974 
C  2.851165744762623  -1.0486719252686145  -1.0325822354148089 
C  1.3291031255245696  -1.7412980711360344  2.7791690032128384 
H  0.4522524924121121  -1.4485496560777447  2.676949715937059 
C  0.8156308395122669  3.7573244642445065  -1.8643990683139573 
H  0.6817493065603821  2.8027410741556458  -2.0413509727990706 
C  0.9452078027049065  -3.1532799316734215  -0.7273585264995477 
C  -0.05566920888761917  -5.498664199933413  -1.843905794582829 
H  -0.4178195727184997  -6.274953473462431  -2.20819290357921 
C  3.4496564449349436  -2.370904747291319  1.8189569938740473 
H  3.995940129120501  -2.508282271880228  1.0788273281191219 
C  1.166218346234295  3.0429697241293017  1.9396493179984122 
C  3.245703603638078  0.302395970948725  -1.1574262489870573 
C  1.3744389439846305  1.846960194951329  2.85950784450561 
H  0.9424611886948815  1.0763551529555284  2.436062938779386 
C  0.9007078232329316  3.9446259296135766  -0.34679672319463845 
C  4.35900944834141  0.6178860005044305  -1.9432941688058842 
H  4.626391576976743  1.5063254401827033  -2.0177579137466872 
C  0.545610534021572  -5.539209907635565  -0.5951024768605913 
H  0.6170224444079708  -6.344316758972808  -0.13642263537548782 
C  1.0490931907466772  -4.356636724554427  -0.029285102583005783 
H  1.449491640143446  -4.3765974365865565  0.8095918023892628 
C  3.6087638293057562  -2.0231251173715727  -1.6884111831572477 
H  3.3764951546950748  -2.91820762006114  -1.5968604531023796 
C  -0.1230109522509566  -4.3028741965304045  -2.557106535862748 
H  -0.4929289415557023  -4.290481078887011  -3.409483554454766 
C  0.3630628759020418  -3.1326938456860645  -1.9962257998033532 
H  0.301718915429114  -2.332025506258586  -2.4650148448459515 
C  0.9977168884803084  4.334318983748275  2.4418814340725628 
H  1.017209100391598  4.476486009431179  3.3616736336933744 
C  0.7422200435348436  2.0160040592405917  4.247911836094616 
H  0.9080980364269302  1.230855620648629  4.773664630524345 
H  -0.2048926063754961  2.1448474006755305  4.15474273312737 
H  1.1273524073442074  2.780427072496593  4.683858414417249 
C  2.859609388639707  1.5090624124798335  3.0070926639280877 
H  2.958301851683134  0.750838789663785  3.5874811506696114 
H  3.322955586807414  2.26179914157986  3.3810405577680025 
H  3.2264557812465844  1.3026631186312265  2.1446445868068293 
C  1.8257959946745879  -1.9878919295754094  4.047984952322686 
H  1.2785432113055837  -1.8658425190027785  4.790733733859079 
C  4.697643968022094  -1.681906135299853  -2.472845265654346 
H  5.179318401611387  -2.346858981781739  -2.909283682824622 
C  0.7260977848966259  5.208785926607464  0.21964410911939033 
H  0.5566022639100293  5.935740470051002  -0.3363511068063284 
C  5.069765669605989  -0.36050555891262526  -2.610906541665463 
H  5.7943218473498  -0.1297673100162275  -3.148741112690348 
C  -0.36964213564982895  4.514879206236843  -2.4646060636375196 
H  -0.3922699329434032  4.372958893712334  -3.41395291157173 
H  -0.2759814765845148  5.452345714103171  -2.2846369341392445 
H  -1.1852297330520458  4.19397320169413  -2.0733744478239045 
C  3.9438140122302547  -2.6104576730965765  3.108645353293868 
H  4.820340244413495  -2.9036798727255513  3.2219651041733472 
C  3.1348335537898238  -2.4137332523978716  4.216245127653865 
H  3.469561446905826  -2.567142447124842  5.069705259437049 
C  0.7980993896393993  5.411660631777441  1.5869236061194372 
H  0.7128521134060763  6.269905406819242  1.9336477587368337 
C  2.109926217517618  4.184676992547876  -2.5680873665161794 
H  2.0194829832732943  4.049173104542494  -3.514768504869995 
H  2.842704716683815  3.6568934713529577  -2.238468850781923 
H  2.278812296168932  5.113021985759554  -2.3924764727293746 
C  2.614944402393224  1.448642573562401  -0.49652366460598263 
H  3.1470992962127244  2.2068702690343773  -0.4262474893883235 
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