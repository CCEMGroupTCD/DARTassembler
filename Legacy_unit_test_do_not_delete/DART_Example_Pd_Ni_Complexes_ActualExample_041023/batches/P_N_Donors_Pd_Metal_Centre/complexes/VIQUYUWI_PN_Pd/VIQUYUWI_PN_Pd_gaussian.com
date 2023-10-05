%chk=VIQUYUWI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.4356394176463652  -1.5591790989171193  2.117546749519919e-15 
C  2.1331141832509997  -1.9274494890209832  1.6453574630205978 
N  1.4356394176463674  1.5591790989171193  -3.019666717411461e-16 
C  1.1562791818156906  2.8836958727978708  0.5424603067000967 
C  2.8511657447626226  -1.0486719252686145  -1.0325822354148053 
C  1.3291031255245676  -1.7412980711360344  2.7791690032128393 
H  0.45225249241211163  -1.4485496560777447  2.676949715937056 
C  0.815630839512269  3.7573244642445065  -1.8643990683139622 
H  0.6817493065603815  2.8027410741556458  -2.041350972799074 
C  0.9452078027049083  -3.1532799316734215  -0.7273585264995488 
C  -0.055669208887607624  -5.498664199933413  -1.8439057945828325 
H  -0.41781957271849546  -6.274953473462431  -2.2081929035792127 
C  3.449656444934944  -2.370904747291319  1.8189569938740522 
H  3.9959401291204864  -2.508282271880228  1.0788273281191152 
C  1.1662183462342894  3.0429697241293017  1.9396493179984144 
C  3.245703603638077  0.302395970948725  -1.1574262489870497 
C  1.3744389439846263  1.846960194951329  2.8595078445056004 
H  0.94246118869488  1.0763551529555284  2.4360629387793895 
C  0.9007078232329329  3.9446259296135766  -0.34679672319464167 
C  4.359009448341418  0.6178860005044305  -1.9432941688058842 
H  4.626391576976741  1.5063254401827033  -2.017757913746682 
C  0.545610534021575  -5.539209907635565  -0.5951024768605938 
H  0.6170224444079684  -6.344316758972808  -0.1364226353754877 
C  1.049093190746675  -4.356636724554427  -0.029285102583007365 
H  1.4494916401434454  -4.3765974365865565  0.8095918023892629 
C  3.60876382930576  -2.0231251173715727  -1.6884111831572455 
H  3.3764951546950743  -2.91820762006114  -1.5968604531023738 
C  -0.12301095225094505  -4.3028741965304045  -2.557106535862743 
H  -0.4929289415557019  -4.290481078887011  -3.4094835544547704 
C  0.3630628759020462  -3.1326938456860645  -1.9962257998033532 
H  0.30171891542911755  -2.332025506258586  -2.465014844845953 
C  0.9977168884803022  4.334318983748275  2.441881434072555 
H  1.017209100391592  4.476486009431179  3.361673633693368 
C  0.7422200435348347  2.0160040592405917  4.247911836094623 
H  0.9080980364269311  1.230855620648629  4.773664630524357 
H  -0.20489260637549078  2.1448474006755305  4.154742733127345 
H  1.1273524073441998  2.780427072496593  4.683858414417253 
C  2.8596093886396936  1.5090624124798335  3.0070926639280735 
H  2.9583018516831343  0.750838789663785  3.5874811506696225 
H  3.322955586807405  2.26179914157986  3.381040557767996 
H  3.226455781246573  1.3026631186312265  2.144644586806825 
C  1.8257959946745892  -1.9878919295754094  4.047984952322693 
H  1.2785432113055752  -1.8658425190027785  4.790733733859096 
C  4.697643968022101  -1.681906135299853  -2.472845265654347 
H  5.179318401611397  -2.346858981781739  -2.9092836828246242 
C  0.7260977848966261  5.208785926607464  0.219644109119393 
H  0.5566022639100315  5.935740470051002  -0.3363511068063301 
C  5.069765669605966  -0.36050555891262526  -2.610906541665435 
H  5.794321847349796  -0.1297673100162275  -3.148741112690343 
C  -0.36964213564982407  4.514879206236843  -2.464606063637519 
H  -0.39226993294339074  4.372958893712334  -3.413952911571728 
H  -0.2759814765845128  5.452345714103171  -2.28463693413925 
H  -1.1852297330520472  4.19397320169413  -2.073374447823907 
C  3.9438140122302485  -2.6104576730965765  3.108645353293873 
H  4.820340244413498  -2.9036798727255513  3.2219651041733584 
C  3.1348335537898153  -2.4137332523978716  4.216245127653867 
H  3.4695614469058134  -2.567142447124842  5.069705259437064 
C  0.7980993896393958  5.411660631777441  1.5869236061194385 
H  0.7128521134060797  6.269905406819242  1.933647758736839 
C  2.109926217517624  4.184676992547876  -2.5680873665161776 
H  2.0194829832732997  4.049173104542494  -3.5147685048699917 
H  2.8427047166838193  3.6568934713529577  -2.2384688507819184 
H  2.2788122961689314  5.113021985759554  -2.3924764727293755 
C  2.6149444023932285  1.448642573562401  -0.49652366460598424 
H  3.147099296212721  2.2068702690343773  -0.42624748938832485 
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

-Pd 0
lanl2dz


