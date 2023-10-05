%chk=EFELEMEJ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.521536413465021  -1.4754751582117542  1.9318691422499145e-16 
N  1.521536413465021  1.4754751582117542  -2.0844916858817484e-16 
C  2.6867098737026502  1.4442655341501665  -0.5268965998408409 
H  3.238437384323517  2.179091724862893  -0.3840295866196863 
C  3.228550872880665  0.34060970779683825  -1.3412326475300853 
C  2.8134524597708515  -1.0032742125554333  -1.2066207756452922 
C  3.4589243140562322  -1.9882390131565904  -1.9531024478262937 
H  3.198467990693815  -2.8771125110943836  -1.8632610971888124 
C  4.491640866102917  -1.6583235501338374  -2.835602561360155 
H  4.894835788956428  -2.3215972247554397  -3.3481934086634304 
C  4.9158830998973375  -0.343796816691102  -2.9487249911850264 
H  5.607960885605101  -0.12380704716262246  -3.528344711202104 
C  4.299338973182595  0.6430721061071407  -2.1891769156634426 
H  4.603851983051542  1.5204553800952836  -2.245272257238349 
C  1.2575508110783606  2.5393929417565118  1.0175634238659454 
H  0.2901476079032499  2.555130539825017  1.177552362526078 
C  1.675069545716415  3.9557274007742795  0.5987242107720757 
H  1.2340968072601117  4.198606813555671  -0.23023021185038337 
H  2.633490803188229  3.987605250250923  0.45384896810128755 
C  1.280016087307438  4.939586604991167  1.7084894838261817 
H  1.537080300850847  5.838664198820514  1.450520473304871 
H  0.3170445542781597  4.924772617864138  1.826474959036755 
C  1.9579292384964984  4.575705354932908  3.0251943920349142 
H  1.6750086592247855  5.195129524776394  3.7159766903303235 
H  2.9191792293436114  4.654469368242931  2.9266084127896987 
C  1.604296106011537  3.146833846479202  3.4453267520574533 
H  2.099966055529917  2.9157544832636075  4.24678510346464 
H  0.6579436876733722  3.098096129954611  3.6541560251336254 
C  1.932697726673041  2.1423330036481403  2.3382661564717417 
H  1.6305001265194219  1.259819480629459  2.603162574247728 
H  2.8943118765662175  2.105140784792684  2.2108949447248567 
C  1.16536235524223  -3.203641739538945  -0.4255167534104976 
C  0.38150463799382694  -3.441860206345252  -1.560116692235907 
H  -0.008383592420046382  -2.730323353856049  -2.0162156237568025 
C  0.18678765891383486  -4.747571650741873  -2.00530776585448 
H  -0.3248159741871586  -4.906998746949371  -2.7652981431571964 
C  0.7563148833096447  -5.812175226923916  -1.3151907784875256 
H  0.6340365643905335  -6.681205713722516  -1.621818069733619 
C  1.5060893540087843  -5.582261603429812  -0.17110970670123954 
H  1.8696035449654915  -6.298028313646722  0.2996677242768855 
C  1.7152520687756239  -4.27609621713097  0.27443443520410543 
H  2.221911409896073  -4.121397967059014  1.0381933654170918 
C  2.3620355002057494  -1.4882304051539046  1.6118482183082568 
C  3.711281565803312  -1.1629876809174535  1.761016696214251 
H  4.23419990856387  -0.9983619780392078  1.0094604676082868 
C  4.277492457955057  -1.0842713064982696  3.031690005136356 
H  5.173760553569668  -0.8544168503845958  3.12738655569088 
C  3.5026546942550567  -1.3492678520679555  4.15533827785363 
H  3.8832397090435116  -1.3020886470419828  5.00348256931093 
C  2.1597440085399313  -1.6843769403153357  4.0151752329161035 
H  1.6468546743968808  -1.8689604004191842  4.76852503691048 
C  1.584988032867497  -1.7432133290656648  2.753722971517497 
H  0.6832097863148842  -1.9519324649223757  2.6631641708864477 
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