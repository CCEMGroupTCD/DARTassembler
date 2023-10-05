%chk=EVIYABUK_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5724445483068719  1.4210975133677501  8.806317839934983e-16 
C  3.1556309928357056  0.4144431447946597  0.11666307939064147 
C  2.7691872352581357  -1.0272913637293901  -0.2345006561519899 
C  3.543964531129623  0.3972385337317226  1.6199544398375683 
N  1.572444548306871  -1.4210975133677497  2.220446049250313e-16 
C  4.307306891502943  0.9651056099941857  -0.7002099669571332 
C  1.952002976589788  2.970444851807734  0.8517063721661129 
C  1.512609807670431  3.196009021039413  2.181360243918413 
C  1.8516885573652044  4.312495757575847  2.8464324419869076 
C  2.6118491580520855  5.324289263281543  2.2514686674802764 
C  3.074767625749349  5.107835095213204  0.9284108700312456 
C  2.7025888339131607  3.966083509036075  0.2605629589247 
C  1.293956052958778  1.9660343334771535  -1.7212571423375078 
C  0.46311773297990166  3.045211452270081  -1.9089725980903511 
C  0.1746486641994469  3.4806594839350717  -3.218427380843867 
C  0.6896173194045807  2.7823054789040462  -4.307583649619549 
C  1.501521374499855  1.7109140081021414  -4.105132614520197 
C  1.7810231928678655  1.2788105358043227  -2.8219209707306003 
C  3.851631821081696  -1.8849655335416848  -0.8135995916181643 
C  3.912067442429534  -2.151437697793516  -2.213687567907968 
C  4.894907330427545  -2.8729394317386854  -2.750990548054337 
C  5.9393916326906595  -3.3747527579122143  -2.0016501861355174 
C  5.918719329217225  -3.1223736905883617  -0.6519339595731015 
C  4.940815824028347  -2.385213411194563  -0.07348669575092348 
C  1.2026818172916345  -2.834016796606549  -0.17436052249323672 
C  0.5958921988539545  -3.208106727099888  -1.3729542413403537 
C  0.1750933641525394  -4.531235184023326  -1.468212402110706 
C  0.27756447094504755  -5.427009128475074  -0.4381181429702633 
C  0.886935962189815  -4.996262062040941  0.7009362369504426 
C  1.3524889324952187  -3.6970256108951642  0.8765788452604726 
C  0.36441757443502576  -2.2836478967648377  -2.5039802088873895 
C  -0.23637665308207545  -6.862845051630317  -0.5464969361480937 
C  1.9814076131671443  -3.2520007417000385  2.1854745091312795 
H  3.0699616880166776  0.1704410747479237  1.9191920666817766 
H  4.028354647309503  0.0562420500691696  1.7154151829574922 
H  3.66107097844129  0.9556671760557474  1.8146020738591144 
H  4.163036346624047  0.9750850624009324  -1.2781257051145305 
H  4.425935189005465  1.5315633002144105  -0.5101903305454722 
H  4.791381137223825  0.6072681708202241  -0.6195455897192654 
H  1.1735831791611635  2.7671676929339473  2.4480078962058034 
H  1.672978491759031  4.390933235617618  3.4206955563251746 
H  2.7166067538326226  5.8386083604052175  2.5354874035879176 
H  3.445247077735571  5.501033475158158  0.6659869769409617 
H  2.8685933945707847  3.896821323243671  -0.31423321839395457 
H  0.24154233534155067  3.3251647322228317  -1.4273963484630328 
H  -0.18713667830356995  3.957439475192902  -3.310102818591205 
H  0.5905872172793516  2.9736388620408842  -4.86016908678788 
H  1.726200370349726  1.4055422085473142  -4.579177161857819 
H  2.1078524082928523  0.7720261995955582  -2.7386561491182944 
H  3.473802795654417  -1.9331325378155748  -2.565189597082912 
H  4.908665032896289  -2.9605817601703524  -3.3461782561682334 
H  6.370402235537935  -3.7171069246203374  -2.2559857612539638 
H  6.348960990287363  -3.353382809886598  -0.3087920211851072 
H  4.96727132713506  -2.248396438276755  0.5136206387944569 
H  -0.06660121713627198  -4.721774937596119  -1.9882376819616534 
H  0.957601621991416  -5.380842169921667  1.152312970009711 
H  0.5814343555464243  -1.7447545148566932  -2.358844480509946 
H  0.6419450776865703  -2.485640682129537  -2.990357829260267 
H  -0.22525404445831065  -2.2480568064109194  -2.6140341327572427 
H  -0.1364003614079321  -7.165255365459914  -0.03468551267985942 
H  -0.8288013833773942  -6.855236995252555  -0.656012169849521 
H  0.04204339059859574  -7.144750885840534  -1.0038268716597523 
H  2.0170307086073658  -3.713195321125356  2.5613691818522 
H  2.5312573123471953  -3.0386787935941793  2.0894504768896827 
H  1.6347341525481687  -2.8263832458488016  2.425720568261889 
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

