%chk=MUBUHAWE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.519864348387711  -1.4771974690270755  4.497347850213488e-16 
N  1.5198643483877106  1.477197469027076  -6.249937250653206e-16 
N  1.8253427767790174  2.5732620839132547  -0.8004056720076441 
C  2.6466694868380314  3.682805814095449  -0.27518578838075863 
H  2.6066818403790304  3.668221311717314  0.7139950699161826 
H  2.2560676341479313  4.540180357144295  -0.5790003792200601 
C  4.085392935541245  3.628194374737153  -0.702090507074201 
C  5.109943142079663  3.2638528013966055  0.19900022210785356 
H  4.877770228948285  3.037145711097725  1.0917771335567907 
C  6.425187725903079  3.2282390320102765  -0.17258950346783702 
H  7.094478741726846  2.9880678785668113  0.45811521602554933 
C  6.783935400011999  3.551487649755673  -1.5125526997710554 
C  8.159038876528323  3.5144864389182353  -1.944864727442377 
H  8.839041546546866  3.2741087998786353  -1.3264639840587664 
C  8.494305813495037  3.8191578958064043  -3.23161923951318 
H  9.406664395898906  3.8176962327279944  -3.4978566348112263 
C  7.502172429040563  4.133663132744014  -4.15037431753439 
H  7.743701952122281  4.3154439580096176  -5.0510274311390235 
C  6.160499087160019  4.18920811350684  -3.779557612410351 
H  5.500806133247427  4.418672384267641  -4.424060153136371 
C  5.777294776620703  3.906737727015573  -2.446892963243263 
C  4.420768449639799  3.9489195620826476  -2.0167605336380374 
C  3.3011492574861028  4.306679405354023  -2.9274714916832303 
C  3.1958417788384885  5.567170963890446  -3.615047651682475 
C  4.121639854278286  6.635725329908704  -3.378133362591863 
H  4.875717479642163  6.479539853804443  -2.821693822318449 
C  3.9512615231285046  7.85591723234078  -3.923228359371157 
H  4.566182059745955  8.550663426751697  -3.720223164862087 
C  2.8604313635877485  8.11616358264979  -4.797576007728007 
H  2.7624424775668803  8.971968934680778  -5.197795781916519 
C  1.9632678211275871  7.13637263029102  -5.059344469839533 
H  1.2353798948174857  7.31909725625303  -5.641674417377069 
C  2.0823700921494845  5.827888865914512  -4.477706036274058 
C  1.0875160744350687  4.825454116296851  -4.66527254766796 
H  0.380385462893605  4.959581390513274  -5.286390783217606 
C  1.1625048699928855  3.673744821722631  -3.944375875923572 
H  0.4745785447948794  3.025699103345859  -4.036957882777532 
C  2.2432984637867803  3.4142247419512426  -3.054350221344113 
C  2.175466457273516  2.1753573099583328  -2.1995513181863737 
H  1.4885971655120165  1.5568616997527802  -2.556106455846068 
H  3.048905069961619  1.7103330431266055  -2.209790025287038 
C  1.7815285020259322  1.4591187181762584  1.2511182396034697 
H  2.1413249700437555  2.2553834895847173  1.6267844550337969 
C  1.5782078763223404  0.3172474075579477  2.170308885307793 
C  1.541779951048537  0.637607156286522  3.5204965676893223 
H  1.577764469280068  1.5479654483879215  3.7887267569855285 
C  1.454592163826088  -0.3668107100514296  4.491881125699404 
H  1.4016616850917856  -0.14103386437373944  5.413205846075872 
C  1.4429114492682726  -1.692252846971387  4.09227899534411 
H  1.4146847983686321  -2.382407339770151  4.744655330176009 
C  1.474069751060853  -2.0218989607849456  2.7341110565835716 
H  1.4598350746468947  -2.934876273079053  2.475011864721407 
C  1.527789969687099  -1.0363444391143473  1.7551031082371946 
C  1.0078083992487397  -3.20758321541177  -0.05367423646474589 
C  1.8483278760141877  -4.252846693974972  -0.39224642900262424 
H  2.74374321772806  -4.075386875300456  -0.6563743309354764 
C  1.377560390830217  -5.584628335393383  -0.34668674864321636 
H  1.9562273039992073  -6.3006522199297095  -0.5844181524942597 
C  0.08129653559124694  -5.851746732407619  0.03921581623239228 
H  -0.233369084560058  -6.747648666812936  0.053511319847402365 
C  -0.7603274593600833  -4.822453119413938  0.4037411094187198 
H  -1.6467827561933766  -5.015199228833911  0.6868631396806988 
C  -0.3141494919200791  -3.4854971325472963  0.35919225362487495 
H  -0.8981873283675006  -2.776650030548365  0.6049633698804674 
C  3.2149775076260223  -1.3762981491154547  -0.6069825076709756 
C  3.4975192951291145  -1.7819034609053692  -1.9145206789540814 
H  2.8106167793638424  -2.17525941850981  -2.43953956856896 
C  4.751197599425907  -1.6229573450860397  -2.4549875619372092 
H  4.93598908032506  -1.9494008268214935  -3.3288239206563683 
C  5.736180130678617  -0.9941200640339605  -1.7379253421221381 
H  6.586830142016495  -0.8410727059519671  -2.131115581567325 
C  5.487797670519132  -0.5827610698974335  -0.4312677927129011 
H  6.175745239670061  -0.1676963728235712  0.07666488918729454 
C  4.2371882785366415  -0.7778322476119393  0.1243813720577812 
H  4.072704487580202  -0.5009599820038927  1.0178402245236828 
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

