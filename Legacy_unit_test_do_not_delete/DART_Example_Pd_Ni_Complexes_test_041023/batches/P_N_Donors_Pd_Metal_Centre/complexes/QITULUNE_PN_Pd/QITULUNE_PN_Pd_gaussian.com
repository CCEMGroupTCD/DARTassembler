%chk=QITULUNE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5572926715938815  -1.4376854784687787  -3.949119182438976e-15 
N  1.5572926715938777  1.4376854784687787  -6.504338947622643e-17 
C  3.0727663985179694  -0.4442337942233517  -0.25404840891605934 
H  3.814600381960844  -0.8191729502994869  0.2848224942232793 
H  3.3368115344220635  -0.47047133057252316  -1.2073920237222473 
C  2.8127791571516756  0.9702231996428232  0.16061684985276747 
C  3.8057555434251147  1.7746093343846805  0.6849216300507135 
H  4.685123497560955  1.4328017086142777  0.7987084843184941 
C  3.5226373983914243  3.062667190385496  1.0398077181793837 
H  4.2016064963326265  3.6264908827993194  1.3887791466857475 
C  2.2411356915291116  3.5302750712941986  0.881873602102155 
H  2.0249796559287687  4.424701289168753  1.1192891008636354 
C  1.2740183987824072  2.693912951061437  0.37625715071822907 
H  0.38338841279235836  3.0138443803542216  0.29169457858341585 
C  1.7733480576465328  -2.9695359095836578  -0.9557715819435157 
C  2.0668617523550994  -2.950534326546807  -2.310249234347085 
C  2.3444899795237943  -4.207361544363933  -2.921671651365981 
H  2.5890522254017716  -4.2365408507061515  -3.8400782802469617 
C  2.2664148062776794  -5.34809914278965  -2.2244646169346534 
H  2.4661367657298383  -6.16872914682681  -2.6592649681150133 
C  1.9107974028348473  -5.368685432658503  -0.9126043091891766 
H  1.848892797489627  -6.194789711227794  -0.44562897711416405 
C  1.6462571751418187  -4.206096684958149  -0.28013301999680706 
H  1.3710865296593113  -4.218733496188483  0.6285351190213346 
C  2.0229812318656806  -1.7586505552611442  -3.146788525514847 
H  1.8294152377919097  -0.976004304150362  -2.589192048158251 
H  1.3212372206443248  -1.859844079282093  -3.822763823444586 
H  2.8882981982982807  -1.637884137570871  -3.590728014592086 
C  1.5410708637200607  -1.8577175293198178  1.767504350548324 
C  2.617740047274657  -2.5112782761395933  2.4406449000625337 
C  2.4617069794117628  -2.7128998577310868  3.783562994947641 
H  3.160041799296919  -3.1498962240497748  4.257886044314079 
C  1.3450514779156975  -2.3154130300776283  4.48420467428208 
H  1.279402035865282  -2.4860908095188914  5.416274900372086 
C  0.3165564025225991  -1.663045308454038  3.8216741197596162 
H  -0.44966346982084504  -1.3665586941119445  4.297368501418683 
C  0.41803688991931764  -1.454063372210236  2.4624707203371488 
H  -0.29312773816389215  -1.0259391306998107  1.999360673378104 
C  3.8547741706425898  -2.9995111584439984  1.7884759033598971 
H  4.441209405862232  -3.4065734879761167  2.460265092591754 
H  4.317618199229946  -2.247648288443205  1.3615297124702275 
H  3.6240665959014855  -3.6686808999542526  1.1081414914680878 
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


