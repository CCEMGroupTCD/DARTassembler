%chk=ATODILEC_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5629592621050616  -1.4315230857377055  -9.445241615798354e-16 
N  1.5629592621050614  1.4315230857377055  -4.528667726417115e-16 
C  1.677053052885628  2.7281487311728765  0.3523026575334755 
C  2.3889874321649067  3.64511068804104  -0.4672303919872703 
H  2.3722842331494474  4.5725985247571295  -0.26469383754248266 
C  3.0927585806753406  3.2012659562792947  -1.5379651091846205 
H  3.5550399370343158  3.820217158553055  -2.090385507764163 
C  3.14113771479229  1.8197178592181669  -1.8352801691611802 
C  3.924517140179881  1.2646083587727048  -2.8683579679714493 
H  4.465687170997994  1.831719398118354  -3.4038368627905125 
C  3.912387117467355  -0.07816647954534266  -3.1061541540824487 
H  4.412596587548696  -0.4349402613225377  -3.831741043901553 
C  3.1654866130531056  -0.9357463273529377  -2.2851980035237958 
H  3.1767177642778544  -1.8708849453306868  -2.448873494691008 
C  2.416584681270142  -0.43631217073815515  -1.2467631043710936 
C  2.351955869563983  0.9655994540162601  -1.0301643668958147 
C  1.1459943313260348  3.1704164293383843  1.6590655309787332 
C  1.3333865462643992  2.371807918312618  2.780635463661365 
H  1.7618534800962562  1.5285082907848546  2.696079891259128 
C  0.8956744073636591  2.8057093769676906  4.020153102759732 
H  1.020470929461098  2.253180540457591  4.783246100301445 
C  0.27965523485796684  4.034399348875193  4.156998351268278 
H  -0.030331751963897036  4.32218579752174  5.007647684411232 
C  0.11523669708005801  4.841488722082614  3.0500296217683722 
H  -0.30545030804806883  5.687913458136095  3.143009350274792 
C  0.5593417698222247  4.4245316286576015  1.8037267867369864 
H  0.4635329363481715  4.994517717690542  1.0494586542156081 
C  0.8344063303859941  -2.8100695406717797  -0.936909615316376 
C  -0.03559698339815687  -2.4840707463460068  -1.971826080594024 
H  -0.1877255307870409  -1.5720661388008832  -2.1914260143048767 
C  -0.6782191225094554  -3.475559905089525  -2.6800672714450666 
H  -1.263253569061475  -3.2478930457513457  -3.39328854672135 
C  -0.4694321063615501  -4.797148473132111  -2.3529651181241302 
H  -0.9093277498471029  -5.481941292589035  -2.844043551027329 
C  0.38050646391115306  -5.134607798666533  -1.308285339257569 
H  0.5109346547622673  -6.04700395770362  -1.0769618730772503 
C  1.0374062250761078  -4.142513770180693  -0.6059223829278758 
H  1.6261895150974823  -4.372347548816859  0.10228784122131633 
C  2.830641826944108  -2.1285666227518414  1.0960096907920023 
C  2.589582558704117  -2.108028909952854  2.4695444633946657 
H  1.7779325671452908  -1.747978502708188  2.80411120628093 
C  3.5391932721724344  -2.6138691420905302  3.343717719916002 
H  3.374130238794412  -2.607212738784892  4.278777155089562 
C  4.721398878087481  -3.123216236988585  2.8583859753721255 
H  5.368928056234037  -3.4702439971060524  3.4615382091876823 
C  4.972594642230851  -3.133993859196508  1.5000358029438228 
H  5.79794926613133  -3.4692874879071143  1.1718271233629796 
C  4.021398178047655  -2.657744634302822  0.6197553633601761 
H  4.18179229784033  -2.6919203527997313  -0.3156335218641299 
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
