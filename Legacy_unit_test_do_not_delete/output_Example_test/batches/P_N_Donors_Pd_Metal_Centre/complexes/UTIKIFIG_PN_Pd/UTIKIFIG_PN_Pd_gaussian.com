%chk=UTIKIFIG_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5055490086676027  1.49178489820081  -1.0800584078999045e-15 
N  1.5055490086676022  -1.49178489820081  -1.1102230246251565e-16 
N  1.428483993311091  -2.6182023475439697  0.7502477478962306 
C  2.3458348703894676  -1.4478675780301766  -0.9815880600228494 
H  2.830933489580004  -2.24486063073767  -1.1605459701453942 
C  2.6303196925998575  -0.2967150411790913  -1.8513981731126221 
C  2.4127114842894417  1.0457595289920965  -1.5134970838577635 
C  2.8525742013378492  2.0568107849868937  -2.379564083518981 
H  2.725474223869032  2.9677449529349307  -2.141251605053744 
C  3.4698081814498654  1.7425472436644018  -3.577516323965379 
H  3.773281979248468  2.435239994640958  -4.151997495801519 
C  3.643105546503322  0.4191107030605159  -3.93699894699458 
H  4.025192098043864  0.2002222759194301  -4.779066143597987 
C  3.258284345370593  -0.5855304464829209  -3.0694547672533847 
H  3.4230605647497834  -1.490868978619073  -3.304860416626257 
C  2.313735752362092  -3.786376872639221  0.514165681623493 
H  3.2652659944148232  -3.486856596264997  0.447942328493366 
C  2.1135102572241773  -4.612499198643696  1.7956389796573462 
H  2.0786862164683537  -5.580541593188441  1.5889091352963716 
H  2.8474405662067754  -4.447241264810653  2.4386629969863947 
C  0.7788420173698813  -4.127397991443656  2.357204085675676 
H  0.7044942098354245  -4.320496326429767  3.3258313599873444 
H  0.016106418802309008  -4.542725375939407  1.8827976983180577 
C  0.841571265897583  -2.60198826733095  2.0914059950851813 
H  -0.0799632898969096  -2.212859137052102  2.066276879599137 
C  1.9249849393498812  -4.561726550327339  -0.7403407881451227 
C  0.7117743134797058  -4.393538026430495  -1.3669633521381157 
H  0.09944001405426839  -3.740563877024085  -1.0499307261291513 
C  0.38353277888508375  -5.186114812417534  -2.4707060863300496 
H  -0.4562047312790338  -5.071534030559403  -2.902007243819991 
C  1.266946148224709  -6.128814099490419  -2.9357192936955356 
H  1.0405785576335793  -6.662891794579643  -3.6883626551477953 
C  2.4810830559166077  -6.299985741131229  -2.3130379237490812 
H  3.0919855445201763  -6.95188564384032  -2.639288787756903 
C  2.82107580663973  -5.526426408688742  -1.2100268863396924 
H  3.657782428317632  -5.653256939168164  -0.777123426225675 
C  1.6850668210218362  -1.9049026377243008  3.1506220700583714 
C  3.044861653733042  -1.6762804647747456  2.9878379953833667 
H  3.4618493236851737  -1.8736595945374364  2.15705760789096 
C  3.805450029566403  -1.1553416444815063  4.039036515197156 
H  4.7348321221732315  -0.999890429443945  3.921151038207335 
C  3.2135287830019634  -0.8705154336032108  5.2375880023942365 
H  3.7353152202356936  -0.5319883908579235  5.954473603622943 
C  1.848865149318789  -1.0751244026094873  5.411080989874627 
H  1.4382529406612339  -0.8664750186564789  6.240819809430178 
C  1.0907383479752109  -1.5836768493563227  4.363858448104056 
H  0.15696484159490542  -1.71389928305591  4.479869337867093 
C  0.9287195753306873  3.1921314582573226  -0.32625407359941544 
C  1.3961789911643454  4.251442986283852  0.42286614448214865 
H  2.009374098820262  4.101767949172741  1.1331279857862238 
C  0.9556603530670527  5.559578054787083  0.12575178595842915 
H  1.2765906655287764  6.296669086695347  0.6322689872691429 
C  0.06616811572000891  5.764546720422977  -0.8951459688666172 
H  -0.23112363728454666  6.645382839346552  -1.0913742675269775 
C  -0.4060595302043215  4.6917433075916914  -1.6502539266066747 
H  -1.0213265878067872  4.8407128810657385  -2.358413544212903 
C  0.024157671324325003  3.4029340581356786  -1.3610694560251702 
H  -0.2965647230502737  2.667241461754223  -1.8694335867190621 
C  2.6699558607459615  1.63910815188508  1.3721441330111483 
C  4.044515864620886  1.6783508906141054  1.1628152453618044 
H  4.394511336402229  1.5451185692253377  0.29009296521195316 
C  4.906683294193456  1.9112665880575241  2.2275217829902085 
H  5.846139137388699  1.9402723281291563  2.08895802247096 
C  4.378649993826633  2.1023740160411686  3.495027603360779 
H  4.960987677685001  2.2711675500846438  4.2262496112513865 
C  3.01885310472284  2.051582759143887  3.7063217315481944 
H  2.670701960759461  2.184569485473715  4.57993739737447 
C  2.1562426955481677  1.8050866420136675  2.649978564663845 
H  1.220630311214102  1.748827787200739  2.8005283821193867 
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
