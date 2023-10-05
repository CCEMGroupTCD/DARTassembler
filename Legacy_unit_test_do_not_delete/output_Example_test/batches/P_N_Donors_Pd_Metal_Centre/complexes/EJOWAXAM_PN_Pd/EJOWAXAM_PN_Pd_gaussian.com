%chk=EJOWAXAM_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.570237771007944  -1.423535507987069  2.734753313070508e-16 
N  1.570237771007944  1.423535507987069  -2.853551227954122e-16 
N  3.689121796318001  2.0484792343700993  -0.18802863117197652 
C  2.793019004102819  1.3155018171257375  0.5186343872420048 
C  3.0017837427464933  2.642881920676458  -1.2511609355557443 
C  3.3922749382449133  3.4138883787815892  -2.3248695169680267 
H  4.295253776955787  3.6972723216032444  -2.420226623631513 
C  2.4509471896511132  3.7650228444801748  -3.2412325545933536 
H  2.709920693647489  4.273593256625156  -4.000555483338802 
C  1.1003544690954055  3.38924141718271  -3.0919584038497563 
H  0.4595849527654321  3.6809115389818574  -3.730054679494706 
C  0.6967918743137845  2.6045253958875474  -2.034272194231893 
H  -0.20500575425437662  2.3219078512301  -1.9410962286579865 
C  1.6665416895877374  2.2524464753492652  -1.1127106827261344 
C  5.140555188824176  2.218063771352503  0.09689888932141819 
H  5.3213770000201865  1.7626274184393533  0.9103551948725898 
C  5.462747339491937  3.694749033648819  0.3265289269143811 
H  4.895563882793428  4.038102955807844  1.0045102692723882 
H  5.325310936814379  4.177882282620256  -0.47963572063433774 
H  6.369454781836717  3.779895950827861  0.596827107617737 
C  6.005021666850213  1.583567522448113  -0.973880162815641 
H  5.7639763338846315  0.669930934036571  -1.0748890325807023 
H  6.917871880035446  1.6442679788225276  -0.7197615037532409 
H  5.874165488845987  2.0392364722990597  -1.7971945309763504 
C  3.120944969707182  0.3348104840559631  1.5572464468312854 
C  2.7846799309548054  -1.0071056109638146  1.3130340441164725 
C  3.2779622818497023  -1.9761475296656636  2.210870583857555 
H  3.0584882576342305  -2.8884016446587646  2.062844936171689 
C  4.058889993532788  -1.6547074549835532  3.2802240649449956 
H  4.447956441608775  -2.3450125014586747  3.801341529431206 
C  4.2965399934055695  -0.3082281966541489  3.6237661794654352 
C  5.003424226580315  0.08613850993670749  4.770308582354084 
H  5.3695106792901734  -0.5848913587532107  5.336878520197852 
C  5.181312833093482  1.379889515560891  5.092929721027949 
H  5.67058986339458  1.6169313266760523  5.872112419374741 
C  4.6333549317136  2.385794719165414  4.266536812251072 
H  4.733698880368242  3.2984954088989653  4.5068923882801935 
C  3.950945089181827  2.0500705329273137  3.1066082594554296 
H  3.585162723875283  2.7327287311580974  2.556777780472072 
C  3.801557923819642  0.6947496510982604  2.746258773396154 
C  2.42817479902451  -1.157181740994526  -1.5624266852454365 
C  1.6693851799751098  -0.7703896294450313  -2.665387535169832 
H  0.7264000403318092  -0.6783305769909703  -2.597207499432144 
C  2.342658047063224  -0.5132154648766485  -3.8992625310756392 
H  1.8378594762953155  -0.23535552546473185  -4.654079288322517 
C  3.6708204972531613  -0.6503212018989641  -4.011515071274344 
H  4.102081937278806  -0.4896296675350171  -4.843637268509101 
C  4.393256774657152  -1.0266467146763123  -2.9261685581627193 
H  5.336025388880269  -1.1064137976629298  -3.0010172221966975 
C  3.7807035654774044  -1.2924149535978502  -1.7015106675237834 
H  4.306412835248758  -1.572902119362057  -0.9607569948916165 
C  1.2834164950421125  -3.212513675873509  0.10036682957533899 
C  0.6814315817150168  -3.733009069801351  1.216086767116348 
H  0.4477499589606855  -3.1584679744713773  1.936498731872582 
C  0.4068894232622171  -5.1051064983803816  1.3043456460968257 
H  0.016638757625042633  -5.471577715341931  2.0892096311977255 
C  0.7174820960337509  -5.9233605977439385  0.21799995503186667 
H  0.5226756624770172  -6.853217592931742  0.2519509777736891 
C  1.2976391159343137  -5.392800038655402  -0.8860973529078723 
H  1.4994571944026691  -5.9644100901325245  -1.619438877883571 
C  1.6035251805870845  -4.0601893595705825  -0.9742777557213045 
H  2.0247870091134885  -3.7148863184858762  -1.753196481676877 
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
