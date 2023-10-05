%chk=LUXULUPA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.570237771007944  1.423535507987069  -9.914251097415421e-17 
N  1.570237771007944  -1.423535507987069  1.1102230246251565e-16 
N  3.689121796318001  -2.0484792343700993  0.18802863117197627 
C  2.793019004102819  -1.3155018171257375  -0.518634387242005 
C  3.0017837427464933  -2.642881920676458  1.251160935555744 
C  3.3922749382449133  -3.4138883787815897  2.3248695169680262 
H  4.295253776955787  -3.697272321603245  2.4202266236315126 
C  2.4509471896511132  -3.765022844480175  3.241232554593353 
H  2.709920693647489  -4.273593256625157  4.000555483338801 
C  1.1003544690954055  -3.3892414171827103  3.091958403849756 
H  0.4595849527654321  -3.680911538981858  3.7300546794947054 
C  0.6967918743137845  -2.604525395887548  2.0342721942318924 
H  -0.20500575425437662  -2.3219078512301006  1.9410962286579863 
C  1.6665416895877374  -2.2524464753492652  1.1127106827261342 
C  5.140555188824176  -2.218063771352503  -0.09689888932141846 
H  5.3213770000201865  -1.762627418439353  -0.91035519487259 
C  5.462747339491937  -3.694749033648819  -0.32652892691438157 
H  4.895563882793428  -4.038102955807844  -1.0045102692723886 
H  5.325310936814379  -4.177882282620256  0.47963572063433724 
H  6.369454781836717  -3.779895950827861  -0.5968271076177375 
C  6.005021666850213  -1.5835675224481132  0.9738801628156408 
H  5.7639763338846315  -0.6699309340365711  1.0748890325807023 
H  6.917871880035446  -1.6442679788225276  0.7197615037532407 
H  5.874165488845987  -2.0392364722990597  1.7971945309763502 
C  3.120944969707182  -0.3348104840559629  -1.5572464468312854 
C  2.7846799309548054  1.0071056109638148  -1.3130340441164723 
C  3.2779622818497023  1.9761475296656639  -2.2108705838575546 
H  3.0584882576342305  2.888401644658765  -2.0628449361716887 
C  4.058889993532788  1.6547074549835536  -3.2802240649449956 
H  4.447956441608775  2.345012501458675  -3.8013415294312054 
C  4.2965399934055695  0.30822819665414936  -3.6237661794654352 
C  5.003424226580315  -0.08613850993670691  -4.770308582354084 
H  5.3695106792901734  0.5848913587532114  -5.336878520197852 
C  5.181312833093482  -1.3798895155608903  -5.092929721027949 
H  5.67058986339458  -1.6169313266760517  -5.872112419374741 
C  4.6333549317136  -2.3857947191654136  -4.266536812251072 
H  4.733698880368242  -3.298495408898965  -4.5068923882801935 
C  3.950945089181827  -2.0500705329273132  -3.10660825945543 
H  3.585162723875283  -2.732728731158097  -2.5567777804720726 
C  3.801557923819642  -0.69474965109826  -2.746258773396154 
C  2.42817479902451  1.1571817409945258  1.5624266852454367 
C  1.6693851799751098  0.770389629445031  2.665387535169832 
H  0.7264000403318092  0.67833057699097  2.597207499432144 
C  2.342658047063224  0.513215464876648  3.8992625310756392 
H  1.8378594762953155  0.23535552546473126  4.654079288322517 
C  3.6708204972531613  0.6503212018989637  4.011515071274344 
H  4.102081937278806  0.48962966753501647  4.843637268509101 
C  4.393256774657152  1.0266467146763119  2.9261685581627193 
H  5.336025388880269  1.1064137976629294  3.0010172221966975 
C  3.7807035654774044  1.29241495359785  1.7015106675237837 
H  4.306412835248758  1.5729021193620567  0.9607569948916167 
C  1.2834164950421125  3.212513675873509  -0.1003668295753386 
C  0.6814315817150168  3.733009069801351  -1.2160867671163476 
H  0.4477499589606855  3.1584679744713777  -1.9364987318725815 
C  0.4068894232622171  5.1051064983803816  -1.304345646096825 
H  0.016638757625042633  5.471577715341931  -2.0892096311977246 
C  0.7174820960337509  5.9233605977439385  -0.21799995503186595 
H  0.5226756624770172  6.853217592931742  -0.25195097777368824 
C  1.2976391159343137  5.392800038655402  0.886097352907873 
H  1.4994571944026691  5.9644100901325245  1.6194388778835718 
C  1.6035251805870845  4.0601893595705825  0.9742777557213049 
H  2.0247870091134885  3.7148863184858762  1.7531964816768775 
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


