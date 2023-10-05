%chk=UYEFATIT_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.4485495020882095  -1.5471924056173492  2.0250751904671284e-16 
N  1.4485495020882095  1.5471924056173492  3.2568182204593296e-17 
C  3.2000918211492606  -1.065507718377282  -0.0393722823410207 
C  4.161997876510274  -1.9903645395488152  0.36292708503647325 
C  5.5002917442395045  -1.6593896988368269  0.39029601658414215 
C  5.8751256485951195  -0.39633760337346025  0.014603538230972817 
C  4.931061821847928  0.5508757649698512  -0.3738268105566052 
C  3.576563211766826  0.23512780225601276  -0.4083405662265653 
C  2.644792735163494  1.3100478729348748  -0.8794106253283808 
C  2.2208181395317377  1.1016253400752132  -2.3225105837703826 
C  1.933998606867823  1.8068654356434317  1.3840679225189052 
C  0.8435077939836679  2.835672825222665  -0.44975646715491013 
C  1.1591186618011955  -2.8544681346453453  -1.2651376276852688 
C  1.1302664514079672  -2.3123878885460596  -2.693354392306254 
C  0.6816814527968034  -3.3837760584864234  -3.678571453837803 
C  1.5133218904890504  -4.635195302124344  -3.5873721717805673 
C  1.5701713977371008  -5.139285050217957  -2.1807517825523086 
C  2.076629009409563  -4.072529723596871  -1.191295921634505 
C  1.2765719689814525  -2.432804644292517  1.6151499189905962 
C  1.6026162987420827  -1.5297640496468181  2.804429662975563 
C  1.5517006928178025  -2.34334518243278  4.095433257644491 
C  0.21817357543341998  -3.045778716913808  4.2942854058275834 
C  -0.08822677792017863  -3.919861394084388  3.0940401155076684 
C  -0.08424612797709408  -3.098374819185608  1.8073451043697717 
H  3.9049721987547237  -2.816020518315259  0.6087402777286754 
H  6.114574798710979  -2.261595673101694  0.6485787174158574 
H  6.746428199154426  -0.17366728858314517  0.01973737695720492 
H  5.199241472607939  1.3772993341120532  -0.6067928365872256 
H  3.1599224209514345  2.1040890054660277  -0.8210056026962338 
H  1.6317318908616534  1.7996615029201743  -2.583520175526118 
H  2.988352056317576  1.1078906629689944  -2.8809817220689 
H  1.7762378702874735  0.26533321163517204  -2.405163352072039 
H  1.1921374774445872  1.9595345213643838  1.9568324770402563 
H  2.4228479631042665  1.0541271299410475  1.693599374936139 
H  2.4967460364764396  2.5719756609147737  1.383490294493819 
H  0.07841103532814575  3.0278031762512354  0.0810727314659535 
H  1.4786258103438206  3.533827072928525  -0.35520610118511786 
H  0.5834819539742926  2.76122582717523  -1.3614029265947314 
H  0.29236763531076404  -3.166221086798404  -1.0315681414469031 
H  1.9988490948717677  -2.01720875740382  -2.931013332831452 
H  0.5232288175724313  -1.5819331655931934  -2.7353199594204254 
H  0.752865626683227  -3.0354981015603157  -4.558831783574163 
H  -0.22298589790218326  -3.604990070885504  -3.4966016066169474 
H  1.1299214569892277  -5.302627688846772  -4.141504255986275 
H  2.393290754699763  -4.444426591698683  -3.882497259060981 
H  0.6959958268649052  -5.4077192788840955  -1.9181064586443808 
H  2.152428184555887  -5.887090200889144  -2.1481826125734003 
H  2.0666076137788334  -4.425161890106441  -0.3108001881000369 
H  2.9605631195810105  -3.8192985814470566  -1.423426870714794 
H  1.9180963753588163  -3.1303037633949016  1.5789892376699628 
H  0.965204703663523  -0.8254029203888428  2.8494163059913733 
H  2.470899773706857  -1.1650941495416207  2.695593404952864 
H  1.6985242662556854  -1.757316116682412  4.829727736977364 
H  2.2377145211227503  -2.9973275289892594  4.0699117418803885 
H  -0.4682653541151094  -2.395904404999254  4.393948160725635 
H  0.2596830496648419  -3.5851182890124957  5.0741271074039105 
H  0.567837841374371  -4.6007370926563835  3.02972372254279 
H  -0.9462165513948197  -4.3143191407642405  3.205833383181511 
H  -0.26281012924056224  -3.669664994079735  1.0685008626501566 
H  -0.7575107827904093  -2.4295978162458005  1.8579476232593457 
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

