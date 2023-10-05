%chk=IXIYETAX_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.3681252811785884  1.6187443328086126  1.6289042896069139e-15 
N  1.368125281178591  -1.6187443328086126  -2.220446049250313e-16 
C  2.2116202791244657  -1.7555464797021725  1.0537071704442607 
C  3.0796093253951025  -2.84101579340864  1.1484952211459603 
H  3.657269216876159  -2.9265081050696726  1.8984490426716478 
C  3.0976621374781725  -3.791085380575393  0.14589631554276872 
H  3.6836222710997495  -4.537089156789248  0.20001089473247702 
C  2.2481960637147878  -3.6420463877522278  -0.9400342292635399 
H  2.2466002365978817  -4.27752694150662  -1.6475289989515125 
C  1.4078741851138479  -2.551741029315697  -0.9689184067519966 
H  0.8250028098559824  -2.452840686157597  -1.712122774357735 
C  2.22895068536903  -0.6759948917327748  2.0987080728037695 
H  2.46168403010449  -1.0708802556808426  2.9764370740267525 
H  1.3264236246269945  -0.27513068986457245  2.1738274863327227 
C  3.251693354849646  0.42641786044270225  1.7439255818820854 
H  3.1802622727127696  1.1513341064591391  2.414452156125205 
H  4.162654282206511  0.0454660927660524  1.8108867611816206 
C  3.0777776258325886  1.0347665323366666  0.3432067583470507 
H  3.325933325853505  0.35623575836677834  -0.3333751785408452 
H  3.7013740050129718  1.7974090942154985  0.24768598215652415 
C  1.062884380368339  2.783596639490287  1.3625629108745676 
C  0.06191976867448501  2.5385512781194994  2.3022352404699946 
H  -0.4950484625818381  1.7751554193949066  2.210071460961153 
C  -0.12461889484106403  3.4069837269676193  3.373034039872344 
H  -0.8189173167325297  3.239803988008808  3.9990189408409877 
C  0.6882640178723483  4.503580712064139  3.5292213712913965 
H  0.5530204355407309  5.090187327917086  4.264046787070991 
C  1.7034842652931625  4.757634268678888  2.6222244034914706 
H  2.2733258581490294  5.507891684495021  2.7394035358674715 
C  1.8823469116854712  3.900884606677774  1.5347462380156813 
H  2.5689906989825806  4.081521581939123  0.903265032804684 
C  1.5689526132375908  2.6572125862129945  -1.4891273852056346 
C  0.8721330829301179  3.8529136742543284  -1.6503391721160605 
H  0.3442862497886958  4.191553025070705  -0.9368862368629698 
C  0.9453578857908695  4.552744115071656  -2.850310135244552 
H  0.4596991777420635  5.362543589535658  -2.9573142494373945 
C  1.7254411424700773  4.071974802570214  -3.8899466394159266 
H  1.7768652765562198  4.553531002768608  -4.706534873769064 
C  2.4320036222595958  2.8880013138879854  -3.7373620391073166 
H  2.968259625939891  2.561205689478195  -4.4496351457058125 
C  2.3566679626806204  2.1796734255607517  -2.5446785769038307 
H  2.8401925416295155  1.3684869223025804  -2.444206050887602 
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


