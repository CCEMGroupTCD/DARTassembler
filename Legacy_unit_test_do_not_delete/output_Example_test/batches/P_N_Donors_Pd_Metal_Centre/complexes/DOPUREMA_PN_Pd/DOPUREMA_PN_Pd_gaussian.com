%chk=DOPUREMA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5197065555560383  1.4773598021470606  -1.1382186207774173e-15 
N  1.5197065555560385  -1.4773598021470606  -2.498001805406602e-16 
N  2.018498152210419  -4.313326595010254  -0.4615477796330938 
C  1.4978394778987514  -2.1078002091802053  -1.3543123657763187 
H  1.7849139906731766  -1.436987700568745  -2.0239315205245463 
H  0.5678648487367878  -2.3709827591825348  -1.5702601013267363 
C  2.3931591411536637  -3.3287535256157414  -1.469051677584443 
H  2.3034114232889262  -3.7241270014278696  -2.3722865704901097 
H  3.339196502696632  -3.0626155378157587  -1.341599953333288 
C  2.1953519331251856  -3.74410967928136  0.8659409785205409 
H  3.136736995179231  -3.460625774824678  0.9871265009871856 
H  1.984211304066786  -4.424221557606302  1.5539790644737501 
C  1.2731272323026535  -2.5438067859557942  1.0245106334583232 
H  0.333487595724717  -2.8492187575988988  0.9515404960816617 
H  1.396652350997119  -2.1587559512347787  1.9285580988782465 
C  2.8519370912532587  -0.8225930089838941  0.2733513605723642 
H  3.5805914005998414  -1.4307310662765662  -0.008912152377790243 
H  2.9475270880609825  -0.6617823251851812  1.2451530412081404 
C  2.980120957064452  0.49468324773551187  -0.47243783218187163 
H  2.995991764461838  0.3427985919638481  -1.4510648442828549 
H  3.812423060122561  0.9623267723501261  -0.2095757331050493 
C  2.7550087019982583  -5.560492391942229  -0.6192161182881946 
H  2.608412176248282  -5.915229358684087  -1.5217113000287095 
H  2.4404894540015905  -6.211745487905528  0.04317490492832923 
H  3.7114595552780134  -5.394765830515434  -0.48485030381854943 
C  1.3026374663594682  2.768669334345975  -1.2470955025849593 
C  1.133636057845619  4.110211853046991  -0.9127127362262938 
H  1.2052628148030355  4.389000902540163  -0.006586076442753692 
C  0.8562001405006674  5.038262665724727  -1.9072582660313238 
H  0.7494055935853444  5.955297296036904  -1.6821687248319208 
C  0.7359065508614413  4.632971230199127  -3.226135305251179 
H  0.5396673982954926  5.272378061647468  -3.900578808368674 
C  0.9003157821161833  3.3029048393684937  -3.5679933908086223 
H  0.8111834777371677  3.028120948260841  -4.473559233530107 
C  1.1933822055824574  2.366268751809465  -2.5780020166447626 
H  1.319016659255716  1.453822164295088  -2.8103184366059644 
C  1.8901155056068824  2.2483412268106098  1.5904093914920325 
C  2.9921554348338564  3.0886124614792996  1.7483447788432145 
H  3.5761154366827625  3.2479386693928465  1.0155439766380525 
C  3.2379232045644337  3.697994910962234  2.9750658733157382 
H  3.9804686687948836  4.2823970239378655  3.0782043428623713 
C  2.3981795317957824  3.4479694165206096  4.045573638899111 
H  2.5615598952384113  3.8699834574427383  4.88105040464783 
C  1.325916492457267  2.588293071493335  3.9085743094514687 
H  0.7678995552409598  2.402370996450536  4.65572597396915 
C  1.0618406403773821  2.0003142638003313  2.685744003950231 
H  0.3129443655778721  1.4241079674310577  2.5886681901971533 
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
