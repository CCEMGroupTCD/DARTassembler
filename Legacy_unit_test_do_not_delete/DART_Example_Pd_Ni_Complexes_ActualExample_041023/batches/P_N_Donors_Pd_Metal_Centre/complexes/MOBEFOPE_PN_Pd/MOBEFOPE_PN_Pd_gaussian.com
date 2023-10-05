%chk=MOBEFOPE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5251001229427508  1.4717912946474443  -9.362913210253328e-16 
N  1.3546150911655663  -1.0286024251574615  -2.310744752217388 
N  1.5251001229427528  -1.4717912946474443  4.440892098500626e-16 
C  1.6933312980494097  1.4700145243882976  -1.8217251507320635 
C  1.8895800895828119  2.6806819821183203  -2.50091116988676 
H  2.1108853526130726  3.43573469852007  -2.005481017383456 
C  1.7710868128843744  2.805519617500592  -3.877114468156841 
H  1.9196079804487667  3.6298410035641093  -4.283323571038133 
C  1.4316950443285232  1.7112300630934913  -4.64615844515172 
H  1.2857555874485402  1.8019466344481159  -5.560151671935838 
C  1.3132453102606123  0.4762956185088174  -4.031333439090947 
C  1.0237154646669475  -0.833368850910218  -4.568748824106937 
C  0.6964537034913756  -1.295167746903016  -5.848916504269255 
H  0.6337693088852878  -0.704533671839926  -6.564674910470169 
C  0.46938586480542344  -2.6405889130706948  -6.026351033779449 
H  0.2653270232028111  -2.9661361029083486  -6.873778591953327 
C  0.5414165860967433  -3.51761171154376  -4.9516708878063405 
H  0.38613599800478116  -4.422630096177615  -5.097079011501462 
C  0.8382203334580605  -3.0893424913816783  -3.675979839544994 
H  0.8733297939817357  -3.6831552746060514  -2.9619328143001096 
C  1.0803823045526253  -1.735801787448382  -3.5037395639084665 
C  1.4830762201420942  0.34086136521154775  -2.639288136278339 
C  1.9777916835038636  -1.7003479461396513  -1.234929628786785 
C  3.007060554707963  -2.602192891701685  -1.465924011693486 
H  3.3223256072930747  -2.74537361016794  -2.328808058641861 
C  3.551753279033732  -3.2776870827704996  -0.4104647211953305 
H  4.246407634456255  -3.8793000921830214  -0.5514887737470773 
C  3.0701299314706336  -3.0659634497165045  0.8651451385535502 
H  3.4245668196425463  -3.52846764271967  1.589321728728548 
C  2.061916073803583  -2.1618723602806607  1.0348554187691743 
H  1.7315376078943576  -2.015582495146232  1.8916848852457002 
C  1.1358026750589492  3.2274989003368675  0.4311835808509892 
H  1.86780632909232  3.7775672306390446  0.08030366908574488 
C  -0.16084601582817615  3.72795727556896  -0.23490945850673917 
H  -0.10450212936544734  3.6004423923905287  -1.1943138211329853 
H  -0.9150174907240765  3.2153544848446667  0.09361668260690181 
C  -0.37465228308206533  5.206349329464901  0.07284600458093349 
H  0.34216424429767334  5.7253071448683155  -0.32353794048620865 
H  -1.2093840022075295  5.499383098293729  -0.3237908159637597 
C  -0.4120607990375935  5.460999215194644  1.5674341514066252 
H  -1.185954562781202  5.019523666070917  1.9504473260686508 
H  -0.4993556881658592  6.4133140271501965  1.730226630539248 
C  0.8504213729791149  4.945303880020228  2.2449877335086836 
H  0.7734232236780556  5.071786368180748  3.2036334215462077 
H  1.611691541637252  5.460042398147879  1.935332395645807 
C  1.0830050748561202  3.4652609806485493  1.947078724513213 
H  0.36518566400360375  2.9393044215798745  2.3338396860592185 
H  1.9171497255222876  3.181049746677228  2.3507239849208204 
C  3.2408137801894785  1.116943095789805  0.5997600120920326 
H  3.5099028910707535  0.2915243530522311  0.14220016951114145 
C  3.3388371064807916  0.8012951967367785  2.0929084445414734 
H  3.161341736762482  1.6054025245243726  2.6077795440672045 
H  2.672050024259305  0.13899348755762597  2.330180973505394 
C  4.726463659073931  0.2756199479337009  2.4298921434563794 
H  4.7896909899261075  0.11760791467479526  3.384842970451524 
H  4.868190844347657  -0.5715260533603661  1.977206328076813 
C  5.803523949240331  1.2526059718552807  2.011423099907653 
H  6.674991024561501  0.8528093993793231  2.164766789923376 
H  5.741191291551423  2.050412562552144  2.5595679198467622 
C  5.687239541841131  1.6431059952554539  0.5509529843899466 
H  5.892396157957841  0.8734114916122864  -0.0011017537701887292 
H  6.337157060522886  2.334906516635694  0.3520816000014393 
C  4.282296241810709  2.1544443030284155  0.2030898378745003 
H  4.112027398112871  2.9856414024059084  0.67369174272853 
H  4.225124116725819  2.3288431739372726  -0.7479640627741433 
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


