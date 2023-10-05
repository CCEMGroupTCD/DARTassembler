%chk=TESOJAGA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.4419963540175806  -1.5533018106601177  -2.2286051768574793e-15 
N  1.441996354017582  1.5533018106601182  -1.9022460905347018e-16 
C  1.1253104140887498  2.673809453189894  0.7196396031929929 
H  0.2987453002354836  2.6975172733619104  1.1885126263276735 
C  0.6205352857928879  -3.4057985245969493  1.8928940113605104 
H  0.23146479351851457  -3.873461981894593  1.1628554557007367 
C  1.991299497798466  -1.6095952363162938  2.7125710527948415 
H  2.525068821981853  -0.839065126336341  2.555091536110607 
C  1.9639102956275427  3.771992641603914  0.7904504312346825 
H  1.7186838020787876  4.528642850812465  1.3090607934204415 
C  6.609777342181284  -5.933316907087951  -2.8115581486852976 
H  7.19484518833961  -6.659911489193254  -2.989430261864678 
C  3.1522548408324207  3.759522550558459  0.10483494226788502 
H  3.715291056091863  4.5254430444606735  0.10298093932212227 
C  4.805240114902084  -4.4900521341927995  -3.5374216418322066 
H  4.170178704505117  -4.239484979549109  -4.198631015660491 
C  3.524799440215655  2.6063097314484738  -0.5930661791746806 
H  4.363229107985248  2.5610602863518173  -1.0381358253226207 
C  2.6487346941998187  1.5278151764926247  -0.62385961546895 
C  3.0225201382612896  0.2672584568890551  -1.3483533548967492 
H  3.885585956384991  0.3985385781124635  -1.8165580893432043 
H  2.3334455128330642  0.05806224159467763  -2.026668174208165 
C  3.146922718735412  -0.9140352636276334  -0.35644650753695617 
H  3.5091439297976716  -0.5442919909046543  0.4992046810514226 
C  4.068053538299674  -2.045300387340604  -0.7814738162120197 
C  5.018997449079372  -2.478762139834065  0.17260821258652034 
H  5.069094721310384  -2.0423436502992383  1.0155538011669303 
C  5.861294568572383  -3.5042962241103375  -0.09766647080317538 
H  6.489109372537924  -3.777970130934243  0.5600925843572249 
C  5.822708608026503  -4.171228946322092  -1.3294030163728927 
C  6.697675037778662  -5.264704929768077  -1.62894311457949 
H  7.35136665213379  -5.5276945092538545  -0.9911211061761256 
C  1.387715256305944  -2.278642515388447  1.650485711690801 
C  5.648903782226922  -5.548580237820806  -3.7724642768278436 
H  5.5857137632086795  -6.02737066341359  -4.591642261063962 
C  0.9894970149320421  -3.1880458427428935  4.256853137905669 
H  0.8272850856503826  -3.47658305548813  5.147528514718088 
C  4.872817603258959  -3.7632506444551264  -2.309468215360939 
C  3.9787292010601236  -2.670962659366792  -2.0174983767041574 
C  2.8852616176463526  -2.295671105594469  -2.9780240628211065 
C  3.1809912319248097  -1.6873370950917528  -4.252451085802177 
C  4.503054177941497  -1.412068356546802  -4.682555234061224 
H  5.2372877537313585  -1.6625916111579742  -4.133554637447832 
C  4.735824666732398  -0.7874437096549829  -5.884558133410208 
H  5.62734135046683  -0.6072575873707008  -6.156076998241582 
C  3.661340276160947  -0.41014702077071513  -6.721157842336328 
H  3.828354696671371  0.024161582599357044  -7.549512935506105 
C  2.360626024280031  -0.6786133626870875  -6.324649685784738 
H  1.6370909017834845  -0.4323167238919425  -6.888264618369159 
C  2.1036220332742737  -1.3060709715762475  -5.107152216524847 
C  0.7683403739529127  -1.5311876311978545  -4.673973235224009 
H  0.04346927515199095  -1.299372147864982  -5.242079415273384 
C  0.5143286175050517  -2.0872458995088494  -3.429040739130623 
H  -0.3815462908041847  -2.2268251676923208  -3.1486850604734733 
C  1.5785969485162283  -2.4466843558142144  -2.576710306882941 
C  1.2374083878397673  -2.9435469518038575  -1.1846417608802626 
H  1.8371968353239452  -3.6910860533687586  -0.9360879238628599 
H  0.3030879771647168  -3.2720043748871266  -1.1621172982336718 
C  1.800959411282889  -2.092059809432007  4.023962682764158 
H  2.232496816633252  -1.6621831361303574  4.752630016692728 
C  0.41225656618984585  -3.8662924512061108  3.1935080657658292 
H  -0.12077442913940861  -4.636844025140958  3.350160291810799 
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
