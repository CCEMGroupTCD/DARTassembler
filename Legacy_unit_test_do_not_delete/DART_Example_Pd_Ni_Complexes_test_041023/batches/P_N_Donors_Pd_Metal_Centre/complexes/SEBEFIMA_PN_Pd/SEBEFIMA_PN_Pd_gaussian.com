%chk=SEBEFIMA_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5235621705398175  1.4733832877089377  -7.397292003860191e-17 
P  2.737373268386668  -1.0680949693321642  0.9715889757151474 
N  1.5235621705398175  -1.4733832877089381  8.326672684688674e-17 
C  2.377087452308677  0.6492318940222459  1.3994377634127677 
H  3.202812616212962  1.115817782252984  1.601271488541111 
H  1.812566083740042  0.6781343563721598  2.1878381579455883 
C  1.0760111586262862  3.0980511185231534  0.6777060655273915 
C  0.37476912257476935  3.155657815781123  1.8669159609223893 
H  0.13266862775979238  2.367897971019476  2.299937239142588 
C  0.03412403173822898  4.375121284647982  2.4181487522076903 
H  -0.42558705167687405  4.413825531177528  3.22584005529688 
C  0.3778718890484092  5.5233458295528255  1.7640121112275908 
H  0.14551338161581873  6.346255720084663  2.1321739453285393 
C  1.0476479453681231  5.484201062891543  0.594542004032071 
H  1.2704510820668033  6.27781902690295  0.16371717945063594 
C  1.4042955630919576  4.276777468609922  0.034593103913527046 
H  1.8653762084183485  4.255124576069239  -0.7735625560617297 
C  2.768981010282829  1.790051749942894  -1.2591516662577318 
C  3.978515674062256  2.3969877132186315  -0.9237790919488794 
H  4.157545581497072  2.6248151277143994  -0.03989334059406774 
C  4.904616794199702  2.655775642456583  -1.905400349339406 
H  5.714355538547364  3.0573498084836923  -1.685072831580124 
C  4.642963546967037  2.32999068615618  -3.1945918170667693 
H  5.281271569991239  2.500401683141507  -3.8490010008726374 
C  3.453482147719847  1.7559477727774313  -3.5474181545822736 
H  3.275089356170386  1.55562496304952  -4.4381880035477295 
C  2.5121327383957803  1.4725811983884887  -2.563692169317871 
H  1.7082684089931166  1.0665606265587257  -2.793176145152902 
C  4.386604657304696  -1.1340171097624598  0.2281580729348842 
C  5.389152221883636  -0.2180934113890125  0.46179338600796715 
H  5.225581236537815  0.4956938792455443  1.0358328415174134 
C  6.631831380723293  -0.32994494817158193  -0.14160910976501315 
H  7.295878013978054  0.29595760053587683  0.04135123671618815 
C  6.878674891818431  -1.3505988320217117  -0.994279145835192 
H  7.712681659874605  -1.4187344969484554  -1.3996170320503951 
C  5.91762077797277  -2.278973611626388  -1.2667206970977387 
H  6.088707998529552  -2.9645811473075954  -1.8702699411361259 
C  4.690509538952718  -2.1995372721644455  -0.6368337746714309 
H  4.054796733901089  -2.861762989762477  -0.7885558180440772 
C  2.8260933982959666  -1.9393861979683882  2.545190155636699 
C  1.7099939528667822  -1.9206702126371842  3.3646150697110855 
H  0.9230844125684559  -1.5238821467775077  3.0698685749417267 
C  1.7657853708671478  -2.488154941135865  4.610309554967556 
H  1.0130472745309067  -2.4774367898362417  5.156565190688455 
C  2.92386617140694  -3.0740275020648475  5.062898785696335 
H  2.9643949777961067  -3.427647282706017  5.9220454126419195 
C  4.002121648708378  -3.1347329415429073  4.250781389449146 
H  4.77343168005209  -3.564112258708731  4.546713488990203 
C  3.9779279123391778  -2.5698002339822725  2.9976348102139894 
H  4.733111357744757  -2.6089819835525203  2.4548811512802136 
C  1.4929676084509311  -2.724723092076406  -0.694076304721744 
C  1.6398687625216017  -3.952551264888811  -0.060560528295559954 
H  1.812873753855445  -3.9825362256037504  0.8516621172819936 
C  1.532243117009995  -5.125163299661194  -0.7773430101791498 
H  1.6429552843273332  -5.932293440148388  -0.3321995587783402 
C  1.2673382495956547  -5.14023278645748  -2.121694281903886 
C  1.1361531378629826  -3.926259711751788  -2.745183480165814 
H  0.9664933597840071  -3.902767060322423  -3.6592607061325446 
C  1.247863129252508  -2.738229194097519  -2.0538129738915845 
H  1.1568399464697605  -1.9327140407082855  -2.512507281526223 
C  1.1061790841147683  -6.452899474164656  -2.871078702855772 
H  0.9264611511263652  -6.271034909638026  -3.796655087997451 
H  0.3755265567626629  -6.94854038178995  -2.4936156649800365 
H  1.913542116766572  -6.96498428352422  -2.7963639684958244 
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

