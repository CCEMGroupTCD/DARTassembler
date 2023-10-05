%chk=ZIWIVEDU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.572948640293128  1.4205395365846032  -1.5618289360288314e-16 
N  1.5729486402931279  -1.4205395365846036  0.0 
N  2.9061722743841463  -1.3204724343751368  -0.3592137616469815 
C  1.506104344129938  -2.329101788292923  0.9793371813840086 
C  3.657128622176427  -2.1575725004944006  0.41409751656742844 
C  2.7890866229438505  -2.8004859544879466  1.256138933744996 
H  3.014731594730394  -3.436237668355302  1.896129023529336 
C  3.385369518110556  -0.30308226832964236  -1.237022318018435 
C  2.9831951055347163  1.030569081181508  -1.0928163895121492 
C  4.317283023656062  -0.6735353401710946  -2.2012137044895748 
H  4.554319528828783  -1.5679490275611534  -2.298944718117939 
C  4.893210396322497  0.29132774156712227  -3.019910269104462 
H  5.506630609375479  0.04571719640844396  -3.6755820027676647 
C  3.6073175614711817  1.977777490625071  -1.9142799913078798 
H  3.3806148615460545  2.875178338802082  -1.8250630448355725 
C  4.545664248976107  1.6132177640826368  -2.8491430186182587 
H  4.951493389494225  2.265855389522487  -3.37298281146199 
C  1.3512985958920614  3.24049399995411  -0.12531164781622586 
H  2.2467726357822113  3.6272731679039194  -0.025494970237344705 
C  0.5057270050125888  3.831050481061976  1.022336690011659 
H  -0.41903824362900144  3.551924265431732  0.9261195883467344 
H  0.8336976834341158  3.502687128274561  1.8746098578645158 
C  0.9008494045587376  5.2766207606032225  -1.5041680727552915 
H  1.8300816911595201  5.553354617617266  -1.455667867386368 
H  0.5364945532788068  5.603229262549952  -2.3418288327320758 
C  0.5928368452169315  5.347877876222386  0.9849448589457767 
H  0.03991827149885996  5.718896754209531  1.6911972380672298 
H  1.509080162969878  5.621238226286019  1.1461265008476391 
C  0.8224454864170675  3.74521013881203  -1.4753288098460375 
H  1.3540924360861126  3.3746720112597806  -2.197665245959903 
H  -0.09671194052915322  3.4612577322709144  -1.5962422719133433 
C  0.13395974065938843  5.890030322645428  -0.34983085919957324 
H  0.2540371915902735  6.851920836643538  -0.36102661100077577 
H  -0.8120922309533325  5.706962457723311  -0.4611160955992689 
C  2.1905662539283566  1.0370159520810103  1.6994009160457373 
H  2.614517994295799  0.1575920621950715  1.6153651871127008 
C  1.0932624667352795  0.8399580541551059  2.7652431949391914 
H  0.6858729613662152  1.6942962114146387  2.977708236688479 
H  0.4025014218534042  0.25223731953387984  2.4192597315049036 
C  3.300857445275337  1.950166414699439  2.2136275616379963 
H  3.9833700864107913  2.0572096585401516  1.5338499887228814 
H  2.938850481219938  2.8252996989084296  2.4241208937057555 
C  3.902727796256271  1.3179717436808516  3.464261796154953 
H  4.3251221192420894  0.4783922815931605  3.226288787706748 
H  4.587567983578615  1.905285309064543  3.8209614930213815 
C  1.6984194658374574  0.23522193250838752  4.0182765981704796 
H  1.0171884259231183  0.1774997516832455  4.706716748145798 
H  2.007796040078215  -0.6631563647041756  3.824406485959983 
C  2.8574039421345274  1.0642903636705339  4.532162675188289 
H  2.52223547087338  1.915565607235696  4.855252717692788 
H  3.2708606727177676  0.6036089196514469  5.278174149008857 
C  5.15870177733385  -2.232486420353233  0.3289284240396477 
H  5.41140730623612  -2.804634797127111  -0.39797322722912254 
H  5.509090515471633  -2.5850196248876456  1.1505446372967347 
H  5.5158093147801655  -1.3532341905962242  0.18239534529504775 
C  0.24152654766208492  -2.688592307105094  1.6573096760802792 
H  -0.48190638947778974  -2.191171322470699  1.269654920522569 
H  0.30936554668094285  -2.481190744239509  2.5920963854065615 
H  0.07807633411991155  -3.628208421704171  1.549442652024051 
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
