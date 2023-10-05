%chk=PARAVAQU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.4438034622136078  1.5516222357584324  0.0 
N  1.4438034622136078  -1.5516222357584328  0.0 
N  2.816479315539609  -0.3219304615442291  1.45312849814658 
C  1.0043773858979423  2.572504099587696  1.4315558786014493 
C  -0.15669618495264426  2.309520561716093  2.1603921303623093 
H  -0.7249272453409525  1.6208762542453972  1.8968644282166054 
C  -0.46333165863935966  3.069906435809536  3.268143893509679 
H  -1.222438713390201  2.8764711149672215  3.7679112354580804 
C  0.35665085025119025  4.120017959381226  3.637320210701785 
H  0.13067430891635334  4.651753704065701  4.366242045051138 
C  1.5023120241010393  4.377202524543726  2.9346818583645575 
H  2.0571323512574757  5.078533863397663  3.191062430396442 
C  1.8334438023861934  3.5994092336335233  1.8457205954975955 
H  2.6241724638931294  3.7671633615400295  1.3848476016676288 
C  1.9091331432342131  2.6777035361762116  -1.3342157294870671 
C  1.6172060295088369  4.024896832860467  -1.2803884906485274 
H  1.1025544327665426  4.360245023939763  -0.5825184062413635 
C  2.0870320352852962  4.878589596777665  -2.2586628791160885 
H  1.9065076401463383  5.789454043953606  -2.2020429035722575 
C  2.8130260956644864  4.394676197822678  -3.3067740834288735 
H  3.1460179965924997  4.979380840545552  -3.949329812291893 
C  3.0534197286138727  3.0539482119394736  -3.416481462884247 
H  3.5122464700953735  2.7199348625196182  -4.153005185320137 
C  2.6041510792213627  2.184051913153682  -2.4146191419878846 
H  2.774901205865054  1.2733203720390538  -2.4816801908176007 
C  3.014154585827488  0.735591557787969  0.4503714257984549 
H  3.413050964361006  0.3499803753679531  -0.3460103389259144 
H  3.62961239967773  1.3980282458633928  0.8018096190307114 
C  2.3678692959077314  -1.535728808588034  1.004536074002324 
C  2.843373616188046  -2.7562646847539893  1.5258853968047676 
H  3.4135961536858828  -2.754381652316928  2.2596425377706186 
C  2.478141489870737  -3.9200334735644686  0.9651791456390337 
H  2.768842202061931  -4.722432067312022  1.3332715289933812 
C  1.6666428963051607  -3.930376870707236  -0.16026641358481683 
H  1.475674432251801  -4.724244467221736  -0.6062658473752635 
C  1.1608147636168895  -2.7401086711172384  -0.5879854840339787 
H  0.5905925199126113  -2.7425750607524257  -1.322554840668898 
C  3.7018229684745347  -0.24587679507797144  2.639874446538618 
H  3.4741877633030764  -1.0063441709681058  3.2167621596798632 
C  5.172599444319319  -0.3974391209974959  2.262569583777121 
H  5.417397587329101  0.2983047882616574  1.6349095707861927 
H  5.307399044038372  -1.2552532553440585  1.8298663496421697 
C  6.078795193299343  -0.30473931837466584  3.509392534770695 
H  5.928368044992471  -1.0789857214708327  4.076628597045385 
H  7.007883048441394  -0.3116671161029241  3.2341138522575563 
C  5.7904485084363415  0.9655105656356637  4.300815178811332 
H  6.342619868814259  0.9803202657485048  5.097552922087632 
H  6.02282315888747  1.7377260394712297  3.761538233451507 
C  4.342134909971485  1.0541919223659892  4.696450967726223 
H  4.184932471098091  1.8777649580897098  5.182549145094821 
H  4.1171008013221355  0.31095221271056483  5.279926686103367 
C  3.464247882640718  1.014244753801663  3.4541030450539587 
H  2.5320183690851708  1.0541470048606887  3.719117449221792 
H  3.6530988882426323  1.7903869920645392  2.9047259581023397 
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


