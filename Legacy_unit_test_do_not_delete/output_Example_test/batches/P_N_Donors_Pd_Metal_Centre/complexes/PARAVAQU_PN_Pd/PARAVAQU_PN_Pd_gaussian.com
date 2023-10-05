%chk=PARAVAQU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.4438034622136078  1.5516222357584328  0.0 
N  1.4438034622136078  -1.5516222357584328  0.0 
N  2.816479315539606  -0.3219304615442289  1.4531284981465824 
C  1.0043773858979377  2.572504099587697  1.4315558786014513 
C  -0.15669618495264803  2.3095205617160937  2.16039213036231 
H  -0.7249272453409525  1.6208762542453976  1.8968644282165998 
C  -0.4633316586393623  3.069906435809537  3.2681438935096794 
H  -1.222438713390197  2.8764711149672224  3.7679112354580733 
C  0.3566508502511814  4.120017959381226  3.6373202107017857 
H  0.13067430891634468  4.651753704065701  4.3662420450511386 
C  1.50231202410104  4.377202524543726  2.934681858364563 
H  2.0571323512574655  5.078533863397663  3.1910624303964363 
C  1.833443802386191  3.599409233633524  1.8457205954975917 
H  2.6241724638931267  3.7671633615400304  1.3848476016676274 
C  1.9091331432342162  2.6777035361762125  -1.3342157294870691 
C  1.6172060295088415  4.024896832860467  -1.280388490648531 
H  1.1025544327665422  4.360245023939763  -0.5825184062413672 
C  2.087032035285296  4.878589596777665  -2.2586628791160903 
H  1.906507640146341  5.789454043953606  -2.2020429035722575 
C  2.8130260956644952  4.394676197822678  -3.3067740834288775 
H  3.1460179965925015  4.979380840545552  -3.9493298122918823 
C  3.053419728613874  3.0539482119394745  -3.4164814628842435 
H  3.512246470095364  2.719934862519619  -4.153005185320114 
C  2.6041510792213662  2.1840519131536826  -2.4146191419878824 
H  2.774901205865048  1.2733203720390542  -2.481680190817593 
C  3.0141545858274954  0.7355915577879695  0.4503714257984589 
H  3.4130509643610063  0.34998037536795334  -0.3460103389259138 
H  3.6296123996777236  1.3980282458633932  0.8018096190307104 
C  2.3678692959077314  -1.535728808588034  1.0045360740023286 
C  2.8433736161880416  -2.7562646847539893  1.5258853968047699 
H  3.4135961536858748  -2.7543816523169284  2.259642537770621 
C  2.478141489870737  -3.920033473564469  0.965179145639041 
H  2.7688422020619257  -4.722432067312022  1.33327152899338 
C  1.6666428963051592  -3.9303768707072364  -0.16026641358481875 
H  1.4756744322517998  -4.724244467221736  -0.6062658473752656 
C  1.1608147636168908  -2.740108671117239  -0.5879854840339774 
H  0.5905925199126127  -2.742575060752426  -1.322554840668901 
C  3.7018229684745236  -0.24587679507797122  2.6398744465386206 
H  3.4741877633030698  -1.0063441709681058  3.2167621596798663 
C  5.172599444319323  -0.3974391209974957  2.2625695837771342 
H  5.417397587329095  0.2983047882616576  1.634909570786202 
H  5.30739904403836  -1.2552532553440585  1.829866349642177 
C  6.078795193299332  -0.3047393183746656  3.5093925347706962 
H  5.928368044992472  -1.0789857214708327  4.076628597045399 
H  7.007883048441395  -0.31166711610292386  3.2341138522575656 
C  5.790448508436348  0.9655105656356642  4.300815178811348 
H  6.342619868814244  0.9803202657485053  5.097552922087631 
H  6.022823158887456  1.7377260394712302  3.7615382334515135 
C  4.3421349099714925  1.0541919223659897  4.696450967726232 
H  4.184932471098079  1.8777649580897102  5.182549145094829 
H  4.117100801322117  0.31095221271056506  5.279926686103347 
C  3.4642478826407173  1.0142447538016635  3.454103045053972 
H  2.532018369085163  1.0541470048606891  3.719117449221786 
H  3.6530988882426305  1.7903869920645397  2.9047259581023455 
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