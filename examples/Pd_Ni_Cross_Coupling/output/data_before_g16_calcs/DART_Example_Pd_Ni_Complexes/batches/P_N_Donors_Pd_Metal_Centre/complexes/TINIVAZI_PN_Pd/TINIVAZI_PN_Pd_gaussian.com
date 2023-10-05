%chk=TINIVAZI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.581344186127739  -1.4111876434408006  -7.295372825141504e-16 
N  1.5813441861277384  1.4111876434408006  -3.948652479786386e-16 
N  2.6049606163874675  1.380535514929717  -0.9383796197275954 
C  4.8551834400248755  2.519368325342869  -1.352904304044249 
H  4.67629650854987  3.054336930339493  -2.1537085231986146 
H  5.230279434531008  1.6532615423923753  -1.6156360977666704 
H  5.493993369038713  2.9922078184474525  -0.7789358046002942 
C  3.576601636575716  2.3016823383461396  -0.6013833176746168 
C  3.1490946961237007  2.925771901554496  0.5412495205371034 
H  3.600345328891051  3.623632653508645  1.000179185419473 
C  1.9173253017053349  2.3445539403483275  0.9132994775934011 
C  1.1418559023208412  2.710117054791543  2.1802771687848623 
C  0.5763081918740507  4.133957650575708  2.016287275787153 
H  1.3043574746750006  4.750915203348187  1.7909234146868818 
H  0.15555835142657082  4.417264638032879  2.854643812483418 
H  -0.08930357688076662  4.139692943849863  1.2974298353261353 
C  2.148501579718604  2.6930438185624803  3.366118146648709 
H  2.513844769408  1.788444496896473  3.4714676174066303 
H  1.6865574109247286  2.9563542909345326  4.188877674127197 
H  2.877460573501307  3.3209223073788894  3.1840712058204654 
C  -0.008980996968735111  1.7565634070859752  2.525574970951326 
H  -0.7152192508721618  1.839034179517226  1.8508470551081841 
H  -0.3701851330265167  1.9871852606096332  3.407113746486528 
H  0.3240283309498475  0.834566172311731  2.5404500033165878 
C  2.70113914009949  0.32808326529215504  -1.91813812252786 
C  3.146100312976091  0.6564924085794424  -3.197190071917568 
H  3.301538655846768  1.5653443893640961  -3.4284100723272077 
C  3.363575198406507  -0.35132152857516147  -4.141717002335126 
H  3.6790638484939366  -0.1275284905139919  -5.00906174220582 
C  3.1211912928752223  -1.677920538070056  -3.8157006985302635 
H  3.2939916313677626  -2.363270841763082  -4.450688859242547 
C  2.6226697125750382  -2.00283154627733  -2.5547679601493543 
H  2.427917431621467  -2.909740493467023  -2.3496877549352226 
C  2.4029561399856028  -1.013915737322727  -1.5821152580919167 
C  2.8253762338362804  -1.1735200507865493  1.3075552797782934 
C  4.168274118421393  -0.9249064616365129  1.013291699179584 
H  4.454145956576979  -0.8794584907035881  0.10910172241162852 
C  5.086626826128743  -0.7442854179470142  2.037549963887847 
H  5.99802734534063  -0.5752095185731769  1.829290741207814 
C  4.682089821747054  -0.8093214501335093  3.367260288207892 
H  5.315239177731872  -0.6871493734423963  4.0651724078139715 
C  3.3483057838400887  -1.052577569612421  3.669036878686148 
H  3.072642564478997  -1.0940830173622238  4.577044002721287 
C  2.4072079318079127  -1.236219972665302  2.652891554088593 
H  1.495581505055823  -1.4012796697811742  2.865374410891774 
C  1.1748210728529735  -3.189862401306405  -0.007614072834332934 
C  1.7582815871226736  -4.064073095918225  0.9195626145961033 
H  2.429135177977714  -3.7475457956147076  1.5134760912907195 
C  1.353057760547938  -5.412947424465202  0.9708191864974426 
H  1.744351679414008  -6.007350248142693  1.5996753165166049 
C  0.3639369341168488  -5.871242678824926  0.08493435168442293 
H  0.0867831988563792  -6.779351399133288  0.11932823385955976 
C  -0.21337852243944333  -5.0072664771716635  -0.8453055725276193 
H  -0.8790862907434172  -5.327663475474155  -1.4436242674540367 
C  0.18924666085826258  -3.6671494319366844  -0.8955052422990671 
H  -0.2030973304490662  -3.07895250514378  -1.5291825953711244 
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
-H 0
6-31g(d,p)
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
-Br 0
6-31+g(d)
****

-Pd 0
lanl2dz


