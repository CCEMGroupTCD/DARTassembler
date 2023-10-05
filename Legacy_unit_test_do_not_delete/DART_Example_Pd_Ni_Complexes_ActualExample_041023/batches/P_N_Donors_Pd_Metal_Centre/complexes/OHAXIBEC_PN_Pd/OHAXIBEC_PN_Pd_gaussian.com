%chk=OHAXIBEC_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
C  1.4227629038847949  2.7323717085757018  -0.21058701026494775 
H  0.540841431836909  3.080810304905272  -0.2635549414005567 
C  2.4770856109501262  3.5928376651443745  -0.3501320816677712 
H  2.3246967793179163  4.525157359979606  -0.44777628587933677 
C  3.780894089865818  3.0847414658772685  -0.3472610891649996 
H  4.525648141578391  3.6578290193568144  -0.4825834764648951 
C  3.9578462636442477  1.7324357826847325  -0.1448545739530192 
H  4.833045503431441  1.3627089552913512  -0.1360197688295926 
C  2.8516954575571267  0.9138896719817049  0.046066218120073915 
C  1.7666489087551958  -2.200077354121198  -1.6238058749687727 
C  2.732680874813934  -1.717339521885326  -2.5213667174084353 
H  3.298418298779535  -0.9891580502051855  -2.290475159688093 
C  2.835259656743184  -2.3517499108559887  -3.780016137169161 
H  3.501075312267697  -2.0668313814933033  -4.3951606197781485 
C  1.9778614376828225  -3.3870932211618747  -4.133978428602858 
H  2.0538066706772984  -3.799207400888589  -4.9866093888978815 
C  1.0151611638097253  -3.8132900397135776  -3.241790071504921 
H  0.42697408916062907  -4.519880847881061  -3.483465711628146 
C  0.8975854970391388  -3.2131307300137824  -1.9851147373138722 
H  0.22406758635303103  -3.501002284821272  -1.3798816070923339 
C  1.6919034802736916  -2.6963860262396766  1.2587977680770888 
C  1.3538421969169612  -2.362324163109622  2.562103369345409 
H  1.0151394951166164  -1.497366976550104  2.756187282156834 
C  1.5082154981343234  -3.2863594068083435  3.5819860892078035 
H  1.2803702993398713  -3.055373144489816  4.4752500154952175 
C  1.9931117713736113  -4.534976924335289  3.293515312310888 
H  2.1060119127357524  -5.165761144344722  3.9961899709314204 
C  2.320266057433496  -4.896074709818602  2.0037402805973126 
H  2.633490505927054  -5.77352833000817  1.8225690652176845 
C  2.1901901373815704  -3.9717384729814844  0.9657355665942956 
H  2.4358033267205093  -4.203053805627422  0.0772518311744092 
C  6.0955222697197184  -1.8734408510458027  2.5841130120335807 
H  6.781597536138075  -2.180373740047693  3.1654054173323387 
C  5.878737228232847  -2.4782259485377347  1.3774138607787212 
H  6.41031517870663  -3.2184602768026886  1.1094828961102015 
C  4.8682811036029845  -1.9968994346421964  0.5435298237452819 
H  4.7071746248721364  -2.3791545306602426  -0.31200304767416986 
C  4.109966045730439  -0.9385339779566664  1.017908673504231 
C  5.301667726605753  -0.8121443981645263  2.936314677956604 
H  5.4809968902909425  -0.38177474919952137  3.763627458168391 
N  2.977730600609041  -0.4676786815393178  0.24517241406859874 
N  1.5820488646056419  1.410397599969597  -1.378943933085085e-15 
N  4.282462348779418  -0.3391346383532634  2.18702137428096 
P  1.5820488646056419  -1.410397599969597  0.0 
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


