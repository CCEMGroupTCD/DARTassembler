%chk=NAWOVERE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
C  2.98029188410496  0.4319071734416376  0.48085811748977986 
H  3.809989544939482  0.8219807256699415  0.10786789582473934 
H  3.0631215979907473  0.4086186815404202  1.4674147983969603 
C  2.7880195486263553  -0.9591401262521897  -0.04756661605712819 
C  3.8537323105842716  -1.6797348986287326  -0.5568395026549581 
H  4.723425992528965  -1.298652772519552  -0.5828487845672293 
C  3.6376367673558825  -2.9596498451624815  -1.02563743976255 
H  4.356518176522107  -3.4727295005586787  -1.3731855325119144 
C  2.362116758097038  -3.4825471946376916  -0.9832716765550618 
H  2.189174792455087  -4.359837739720291  -1.3038215104788562 
C  1.3486353873594106  -2.711763116284106  -0.46871150930110017 
H  0.4723191150356958  -3.0787400216737937  -0.44212426101769614 
C  1.8844760489799393  1.8535553390081485  -1.7793262685351225 
H  1.8233566053089127  1.0126788562477076  -2.298042698143612 
H  1.1710547307216677  2.4560702739521107  -2.108814820421416 
C  3.2077550231874654  2.4909420275579293  -2.0720168657189943 
N  4.235899495492735  1.6592390754636728  -2.286694748943718 
C  5.418726117284569  2.215307695894803  -2.57481604658139 
H  6.1590633002395485  1.6386430494810744  -2.725790187118294 
C  5.628133186284632  3.577055795028477  -2.665980973821113 
H  6.48365505378556  3.9233023524526356  -2.88930404853502 
C  4.566951455415696  4.423428715956293  -2.426241734666937 
H  4.67891299203132  5.366111483642802  -2.4710248098803995 
C  3.3384047488105564  3.874585167389982  -2.119168862392593 
H  2.5924874725392684  4.434343729728948  -1.943350066900604 
C  1.8317341816594714  3.025986282151485  0.8783283326107298 
C  0.7326563990800358  3.8202834314837264  1.1776389010893147 
H  -0.13321292205508284  3.5520576750635415  0.8956080653464087 
C  0.8909167849716484  5.004625899913725  1.8860396053030943 
H  0.1325636164976467  5.542227139020646  2.0846384713023127 
C  2.1446731039251916  5.402624637241924  2.3033208027817906 
H  2.249082067244535  6.206164489265223  2.7985067522886062 
C  3.2453514431774555  4.6291717773342285  1.9994710915732667 
H  4.109204028808398  4.908586013868957  2.278262993752075 
C  3.097323644443011  3.445923228993782  1.2892990887091216 
H  3.8600369451376757  2.919474948584825  1.0825343939522338 
N  1.5345508243456778  -1.4619349395578443  1.7903539442911932e-16 
P  1.5345508243456778  1.4619349395578447  -1.7903539442911937e-16 
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


