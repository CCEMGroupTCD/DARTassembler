%chk=BOKALUBU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5505144275368736  1.4449930830284272  3.4195460656297264e-15 
N  1.550514427536874  -1.4449930830284277  -7.216449660063518e-16 
N  3.645315135015532  -1.9884173036857862  -0.4138004499151323 
C  1.6756058280427144  -2.8127708621841485  0.06165142632035099 
C  2.9461526433925664  -3.1507003067850894  -0.20051891997396157 
C  5.052907129963877  -1.8649307540989672  -0.7535026545422423 
C  2.7642654195214638  -0.9686164705139035  -0.2755033802073009 
C  3.0686635483682743  0.4570228520512636  -0.43216515317854837 
C  1.9582826267054354  2.1279885254229303  1.6582283713093409 
C  2.082572863746561  1.0037718498694457  2.6638250963357013 
C  0.9269224453445563  3.1708059324390954  2.08350979804622 
C  1.5377303354957355  2.839959345272014  -1.1833248346049694 
C  1.4974286932941858  2.3124630762916007  -2.6011883237489855 
C  2.723090099267295  3.7919653279698124  -0.9959027992433869 
H  0.9640906409346944  -3.3933705618323264  0.2280500659875611 
H  3.289057191503475  -4.008309969111495  -0.25934415239381203 
H  5.434737934050472  -2.7071287163545827  -0.8211299376687569 
H  5.124641381619847  -1.424349551685366  -1.6366035441461255 
H  5.484042580512787  -1.3196848499810307  -0.13178898465664404 
H  3.3154153696527713  0.6592277784496929  -1.3598187234476438 
H  3.7737192344713417  0.7193887377103458  0.11995987028631518 
H  2.7909468845087853  2.588780246726209  1.5637658212951553 
H  2.311235338284136  1.3877076531033798  3.477899872835685 
H  1.244202887306403  0.5725229827856009  2.6981997562271176 
H  2.7482648722109344  0.4244027858698982  2.342438090458117 
H  1.1665825850046712  3.5546356659611416  2.8985119962779358 
H  0.8835270073766801  3.8899507243450353  1.4104589052913319 
H  0.06580023350101971  2.782796248112482  2.125505406607963 
H  0.7336626709373718  3.36048940789599  -1.0668736228622426 
H  0.7577726168176785  1.739266782324494  -2.7263945021659044 
H  1.3826635118103776  3.064496959888137  -3.235649355682154 
H  2.298457493714826  1.8820658534228436  -2.8236422548078317 
H  3.5267662746978274  3.354611061418236  -1.2035553469630793 
H  2.6324950504693647  4.549722436194295  -1.625994691008849 
H  2.7255259931150815  4.139829913554206  -0.13279087920208402 
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
