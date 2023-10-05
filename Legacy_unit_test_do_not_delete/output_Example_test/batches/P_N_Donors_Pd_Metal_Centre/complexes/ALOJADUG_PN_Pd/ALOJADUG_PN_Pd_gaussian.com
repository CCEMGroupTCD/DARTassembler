%chk=ALOJADUG_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.6062497494163228  -1.3827732071818573  1.6934087821219814e-16 
N  1.6062497494163228  1.3827732071818573  -1.6934087821219814e-16 
C  2.6563439860833586  -0.47147761405407257  -1.1424922258223986 
C  3.546851351755946  -1.014438177724626  -2.0669790823982375 
H  3.639435145670876  -1.96718823330676  -2.1262491650079136 
C  4.26198146001538  -0.16644038682890985  -2.8972130162366563 
H  4.885750076742091  -0.538539001647045  -3.5230154000906313 
C  4.0939999640141895  1.1764822814404212  -2.8077093624029756 
H  4.597384065961452  1.7615550546168834  -3.379718135058038 
C  3.2095150741029066  1.7172427916356823  -1.9388854158338187 
H  3.062417726490138  2.665485730672868  -1.914893773718059 
C  2.5159587120856974  0.8860613664890479  -1.0123695869368141 
C  2.5874069637739194  -1.7873068661323608  1.4811662058966337 
C  3.8751346082433593  -2.341955116486502  1.371173651383477 
H  4.214126533552138  -2.623251762226936  0.5182383363269111 
C  4.615988880857367  -2.516829458673692  2.5368036875416835 
H  5.496169849611277  -2.8948537921544704  2.4833063461137 
C  4.131177823478338  -2.1183039397108954  3.7808785966120806 
H  4.664365550769579  -2.252936010272497  4.568514167826514 
C  2.9016210796720183  -1.5657112795180619  3.8306386776152426 
H  2.5510440564464942  -1.2949291888214225  4.682999851146946 
C  2.115138070190473  -1.3840085182156034  2.731413167389086 
H  1.249228659637937  -0.97480384020973  2.7831080821746697 
C  1.2402161501492854  -2.924350563416857  -0.916577050198808 
C  0.48861828204287394  -2.8222100863910944  -2.0496370078404405 
H  0.21285281754525864  -1.949953758501908  -2.3412328189346665 
C  0.10999663320722441  -3.9208412987329213  -2.7466476794239854 
H  -0.40030647151116483  -3.8371951307284613  -3.5556482170585157 
C  0.47606815790928403  -5.162195416736344  -2.3244490857287117 
H  0.21402738636890684  -5.951337410038278  -2.8041131853421657 
C  1.2137148243909903  -5.297441915604471  -1.1742137011345934 
H  1.4914615373763476  -6.170916740490309  -0.8881757600909427 
C  1.5762273245091711  -4.196773449441579  -0.43840554349246497 
H  2.0482762185799372  -4.289379669561347  0.3930454470065036 
P  1.9577783424207518  2.672054226470174  0.861670100214241 
C  1.4677983311357135  2.4755754159255985  2.5469230116635 
H  1.875572771267655  1.6856237779225771  2.9088253542801548 
H  0.5131177989100575  2.393848398499278  2.5952049545816154 
H  1.7486209457929227  3.2423680030924635  3.0529648749149945 
C  1.1757110402980095  4.1762252786696665  0.22845853945784667 
H  1.4267538859024487  4.3044363468905456  -0.688733927615723 
H  1.464875040327786  4.930433904193636  0.7480667286466143 
H  0.2213622655114338  4.0914795124143115  0.2900775394917237 
C  3.751945003032026  2.976844645326308  0.9476731497215959 
H  4.09894940109764  3.095951346370698  0.05961006553811277 
H  4.1847748044826645  2.2270529697173926  1.3620661228325979 
H  3.917331471960891  3.768659079763787  1.464836767266078 
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
