%chk=ELEMEVIS_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.4790401507396624  -1.5180712211553185  -5.599589096171367e-16 
N  1.4790401507396613  1.5180712211553185  -1.8591010618655744e-16 
C  1.272734120978004  2.761924529014585  0.786416741054989 
H  0.3574734535242905  3.104860349643475  0.6310072704884004 
H  1.9159031032138243  3.453442038450722  0.4897639230454676 
C  1.4659274919561323  2.490623288822073  2.248991203597048 
C  0.3985248781580377  2.2821719094729724  3.096183420572934 
H  -0.48274068987229657  2.297279293553771  2.7422685359602053 
C  0.5764502719367588  2.0526879727872664  4.442998080075848 
H  -0.176023970142984  1.8973086422835983  5.001613707013325 
C  1.8368450073624574  2.0488249195139576  4.976851091428833 
H  1.9594515533880863  1.9065153776496515  5.906860316741133 
C  2.9359750416467287  2.2536270497877844  4.151751039817152 
H  3.814206410237073  2.248418782949037  4.514907737649532 
C  2.7443365934882324  2.4665920974710924  2.786452375182689 
H  3.4947141130759682  2.59734783029431  2.219740935002367 
C  2.5957829195909072  1.430453039513766  -0.6210141805662067 
H  3.1402038449473233  2.2089824052839373  -0.5954804251126966 
C  3.150292841938739  0.2932992512440527  -1.3663384441197837 
C  4.168244943590997  0.6122654323112964  -2.2684291624387267 
H  4.415927677309879  1.5220537978063653  -2.3847627063001156 
C  4.818403815909099  -0.354191252694189  -2.9885584402453564 
H  5.48828834695598  -0.10785471736655072  -3.614227974765647 
C  4.500846477104787  -1.6778325166862524  -2.807259844472962 
H  4.964476058826792  -2.349956710058078  -3.294033282158891 
C  3.4921256174164617  -2.0339439049178125  -1.9034493132212922 
H  3.277088510128576  -2.9510292139591585  -1.7795429436118353 
C  2.800222999959526  -1.064578749651859  -1.1837896549712823 
C  2.3813811244244962  -1.7984033452806456  1.553556937287126 
C  1.9102013732104626  -1.2026930173345265  2.7159727997151464 
H  1.1404588031122607  -0.6461704687040083  2.6800010240509042 
C  2.551263187339036  -1.4114500773351888  3.928476206046896 
H  2.225035623227975  -0.9898621102614621  4.715147784555232 
C  3.6501663995214626  -2.2228413630896116  3.9948213012460525 
H  4.073520123551132  -2.3784338489416523  4.830804209382678 
C  4.150872019269759  -2.8202871102725133  2.843951494740742 
H  4.91826565824531  -3.379298017121907  2.8933508761467235 
C  3.531224422176738  -2.602033605178168  1.6256146640855622 
H  3.8844984271992242  -2.9967237367413033  0.8373286262654869 
C  0.9840244018857636  -3.1807300607577838  -0.5594369414674059 
C  1.081916108010653  -4.324272564434812  0.22972143056956024 
H  1.4784838100604345  -4.262677404275701  1.091496885044767 
C  0.6177019204574695  -5.5423019521421875  -0.21368752147757974 
H  0.6874881227214509  -6.309972716701637  0.3404524547146175 
C  0.04705420163583862  -5.64190291808694  -1.4674758284044702 
H  -0.26926484047840327  -6.481981556841083  -1.7783918850900668 
C  -0.06056970457726374  -4.527393996224639  -2.2674191625855795 
H  -0.44966334788691764  -4.60053502442719  -3.13086834854241 
C  0.39651746878400296  -3.298768205149737  -1.8188196513796588 
H  0.3088383339267762  -2.5322037490273326  -2.3723200614251057 
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