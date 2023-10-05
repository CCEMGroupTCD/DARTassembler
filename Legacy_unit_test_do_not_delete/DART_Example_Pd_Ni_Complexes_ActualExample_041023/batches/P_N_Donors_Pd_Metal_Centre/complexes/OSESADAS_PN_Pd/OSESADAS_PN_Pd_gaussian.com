%chk=OSESADAS_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.4553071651716705  -1.5408377769901667  3.059958487531098e-15 
N  1.4553071651716698  1.5408377769901667  -4.1074281008466443e-16 
C  2.4351220030499294  1.4301936666588246  0.8197150628862548 
H  3.0438054891729474  2.1338155541736703  0.8529721392939696 
C  2.677834688388689  0.2908638837213948  1.7090063695385982 
C  2.4271236963886484  -0.9954076138995054  1.430741670199188 
C  2.996770345963208  -1.9098295641143144  2.4938134700292394 
H  3.5329929964222018  -2.6134063625217676  2.0965922245281194 
H  2.289208089501047  -2.3131980038562783  3.021663652420352 
C  3.8526704777604763  -0.9875952827250214  3.345135391877732 
H  4.793393445649992  -1.154743375045998  3.17651385434886 
H  3.6768726281344373  -1.1424932957178093  4.28658014949114 
C  3.4980308266941753  0.4434210402957696  2.9740565577472626 
H  2.979872376397915  0.8666143768933066  3.6754687120802556 
H  4.296885166898316  0.9677141632084448  2.8100511175829865 
C  1.497517993887715  2.721493711162274  -0.9548870409861626 
C  2.9312484147604128  3.0631983622699823  -1.384367135607097 
H  3.3689239281817245  2.2687493240388705  -1.6985214935349875 
H  3.413801421933124  3.418744311480894  -0.6342746991245641 
H  2.9074183303809438  3.714634431653427  -2.0879352308080867 
C  0.8660301416713242  3.911616268114422  -0.24679150227911523 
H  -0.03324737988425053  3.692587302158486  0.007959731863461375 
H  0.8546767966691285  4.666360915712337  -0.8386745553221706 
H  1.3771483738764119  4.126461419202714  0.5364328421851079 
C  0.7381377097863301  2.3350003055013095  -2.2079972442399654 
H  -0.16866563376266686  2.117979844299226  -1.979365694474824 
H  1.1566596896368737  1.5714604036576991  -2.6140314554520767 
H  0.7456619352082509  3.069325483593838  -2.826420504516603 
C  1.0970734689506376  -3.2779529230300226  0.4799475682725762 
C  1.6923748599751376  -4.411322237489104  -0.1325422201213732 
C  1.4563157472135964  -5.673615837792241  0.4141985428696304 
H  1.844131923592772  -6.414776557015241  0.008290102923514212 
C  0.6690987627160631  -5.867431039794459  1.5351682957137922 
C  0.07500394878503935  -4.756344119700778  2.1095785742230038 
H  -0.4795987987432748  -4.8747382314782985  2.8470322180071226 
C  0.27870855018071183  -3.4617547347646833  1.6229109157731847 
C  2.589690949322278  -4.363726954523459  -1.336012957130882 
H  2.936842650274566  -5.240312822400075  -1.5117516346180748 
H  3.3158668971894896  -3.758920232546298  -1.1705072720611247 
H  2.088619246035338  -4.060351802394786  -2.097690733645375 
C  0.484527771421443  -7.236457112169541  2.1509997437008215 
H  0.9658029393629779  -7.886754077212341  1.6347869290797576 
H  -0.4490042773923799  -7.460771842221883  2.1585426687997256 
H  0.8193258725419582  -7.230482729182915  3.0504332255292796 
C  -0.3595538916637009  -2.358626319572117  2.4438189282359777 
H  -0.0706559227548722  -2.43041193623584  3.356807125750373 
H  -1.3147741663726864  -2.441705136766003  2.402712518534081 
H  -0.0974668966681087  -1.50566479459752  2.0910358447877995 
C  2.482306949078147  -1.3189600880436936  -1.5141597653497136 
C  3.848498587725599  -0.9760912553285601  -1.5156089176583971 
C  4.399774655627819  -0.43640709946345624  -2.6865892756774645 
H  5.293061450673794  -0.1756483804081389  -2.6767159665475884 
C  3.677223548900252  -0.2743151299325248  -3.852102459464078 
C  2.3767143851264825  -0.7360796359857416  -3.8671821109824585 
H  1.8960272615529905  -0.7012808656020494  -4.662374751045496 
C  1.760069066752132  -1.2532359530878265  -2.7303111002655713 
C  4.791400329910392  -1.2430910694190636  -0.3674675962782245 
H  4.454723711912967  -1.971372587277488  0.15967205849995972 
H  4.860763870040759  -0.4563265682474811  0.1793725714095718 
H  5.658839431597182  -1.4704883108628892  -0.7114334949492489 
C  4.297222602144106  0.3553548356463299  -5.081318591161275 
H  5.201878894880514  0.6124113641886785  -4.886732509577326 
H  3.7913764373430627  1.1317081207641786  -5.332925097918947 
H  4.29113039688079  -0.277993012240857  -5.801744727980125 
C  0.3412371885756027  -1.763414650507586  -2.863811954829232 
H  0.2285446236464621  -2.176622500488839  -3.7234822105224006 
H  -0.2720564829609067  -1.029350859438412  -2.7797910461550677 
H  0.16618845192124465  -2.4085143941948433  -2.1736474233656256 
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

