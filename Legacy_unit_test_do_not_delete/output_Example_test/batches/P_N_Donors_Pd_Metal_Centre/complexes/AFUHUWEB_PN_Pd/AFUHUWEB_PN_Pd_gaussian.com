%chk=AFUHUWEB_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.4553071651716685  -1.5408377769901667  -1.97905822090339e-15 
N  1.4553071651716698  1.5408377769901667  -4.1074281008466443e-16 
C  2.4351220030499325  1.4301936666588246  0.8197150628862557 
H  3.0438054891729456  2.1338155541736703  0.8529721392939664 
C  2.677834688388687  0.2908638837213948  1.7090063695385957 
C  2.4271236963886516  -0.9954076138995054  1.4307416701991846 
C  2.9967703459632142  -1.9098295641143144  2.493813470029228 
H  3.5329929964222186  -2.6134063625217676  2.096592224528132 
H  2.2892080895010496  -2.3131980038562783  3.021663652420347 
C  3.8526704777604763  -0.9875952827250214  3.3451353918777222 
H  4.79339344564999  -1.154743375045998  3.1765138543488547 
H  3.6768726281344435  -1.1424932957178093  4.286580149491123 
C  3.4980308266941718  0.4434210402957696  2.974056557747259 
H  2.9798723763979145  0.8666143768933066  3.6754687120802574 
H  4.296885166898313  0.9677141632084448  2.8100511175829874 
C  1.4975179938877135  2.721493711162274  -0.9548870409861626 
C  2.9312484147604168  3.0631983622699823  -1.3843671356071046 
H  3.3689239281817196  2.2687493240388705  -1.6985214935349817 
H  3.4138014219331234  3.418744311480894  -0.6342746991245666 
H  2.907418330380933  3.714634431653427  -2.0879352308080805 
C  0.8660301416713259  3.911616268114422  -0.24679150227911098 
H  -0.03324737988425208  3.692587302158486  0.007959731863466537 
H  0.8546767966691278  4.666360915712337  -0.8386745553221672 
H  1.3771483738764139  4.126461419202714  0.5364328421851103 
C  0.7381377097863281  2.3350003055013095  -2.2079972442399645 
H  -0.1686656337626704  2.117979844299226  -1.9793656944748195 
H  1.1566596896368675  1.5714604036576991  -2.614031455452081 
H  0.745661935208249  3.069325483593838  -2.826420504516601 
C  1.0970734689506396  -3.2779529230300226  0.4799475682725774 
C  1.6923748599751387  -4.411322237489104  -0.13254222012137373 
C  1.4563157472135968  -5.673615837792241  0.4141985428696307 
H  1.8441319235927736  -6.414776557015241  0.00829010292351324 
C  0.6690987627160669  -5.867431039794459  1.5351682957137953 
C  0.0750039487850358  -4.756344119700778  2.109578574223008 
H  -0.4795987987432724  -4.8747382314782985  2.847032218007134 
C  0.27870855018071405  -3.4617547347646833  1.6229109157731831 
C  2.589690949322276  -4.363726954523459  -1.3360129571308852 
H  2.9368426502745635  -5.240312822400075  -1.5117516346180733 
H  3.3158668971894905  -3.758920232546298  -1.1705072720611294 
H  2.0886192460353348  -4.060351802394786  -2.097690733645383 
C  0.48452777142144243  -7.236457112169541  2.1509997437008233 
H  0.9658029393629813  -7.886754077212341  1.6347869290797534 
H  -0.4490042773923715  -7.460771842221883  2.1585426687997273 
H  0.8193258725419634  -7.230482729182915  3.0504332255292814 
C  -0.3595538916637011  -2.358626319572117  2.443818928235984 
H  -0.07065592275486599  -2.43041193623584  3.3568071257503713 
H  -1.3147741663726762  -2.441705136766003  2.402712518534076 
H  -0.09746689666810204  -1.50566479459752  2.0910358447878012 
C  2.4823069490781426  -1.3189600880436936  -1.514159765349712 
C  3.8484985877255955  -0.9760912553285601  -1.5156089176584033 
C  4.399774655627815  -0.43640709946345624  -2.686589275677469 
H  5.293061450673776  -0.1756483804081389  -2.676715966547592 
C  3.6772235489002423  -0.2743151299325248  -3.8521024594640814 
C  2.376714385126471  -0.7360796359857416  -3.867182110982453 
H  1.896027261552971  -0.7012808656020494  -4.662374751045486 
C  1.7600690667521288  -1.2532359530878265  -2.730311100265567 
C  4.791400329910395  -1.2430910694190636  -0.36746759627822706 
H  4.454723711912982  -1.971372587277488  0.15967205849995406 
H  4.860763870040762  -0.4563265682474811  0.17937257140956203 
H  5.658839431597179  -1.4704883108628892  -0.7114334949492603 
C  4.297222602144095  0.3553548356463299  -5.081318591161281 
H  5.201878894880511  0.6124113641886785  -4.8867325095773415 
H  3.791376437343044  1.1317081207641786  -5.3329250979189435 
H  4.2911303968807815  -0.277993012240857  -5.801744727980127 
C  0.34123718857559737  -1.763414650507586  -2.8638119548292353 
H  0.22854462364645856  -2.176622500488839  -3.723482210522398 
H  -0.2720564829609049  -1.029350859438412  -2.779791046155066 
H  0.1661884519212371  -2.4085143941948433  -2.173647423365631 
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