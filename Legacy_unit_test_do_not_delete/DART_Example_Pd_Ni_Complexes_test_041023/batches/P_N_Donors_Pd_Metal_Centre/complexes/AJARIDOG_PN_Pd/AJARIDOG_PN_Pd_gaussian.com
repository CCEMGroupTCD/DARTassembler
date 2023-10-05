%chk=AJARIDOG_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.466961163425944  -1.5297466930835317  1.8372312184652955e-16 
N  1.466961163425944  1.5297466930835317  -2.428510903503574e-16 
C  2.8927865163109505  -0.9934474164063555  -1.0193191104389032 
C  3.7810309831950146  -1.9902563011023633  -1.3987292769354518 
C  4.8093064907220855  -1.7772085615821285  -2.317464620386071 
C  4.890291634632323  -0.5646799590452404  -2.929083729387576 
C  3.9960177281051292  0.48843462017485717  -2.62134318422216 
C  4.0583685659352815  1.6856627908782755  -3.366206990764031 
C  3.2332929645966164  2.7236841980690496  -3.0925617637102034 
C  2.359348036612677  2.635835876116402  -1.9958070166581146 
C  2.2629530944218184  1.5069309712664458  -1.227978447949096 
C  3.029975519242692  0.3225016417623732  -1.5778121845046666 
C  2.383873262176958  1.369341133058829  1.16284299863619 
C  0.781069069165854  2.8197262785417663  0.2340337670775639 
C  2.0543414451423363  -1.7212425344044988  1.7184265440825874 
C  3.3231768665055617  -2.230106994461544  2.0241931458965 
C  3.766884563779679  -2.285426676068058  3.3270380108772626 
C  2.9510220564836454  -1.8417341333295862  4.359888598671501 
C  1.681918008499162  -1.376799593372282  4.0886344577922795 
C  1.2340054335444122  -1.313879065133726  2.7668664288539078 
C  1.1459903790419488  -3.237268563737033  -0.5477444995724765 
C  0.9051487746796171  -3.4989043965252966  -1.9034947540844926 
C  0.5185009962040827  -4.761883165321041  -2.315260530917758 
C  0.35633626719610056  -5.778510971762117  -1.3945064302825352 
C  0.5737294718725641  -5.532940396778374  -0.04328828828521965 
C  0.9692510156583172  -4.2567236814052265  0.37844709534236903 
H  3.684364016807514  -2.8809658814273713  -1.0004261246422543 
H  5.455312571752918  -2.485437619960516  -2.5193630456887086 
H  5.586681613312658  -0.4150419398653564  -3.602902815500127 
H  4.708585656273936  1.764080963339361  -4.095482436435591 
H  3.2441029703611877  3.525304626826557  -3.655375722291544 
H  1.8040199190462298  3.4109188786558664  -1.76877659795939 
H  2.963365600333969  0.5824935847410551  1.0906472527529543 
H  2.9349701545294415  2.1766985704594553  1.2309270166826096 
H  1.822194191510233  1.2890466074649163  1.962105502907873 
H  0.17164095321523742  3.053516903236042  -0.49624435480063755 
H  0.2685339061119747  2.696627169485731  1.0604029860003006 
H  1.4324482117205326  3.5414862176166673  0.3585677503565751 
H  3.901237914122513  -2.5505902138184555  1.2999501056547327 
H  4.65848990529637  -2.6395900564823176  3.524362679767108 
H  3.2774559953344236  -1.855955690466047  5.28315330116968 
H  1.0937449551853256  -1.0948134959678228  4.8205442406385846 
H  0.33463473983693937  -0.976874473242218  2.574284679084655 
H  1.011512829924056  -2.7827745791162517  -2.56409331778988 
H  0.3602358580570548  -4.9361202244401134  -3.2666011025727615 
H  0.08409376137192548  -6.6708153618755865  -1.6918046196594678 
H  0.4553455306107117  -6.250832675915358  0.6117988265375581 
H  1.1253684302663876  -4.083876293536703  1.3304090343730983 
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


