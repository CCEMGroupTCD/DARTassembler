%chk=EDALACUG_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5337779141714096  1.4627458118210412  0.0 
N  1.5337779141714096  -1.4627458118210417  0.0 
H  1.5642414460730854  -1.779476508471077  -0.8146934141280858 
C  2.882755333495579  -0.8527383362985864  0.286828419024796 
H  2.906252797777842  -0.6308912611762129  1.2414544214110528 
C  2.9753634135609346  0.45194602825557695  -0.48501613941590366 
H  2.96087386875608  0.2804883436919354  -1.4401488032356782 
H  3.799203918582628  0.9164270095097014  -0.26719316005226934 
C  1.4024269406600205  2.8121323655536488  -1.1887483082483323 
C  1.4869940378024575  2.5391850518279977  -2.5334040225525327 
H  1.5780605870382736  1.6563330350498482  -2.8141140332059416 
C  1.4405281048815788  3.539051151720776  -3.472824189124995 
H  1.5234426347301409  3.3388336309230944  -4.377889018353052 
C  1.2716554424382496  4.8178481886928095  -3.0697700338172256 
H  1.2561011909093895  5.494967290313709  -3.7051106334010133 
C  1.1257006260225206  5.1303283587306625  -1.7685735094988446 
H  0.9772914580410282  6.011053280032036  -1.511134930692692 
C  1.200951838578392  4.112083713650645  -0.8056429887577288 
H  1.1127336257574427  4.322195304821347  0.09701594804089245 
C  1.924103031551341  2.2327687719946585  1.5814799662176726 
C  3.1945884294144005  2.640132837018487  1.8809642435078853 
H  3.887036043436299  2.4811796011787353  1.279192899348424 
C  3.4551286441501494  3.2873318939902676  3.0732620700429853 
H  4.33175841612103  3.519552651181919  3.2785131914129018 
C  2.5247217820696286  3.5739541537238626  3.904489524605654 
H  2.726151276191596  4.010519314149787  4.700056366297837 
C  1.2321764304920304  3.2307955024096655  3.612579103168846 
H  0.5479734794607735  3.4715384784989958  4.194539240346465 
C  0.9325750148661682  2.520562151129969  2.444755555502076 
H  0.060579912313796225  2.250119975323967  2.2678364478002058 
C  4.022956785748658  -1.8010230476888265  0.01671547149438768 
C  4.774770860997244  -2.2965371104484036  1.0673756971932944 
H  4.568519074229512  -2.048914087314215  1.9401552448850425 
C  5.83263289988453  -3.1577003691774  0.8273422026628015 
H  6.314451331258449  -3.5037797124980106  1.5434899369433488 
C  6.167483869696403  -3.4942404747936626  -0.4220922435279044 
H  6.904375017385144  -4.038811365759546  -0.5740856353928643 
C  5.431894706007707  -3.0388600220745854  -1.4718948037958606 
H  5.6502498266630505  -3.2946327523579324  -2.3384629440635463 
C  4.354339629775562  -2.190935905678795  -1.2495567716898288 
H  3.854944427958996  -1.8861176683632404  -1.9708995297476468 
C  1.2715574023940586  -2.6177555861417607  0.8495419203469845 
C  1.320037513863951  -3.8824577135038103  0.3239552501984895 
H  1.5165165962692035  -4.012533443998041  -0.5757363439122929 
C  1.0718704818852243  -4.967055124194669  1.155711537197016 
H  1.0844563237050455  -5.828427373438305  0.8046281071997465 
C  0.8116617730952023  -4.783337462637709  2.4688999123747664 
H  0.6564562866900698  -5.518911190694046  3.016917957317676 
C  0.7739137622093941  -3.496871287316865  3.010416373278334 
H  0.5927003754885658  -3.371304597249683  3.9139455564366443 
C  1.0077338466354488  -2.4112183890713514  2.187461281763528 
H  0.9872090531226878  -1.547182594344615  2.532862408400661 
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


