%chk=INAVOVIB_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.532480624673605  -1.4641048920757003  1.9834212081897683e-16 
N  1.532480624673605  1.4641048920757003  4.2743467955382577e-17 
C  1.6906836014151225  -1.936550000854311  -1.7390419567512936 
C  0.5426768193121282  -1.888425720018046  -2.538710001529602 
H  -0.28171281764425116  -1.6350465243006818  -2.1901511117325443 
C  0.6935749115222346  -2.2457929536308057  -3.9054871180000768 
H  -0.04839185941733515  -2.2176239609423587  -4.467549195241042 
C  1.8744432247225347  -2.6051612426432533  -4.381582024223415 
H  1.9625085385855843  -2.823610427350618  -5.281970131760863 
C  2.927702556816639  -2.6549456743424487  -3.587916030421795 
H  3.7534348876910064  -2.857127015038832  -3.966013654353761 
C  2.8772804727301313  -2.427633890329372  -2.239140075607502 
H  3.6078654632353264  -2.596957872061674  -1.687085702507488 
C  0.317427757608715  -4.510366006384849  1.3251759821568745 
H  -0.2957033248370262  -3.8443901632076347  1.5375394254405936 
C  -0.07306199022667803  -5.857191483625346  1.3256010602311687 
H  -0.9387170484352962  -6.080170873126026  1.5832994344125981 
C  0.7662644810307123  -6.819371910468123  0.9697829875489518 
H  0.4816803217097281  -7.702496143573376  1.0198130301816841 
C  2.015480513165574  -6.560218500727053  0.5372462893638724 
H  2.562397368390891  -7.251782411537103  0.24099639360766226 
C  2.479471233527471  -5.239313055039897  0.5385558587845102 
H  3.3408755599062534  -5.056942029954766  0.23748431295971212 
C  1.655424245504929  -4.168346597317736  0.995507144990211 
C  2.177705995857292  -2.830042459680384  0.9992190759723819 
C  3.2350641798806876  -2.3860665816417277  1.7339181643062287 
C  3.975087607965422  -3.19717993590173  2.761991368133097 
H  4.68505779126926  -3.697061521545771  2.329018436931252 
H  3.367481833449933  -3.829326587421931  3.1737950273534206 
C  4.567337610419294  -2.2921178483614297  3.833728832389465 
H  3.8600656540478804  -1.9675107573174184  4.404076801768285 
H  5.183769734631669  -2.8083111444090485  4.378740343602275 
C  5.275597498763311  -1.155390873396887  3.2470266123912097 
H  5.643717867017024  -0.6356261733956976  3.9752220264865077 
H  6.033290904135611  -1.5145304697782174  2.7467855867182047 
C  4.578710222934379  -0.2617562492101719  2.378815731762277 
H  5.222574735881723  0.2126912536700274  1.830200733810312 
H  4.1119416003175395  0.39473688881719704  2.91883931704664 
C  3.589889864055856  -0.9441016401771101  1.482699352437836 
C  2.881349146281872  -0.3959135463157976  0.4821182637673263 
C  2.764565500552404  1.0510096703648188  0.11937422850174793 
C  3.868919213094219  1.819360251725922  -0.17522774317735068 
H  4.725626972008939  1.4571200851432553  -0.16724917629509845 
C  3.6534157089904977  3.165271136024855  -0.49095775639226735 
H  4.376642842804118  3.737565351101268  -0.6108170251444611 
C  2.3770249117916977  3.6390660592946906  -0.6210912398236024 
H  2.2187138610506913  4.51247654164549  -0.8973759309178243 
C  1.3345838545731294  2.7883215451949477  -0.3276905002978613 
H  0.465823809492923  3.1196130582724964  -0.35309930654752975 
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
